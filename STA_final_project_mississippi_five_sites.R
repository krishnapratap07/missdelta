################################################################################

# Description:
# This script performs a Bayesian analysis of spatio-temporal data.
# The model includes:
#   - Fixed effects (regression coefficients) with a Horseshoe prior for shrinkage.
#   - Spatial random effects modeled using a Gaussian Process with an exponential
#     covariance function.
#   - Temporal random effects modeled using an AR(1) process.
# Model Structure (Simplified):
#   Y_it = β₀ + X₁_it * β₁ + X₂_it * β₂ + s_i + u_t + ε_it
#   where:
#     Y_it: Centered response variable for site i at time t
#     X₁_it, X₂_it: Centered covariates
#     β₀, β₁, β₂: Regression coefficients (with Horseshoe prior)
#     s_i: Spatial random effect for site i (MVN with exponential covariance)
#     u_t: Temporal random effect for time t (AR(1) process)
#     ε_it: Observation error (Normal distribution)
#
################################################################################

# --- 0. Setup: Load Libraries and Set Seed ---

# Consolidate all library calls at the beginning
library(readr)     
library(dplyr)     
library(geosphere) 
library(horseshoe) 
library(ggplot2)  
library(viridis)   
library(tidyr)    
library(knitr)     
library(kableExtra)

# Set a seed for reproducibility
set.seed(123)

# --- Helper Function: Visualize Spatial Variance Effect ---

#' @param sigma_s2_samples A numeric vector of posterior draws for σ_s² (nsamp).
#' @param s_matrix A matrix (nsamp × n_sites) of posterior draws for s_vec.
#'
#' @return A plot displaying the relationship.
visualize_spatial_variance_effect <- function(sigma_s2_samples, s_matrix) {
  # Calculate the standard deviation of spatial effects for each posterior draw
  s_sd_per_draw <- apply(s_matrix, 1, sd)
  
  plot(sigma_s2_samples, s_sd_per_draw,
       pch = 20, col = rgb(0, 0, 1, 0.3), 
       xlab = expression(sigma[s]^2),    
       ylab = "SD of spatial effects (s_i)", 
       main = expression("Posterior Relationship: "~sigma[s]^2~"vs. Spread of Spatial Effects") 
  )
  # Add a LOWESS (Locally Weighted Scatterplot Smoothing) line to see the trend
  lines(lowess(sigma_s2_samples, s_sd_per_draw), col = "red", lwd = 2)
}

# --- 1. Data Loading and Preprocessing ---

site_files <- c(
  "site_01300_dataframe.csv", "site_01320_dataframe.csv",
  "site_01340_dataframe.csv", "site_01400_dataframe.csv",
  "site_01480_dataframe.csv"
)

# Specify columns to keep from the original datasets
columns_to_keep <- c("Date", "original_stage", "tidal value", "StormSurge", "Flow (CFS)")

# Load and preprocess data for each site
df_list <- lapply(seq_along(site_files), function(i) {
  read_csv(site_files[i]) %>%
    select(all_of(columns_to_keep)) %>%
    # Rename and transform columns
    transmute(
      siteID = as.factor(i), # Use site index as a factor for site ID
      time   = as.Date(Date),
      Y      = original_stage - `tidal value`, # Response: original stage adjusted for tidal value
      X1     = StormSurge,                    # Covariate 1: Storm Surge
      X2     = `Flow (CFS)`                   # Covariate 2: Flow in Cubic Feet per Second
    )
})

# Combine all site data into a single data frame
df <- bind_rows(df_list)

# Center the response variable and covariates by subtracting their respective means.
df <- df %>%
  mutate(
    Yc  = Y  - mean(Y, na.rm = TRUE),
    X1c = X1 - mean(X1, na.rm = TRUE),
    X2c = X2 - mean(X2, na.rm = TRUE)
  )

# --- 2. Spatial and Temporal Structure Setup ---

site_coordinates <- data.frame(
  lon = c(-90.13611111, -90.08404, -90.0274333, -89.7969444, -89.35277778),
  lat = c(29.9347222,   29.9101,   29.9644138,  29.5711111,   29.27583333)
)
# Add site IDs to coordinates 
site_coordinates$siteID_label <- paste0("Site_", seq_len(nrow(site_coordinates)))


# Calculate the distance matrix between sites using Haversine formula (great-circle distance)
# Result is in kilometers (km)
distance_matrix_km <- distm(as.matrix(site_coordinates[, c("lon", "lat")]), fun = distHaversine) / 1000

# Prepare time indexing
unique_times <- sort(unique(df$time))
df$time_idx <- match(df$time, unique_times) # Integer index for each unique time point

# Create a list where each element contains row indices for a specific site
observations_per_site_list <- split(seq_len(nrow(df)), df$siteID)

# Define key dimensions
num_sites      <- length(observations_per_site_list) # Total number of unique spatial sites (n_sites)
num_time_points<- length(unique_times)             # Total number of unique time points (T_time)
num_observations<- nrow(df)                         # Total number of observations (N)
num_predictors <- 3                                # Number of regression predictors (intercept, X1c, X2c)

# --- 3. MCMC Configuration and Initialization ---

# MCMC settings
num_iterations_total <- 10 # Total MCMC iterations
num_burnin           <- 1  # Number of initial iterations to discard (burn-in)
num_samples_to_store <- num_iterations_total - num_burnin # Number of posterior samples to keep

# Prepare predictor matrix (Design Matrix) including an intercept
X_matrix <- cbind(intercept = 1, df$X1c, df$X2c)
colnames(X_matrix) <- c("beta0", "beta1", "beta2") # Assign names for clarity

# Response vector (centered)
Y_vector <- df$Yc

# Initial values for model parameters
beta_coeffs   <- rep(0, num_predictors) # Regression coefficients
tau_Y_precision <- 1                     # Precision (1/variance) of observation error
phi_s_range   <- mean(distance_matrix_km) # Range parameter for spatial correlation (heuristic start)
sigma_s2_variance <- 1                   # Variance of spatial random effects
rho_ar1       <- 0.5                     # Autoregressive coefficient for temporal effects (AR1)
sigma_u2_variance <- 1                   # Variance of temporal random effects (innovations)
s_spatial_effects <- rep(0, num_sites)   # Spatial random effects vector
u_temporal_effects<- rep(0, num_time_points) # Temporal random effects vector

# --- 4. Storage for Posterior Samples ---

beta_posterior_samples     <- matrix(NA, nrow = num_samples_to_store, ncol = num_predictors)
colnames(beta_posterior_samples) <- colnames(X_matrix) # Store beta names
tau_Y_posterior_samples    <- numeric(num_samples_to_store)
phi_s_posterior_samples    <- numeric(num_samples_to_store)
sigma_s2_posterior_samples <- numeric(num_samples_to_store)
rho_ar1_posterior_samples  <- numeric(num_samples_to_store)
sigma_u2_posterior_samples <- numeric(num_samples_to_store)
s_effects_posterior_samples<- matrix(NA, nrow = num_samples_to_store, ncol = num_sites)
u_effects_posterior_samples<- matrix(NA, nrow = num_samples_to_store, ncol = num_time_points)

# --- 5. Proposal Standard Deviations for Metropolis-Hastings ---

prop_sd_log_sigma_s2 <- 0.1 # For log(sigma_s2)
prop_sd_log_phi_s    <- 0.1 # For log(phi_s)
prop_sd_rho_ar1      <- 0.02 # For rho (AR1 coefficient)

# --- 6. MCMC Sampling Loop ---

cat("Starting MCMC sampling...\n")
for (iter in 1:num_iterations_total) {
  cat(sprintf("\n=== Iteration %d of %d ===\n", iter, num_iterations_total))
  
  # --- a) Update regression coefficients (β) using Horseshoe prior ---
  # Residuals after removing spatial and temporal effects
  residuals_for_beta <- Y_vector - s_spatial_effects[df$siteID] - u_temporal_effects[df$time_idx]
  # Sample β using the horseshoe package (Gibbs sampler for Horseshoe)
  hs_fit      <- horseshoe(y = residuals_for_beta, X = X_matrix, nmc = 10, burn = 5, method.tau = "halfCauchy", method.sigma = "Jeffreys") # Small nmc for speed
  beta_coeffs <- as.numeric(hs_fit$BetaHat) # Point estimate (mean of samples) from horseshoe output
  
  # --- b) Update observation precision (τ_Y) via Gibbs sampling ---
  # τ_Y ~ Gamma(shape, rate)
  # Prior for τ_Y: Gamma(a_0, b_0) with a_0 = 0.01, b_0 = 0.01 (weakly informative)
  residuals_full_model <- residuals_for_beta - X_matrix %*% beta_coeffs
  shape_tau_Y <- 0.01 + length(residuals_full_model) / 2
  rate_tau_Y  <- 0.01 + sum(residuals_full_model^2) / 2
  tau_Y_precision <- rgamma(1, shape = shape_tau_Y, rate = rate_tau_Y)
  
  # --- c) Update spatial random effects (s_vec) via Gibbs sampling ---
  # s_vec ~ MVN(m_s, V_s)
  # Spatial correlation matrix R (Exponential)
  spatial_corr_matrix <- exp(-distance_matrix_km / phi_s_range)
  # Add jitter for numerical stability
  spatial_corr_matrix_jit <- spatial_corr_matrix + diag(1e-6, num_sites)
  spatial_corr_matrix_inv <- solve(spatial_corr_matrix_jit)
  
  # Sum of (Y_adj - u_t) for each site
  sum_residuals_per_site <- sapply(observations_per_site_list, function(indices_for_site) {
    sum(Y_vector[indices_for_site] - 
          (X_matrix[indices_for_site, , drop = FALSE] %*% beta_coeffs) - 
          u_temporal_effects[df$time_idx[indices_for_site]])
  })
  
  # Precision matrix for s_vec
  Q_s <- diag(tau_Y_precision * sapply(observations_per_site_list, length)) + 
    (1 / sigma_s2_variance) * spatial_corr_matrix_inv
  # Variance-covariance matrix for s_vec
  V_s <- solve(Q_s)
  # Mean vector for s_vec
  m_s <- V_s %*% (tau_Y_precision * sum_residuals_per_site)
  
  # Sample s_vec from its conditional posterior (MVN)
  # Using Cholesky decomposition for sampling from MVN: m + L * Z, where Z ~ N(0,I), L L' = V
  s_spatial_effects <- as.numeric(m_s + t(chol(V_s)) %*% rnorm(num_sites))
  s_spatial_effects <- s_spatial_effects - mean(s_spatial_effects) # Sum-to-zero constraint for identifiability
  
  # --- d) Update variance of spatial effects (σ_s²) via Metropolis-Hastings ---
  # Prior for σ_s²: Inverse-Gamma(nu/2, nu*sigma0_sq/2)
  # Proposal on log scale: log(σ_s²_prop) ~ N(log(σ_s²_current), prop_sd_log_sigma_s2^2)
  
  current_log_sigma_s2 <- log(sigma_s2_variance)
  proposed_log_sigma_s2 <- rnorm(1, current_log_sigma_s2, prop_sd_log_sigma_s2)
  proposed_sigma_s2    <- exp(proposed_log_sigma_s2)
  
  # Log-likelihood term for s_vec | σ_s², R_inv
  log_likelihood_s_current  <- -0.5 * num_sites * log(sigma_s2_variance) -
    0.5 * t(s_spatial_effects) %*% (spatial_corr_matrix_inv %*% s_spatial_effects) / sigma_s2_variance
  log_likelihood_s_proposed <- -0.5 * num_sites * log(proposed_sigma_s2) -
    0.5 * t(s_spatial_effects) %*% (spatial_corr_matrix_inv %*% s_spatial_effects) / proposed_sigma_s2
  
  # Log-prior (Inverse-Gamma for σ_s²)
  # Hyperparameters for IG prior
  prior_nu_ig        <- 0.1 
  prior_sigma0_sq_ig <- var(df$Yc, na.rm = TRUE) # Empirical estimate for scale
  
  log_prior_sigma_s2_current  <- -(prior_nu_ig / 2 + 1) * log(sigma_s2_variance) - 
    (prior_nu_ig * prior_sigma0_sq_ig) / (2 * sigma_s2_variance)
  log_prior_sigma_s2_proposed <- -(prior_nu_ig / 2 + 1) * log(proposed_sigma_s2) - 
    (prior_nu_ig * prior_sigma0_sq_ig) / (2 * proposed_sigma_s2)
  
  # Metropolis-Hastings acceptance ratio (in log scale)
  # Includes Jacobian for log transformation: log(proposed_sigma_s2) - log(current_sigma_s2)
  log_acceptance_ratio_sigma_s2 <- (log_likelihood_s_proposed + log_prior_sigma_s2_proposed) -
    (log_likelihood_s_current + log_prior_sigma_s2_current) +
    (proposed_log_sigma_s2 - current_log_sigma_s2) # Jacobian for log transform
  
  if (log(runif(1)) < log_acceptance_ratio_sigma_s2) {
    sigma_s2_variance <- proposed_sigma_s2
  }
  
  # --- e) Update spatial range parameter (φ_s) via Metropolis-Hastings ---
  # Prior for φ_s: Gamma(shape=2, rate=1/50) -> mean = 100
  # Proposal on log scale: log(φ_s_prop) ~ N(log(φ_s_current), prop_sd_log_phi_s^2)
  
  current_log_phi_s <- log(phi_s_range)
  proposed_log_phi_s <- rnorm(1, current_log_phi_s, prop_sd_log_phi_s)
  proposed_phi_s    <- exp(proposed_log_phi_s)
  
  if (proposed_phi_s > 0) { # Ensure phi is positive
    proposed_spatial_corr_matrix <- exp(-distance_matrix_km / proposed_phi_s) + diag(1e-6, num_sites)
    # Attempt to compute inverse and determinant, catch errors if matrix is singular
    proposed_spatial_corr_matrix_inv <- tryCatch(solve(proposed_spatial_corr_matrix), error = function(e) NULL)
    
    if (!is.null(proposed_spatial_corr_matrix_inv)) {
      # Log determinant of current and proposed correlation matrices
      log_det_R_current  <- determinant(spatial_corr_matrix_jit, logarithm = TRUE)$modulus
      log_det_R_proposed <- determinant(proposed_spatial_corr_matrix, logarithm = TRUE)$modulus
      
      # Quadratic form t(s) %*% R_inv %*% s
      quad_form_s_current  <- as.numeric(t(s_spatial_effects) %*% (spatial_corr_matrix_inv %*% s_spatial_effects))
      quad_form_s_proposed <- as.numeric(t(s_spatial_effects) %*% (proposed_spatial_corr_matrix_inv %*% s_spatial_effects))
      
      # Log-likelihood part for s_vec related to R(φ_s)
      log_p_conditional_s_current  <- -0.5 * (log_det_R_current + quad_form_s_current / sigma_s2_variance)
      log_p_conditional_s_proposed <- -0.5 * (log_det_R_proposed + quad_form_s_proposed / sigma_s2_variance)
      
      # Log-prior for φ_s (Gamma prior)
      log_prior_phi_s_current  <- dgamma(phi_s_range, shape = 2, rate = 1/50, log = TRUE)
      log_prior_phi_s_proposed <- dgamma(proposed_phi_s, shape = 2, rate = 1/50, log = TRUE)
      
      # Metropolis-Hastings acceptance ratio (in log scale)
      # Includes Jacobian for log transformation
      log_acceptance_ratio_phi_s <- (log_p_conditional_s_proposed + log_prior_phi_s_proposed) -
        (log_p_conditional_s_current + log_prior_phi_s_current) +
        (proposed_log_phi_s - current_log_phi_s) # Jacobian
      
      if (log(runif(1)) < log_acceptance_ratio_phi_s) {
        phi_s_range <- proposed_phi_s
        # Update these as phi_s_range has changed
        spatial_corr_matrix <- proposed_spatial_corr_matrix 
        spatial_corr_matrix_jit <- proposed_spatial_corr_matrix # Already has jitter
        spatial_corr_matrix_inv <- proposed_spatial_corr_matrix_inv
      }
    }
  }
  
  # --- f) Update temporal random effects (u_vec) via Gibbs sampling (AR(1) process) ---
  # u_t | u_{-t}, params ~ Normal
  # Prior for u_1 ~ N(0, sigma_u2 / (1-rho^2))
  # u_t | u_{t-1} ~ N(rho * u_{t-1}, sigma_u2) for t > 1
  
  for (t_idx in 1:num_time_points) {
    obs_indices_at_t <- which(df$time_idx == t_idx)
    residuals_at_t   <- Y_vector[obs_indices_at_t] - 
      (X_matrix[obs_indices_at_t, , drop = FALSE] %*% beta_coeffs) - 
      s_spatial_effects[df$siteID[obs_indices_at_t]]
    
    # Calculate precision and mean for u_t conditional on other parameters
    sum_tau_Y_residuals_t <- tau_Y_precision * sum(residuals_at_t)
    num_obs_at_t          <- length(obs_indices_at_t)
    
    if (t_idx == 1) { # First time point
      precision_u_t <- tau_Y_precision * num_obs_at_t + (1 - rho_ar1^2) / sigma_u2_variance
      if (num_time_points > 1) {
        mean_u_t <- (sum_tau_Y_residuals_t + (rho_ar1 / sigma_u2_variance) * u_temporal_effects[t_idx + 1]) / precision_u_t
      } else { # Only one time point
        mean_u_t <- sum_tau_Y_residuals_t / precision_u_t
      }
    } else if (t_idx == num_time_points) { # Last time point
      precision_u_t <- tau_Y_precision * num_obs_at_t + 1 / sigma_u2_variance
      mean_u_t      <- (sum_tau_Y_residuals_t + (rho_ar1 / sigma_u2_variance) * u_temporal_effects[t_idx - 1]) / precision_u_t
    } else { # Intermediate time points
      precision_u_t <- tau_Y_precision * num_obs_at_t + (1 + rho_ar1^2) / sigma_u2_variance
      mean_u_t      <- (sum_tau_Y_residuals_t + 
                          (rho_ar1 / sigma_u2_variance) * (u_temporal_effects[t_idx - 1] + u_temporal_effects[t_idx + 1])) / precision_u_t
    }
    
    # Sample u_t
    if(precision_u_t > 0){ # Ensure precision is positive
      u_temporal_effects[t_idx] <- rnorm(1, mean = mean_u_t, sd = sqrt(1 / precision_u_t))
    }
  }
  u_temporal_effects <- u_temporal_effects - mean(u_temporal_effects) # Sum-to-zero constraint for identifiability
  
  # --- g) Update variance of AR(1) temporal effects (σ_u²) via Gibbs sampling ---
  # σ_u² ~ Inverse-Gamma(shape, rate)
  # Prior for σ_u²: IG(a_u, b_u), e.g., a_u=2, b_u=2 (weakly informative)
  
  if (num_time_points > 1) {
    ar1_residuals <- u_temporal_effects[-1] - rho_ar1 * u_temporal_effects[-num_time_points]
    shape_sigma_u2 <- 2 + (num_time_points - 1) / 2 # IG prior shape (e.g., 2) + data term
    rate_sigma_u2  <- 2 + sum(ar1_residuals^2) / 2   # IG prior rate (e.g., 2) + data term
    sigma_u2_variance <- 1 / rgamma(1, shape = shape_sigma_u2, rate = rate_sigma_u2)
  } else { # If only one time point, update based on prior or fixed value
    # Simplified: assume u_1 ~ N(0, sigma_u2), so sigma_u2 ~ IG(a0 + 1/2, b0 + u_1^2/2)
    shape_sigma_u2 <- 2 + 1/2
    rate_sigma_u2  <- 2 + u_temporal_effects[1]^2 / 2
    sigma_u2_variance <- 1 / rgamma(1, shape = shape_sigma_u2, rate = rate_sigma_u2)
  }
  
  
  # --- h) Update AR(1) coefficient (ρ) via Metropolis-Hastings ---
  # Prior for ρ: Uniform(-1, 1) or Normal(0, small_variance) truncated to (-1,1)
  # Proposal: ρ_prop ~ N(ρ_current, prop_sd_rho_ar1^2), truncated to (-1, 1)
  
  if (num_time_points > 1) {
    proposed_rho_ar1 <- rnorm(1, rho_ar1, prop_sd_rho_ar1)
    
    if (abs(proposed_rho_ar1) < 1) { # Ensure stationarity condition |ρ| < 1
      # Log-likelihood function for AR(1) process conditional on u_vec and σ_u²
      log_posterior_ar1 <- function(rho_val, u_vals, sigma2_val) {
        if (abs(rho_val) >= 1) return(-Inf) # Penalize non-stationary proposals
        # Contribution from u_1
        log_lik <- 0.5 * log(1 - rho_val^2) - (0.5 / sigma2_val) * (1 - rho_val^2) * u_vals[1]^2
        # Contribution from u_2, ..., u_T
        ar_residuals <- u_vals[-1] - rho_val * u_vals[-length(u_vals)]
        log_lik <- log_lik - (0.5 / sigma2_val) * sum(ar_residuals^2)
        # Assuming Uniform(-1,1) prior for rho, so prior term is constant and cancels in ratio
        return(log_lik)
      }
      
      current_log_posterior_rho  <- log_posterior_ar1(rho_ar1, u_temporal_effects, sigma_u2_variance)
      proposed_log_posterior_rho <- log_posterior_ar1(proposed_rho_ar1, u_temporal_effects, sigma_u2_variance)
      
      log_acceptance_ratio_rho <- proposed_log_posterior_rho - current_log_posterior_rho
      # Proposal is symmetric, so no Hastings correction needed
      
      if (log(runif(1)) < log_acceptance_ratio_rho) {
        rho_ar1 <- proposed_rho_ar1
      }
    }
  }
  
  # --- 7. Store Post-Burn-in Samples ---
  if (iter > num_burnin) {
    sample_idx <- iter - num_burnin # Index for storing samples
    
    beta_posterior_samples[sample_idx, ] <- beta_coeffs
    tau_Y_posterior_samples[sample_idx]    <- tau_Y_precision
    phi_s_posterior_samples[sample_idx]    <- phi_s_range
    sigma_s2_posterior_samples[sample_idx] <- sigma_s2_variance
    rho_ar1_posterior_samples[sample_idx]  <- rho_ar1
    sigma_u2_posterior_samples[sample_idx] <- sigma_u2_variance
    s_effects_posterior_samples[sample_idx, ] <- s_spatial_effects
    u_effects_posterior_samples[sample_idx, ] <- u_temporal_effects
  }
} # End of MCMC loop
cat("MCMC sampling finished.\n")

# --- 8. Post-MCMC Analysis: Posterior Summaries ---

# Calculate posterior means for random effects
s_effects_posterior_mean <- colMeans(s_effects_posterior_samples)
u_effects_posterior_mean <- colMeans(u_effects_posterior_samples)

# Create a data frame of all parameter samples for easier summary
param_samples_df <- data.frame(
  beta_posterior_samples, 
  tau_Y    = tau_Y_posterior_samples,
  sigma_s2 = sigma_s2_posterior_samples,
  phi_s    = phi_s_posterior_samples,
  sigma_u2 = sigma_u2_posterior_samples,
  rho      = rho_ar1_posterior_samples
)

# --- 8.1. Parameter Scaling and Table Generation ---

# Define a scale factor for β₂ (e.g., if X2 was flow in CFS, interpret per 1,000,000 CFS).
beta2_reporting_scale_factor <- 1 
actual_beta2_scale_factor <- 1000000 

param_samples_df_scaled <- param_samples_df %>%
  mutate(
    beta2_scaled = beta2 * actual_beta2_scale_factor # Scaled effect of X2
  ) %>%
  select(-beta2) # Remove original beta2, keep scaled version for summary table

# Summarize posterior means and 95% credible intervals
posterior_summary_table <- param_samples_df_scaled %>%
  pivot_longer(everything(), names_to = "parameter_id", values_to = "value") %>%
  group_by(parameter_id) %>%
  summarise(
    posterior_mean = mean(value),
    ci_lower_2.5   = quantile(value, 0.025),
    ci_upper_97.5  = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    Parameter_Name = recode(parameter_id,
                            beta0        = "β₀ (Intercept)",
                            beta1        = "β₁ (Effect of X1c)", # X1c is centered StormSurge
                            beta2_scaled = paste0("β₂ (Effect of X2c, per ", format(actual_beta2_scale_factor, scientific = FALSE), " units)"), # X2c is centered Flow
                            tau_Y        = "τᵧ (Observation Precision)",
                            sigma_s2     = "σₛ² (Spatial Variance)",
                            phi_s        = "φₛ (Spatial Range, km)",
                            sigma_u2     = "σᵤ² (Temporal Variance)",
                            rho          = "ρ (Temporal AR1 Corr.)"
    ),
    # Format numeric values for presentation
    Mean_Formatted    = round(posterior_mean, 3),
    CI_Lower_Formatted= round(ci_lower_2.5, 3),
    CI_Upper_Formatted= round(ci_upper_97.5, 3)
  ) %>%
  # Select and order columns for the final table
  select(Parameter_Name, Mean_Formatted, CI_Lower_Formatted, CI_Upper_Formatted) %>%
  rename(Parameter = Parameter_Name, Mean = Mean_Formatted, `2.5%` = CI_Lower_Formatted, `97.5%` = CI_Upper_Formatted)


# Print the summary table using kable 
caption_text <- paste0("Posterior Means & 95% Credible Intervals (β₂ scaled for ", 
                       format(actual_beta2_scale_factor, big.mark=",", scientific=FALSE), 
                       " units of original X2)")

kable(posterior_summary_table,
      format    = "pipe", # "html", "latex", "markdown", "pipe" (good for GitHub markdown)
      escape    = FALSE, # Allows LaTeX characters if used in Parameter_Name
      caption   = caption_text,
      align     = c("l", "r", "r", "r") # l=left, r=right, c=center
) %>%
  kable_styling(full_width = FALSE, position = "center") 

################################################################################
# Section 9: MCMC Diagnostic Trace Plots

cat("\n--- Generating MCMC Diagnostic Trace Plots ---\n")

# Set up a 3x3 grid for plots
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1)) # Adjust margins: bottom, left, top, right

# Trace plot for β₀ (Intercept)
plot(beta_posterior_samples[, 1], type = 'l', main = expression("Trace: "~beta[0]), ylab = expression(beta[0]), xlab = "Iteration")

# Trace plot for β₁ (Effect of X1c)
plot(beta_posterior_samples[, 2], type = 'l', main = expression("Trace: "~beta[1]), ylab = expression(beta[1]), xlab = "Iteration")

# Trace plot for β₂ (Effect of X2c)
plot(beta_posterior_samples[, 3], type = 'l', main = expression("Trace: "~beta[2]), ylab = expression(beta[2]), xlab = "Iteration")

# Trace plot for σₛ² (Spatial Variance)
plot(sigma_s2_posterior_samples, type = 'l', main = expression("Trace: "~sigma[s]^2), ylab = expression(sigma[s]^2), xlab = "Iteration")

# Trace plot for φₛ (Spatial Range)
plot(phi_s_posterior_samples, type = 'l', main = expression("Trace: "~phi[s]), ylab = expression(phi[s]), xlab = "Iteration")

# Trace plot for σᵤ² (Temporal Variance)
plot(sigma_u2_posterior_samples, type = 'l', main = expression("Trace: "~sigma[u]^2), ylab = expression(sigma[u]^2), xlab = "Iteration")

# Trace plot for ρ (Temporal AR1 Correlation)
plot(rho_ar1_posterior_samples, type = 'l', main = expression("Trace: "~rho), ylab = expression(rho), xlab = "Iteration")

# Trace plot for Residual Variance (σ²) which is 1/τᵧ
plot(1 / tau_Y_posterior_samples, type = 'l', main = expression("Trace: Residual Variance ("~sigma^2~")"), ylab = expression(sigma^2), xlab = "Iteration")
# Reset plotting layout to default (1 plot per device)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Reset to default margins

cat("MCMC Diagnostic Trace Plots generated.\n")

# --- 10. Visualization of Model Components ---

# --- 10.1. Plot Estimated Spatial Effects (s_i) ---
spatial_effects_df <- data.frame(
  lon = site_coordinates$lon,
  lat = site_coordinates$lat,
  site_label = site_coordinates$siteID_label, # Added site labels
  s_i_posterior_mean = s_effects_posterior_mean
)

plot_spatial_effects <- ggplot(spatial_effects_df, aes(x = lon, y = lat)) +
  geom_point(aes(fill = s_i_posterior_mean), shape = 21, size = 6, stroke=1) + # Use fill for points
  geom_text(aes(label = site_label), vjust = -1.2, fontface = "bold", size = 3.5) + # Add site labels
  scale_fill_viridis_c(name = expression("Mean s"[i])) + # Continuous color scale
  coord_fixed() + # Maintain aspect ratio for geographical data
  theme_minimal(base_size = 14) +
  labs(title    = "Estimated Mean Spatial Effects (sᵢ)",
       subtitle = "Color indicates the magnitude of the site-specific effect",
       x        = "Longitude",
       y        = "Latitude",
       fill     = expression("Mean s"[i]))
print(plot_spatial_effects)

# Alternative spatial plot from user (gradient2)
plot_spatial_effects_gradient2 <- ggplot(spatial_effects_df, aes(x = lon, y = lat)) +
  geom_point(aes(color = s_i_posterior_mean), size = 5) + # Color aesthetic
  geom_text(aes(label = site_label), vjust = -1, fontface = "bold", size = 3) +
  scale_color_gradient2(low = "purple", mid = "white", high = "gold", midpoint = 0, name = expression("Mean s"[i])) +
  labs(
    title = "Estimated Mean Spatial Effects with Station IDs",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 13)
print(plot_spatial_effects_gradient2)


# --- 10.2. Plot Estimated Temporal Effects (u_t) ---
temporal_effects_df <- data.frame(
  date = unique_times,
  u_t_posterior_mean = u_effects_posterior_mean
)

plot_temporal_effects <- ggplot(temporal_effects_df, aes(x = date, y = u_t_posterior_mean)) +
  geom_line(color = "dodgerblue", linewidth = 1) +
  geom_point(color = "dodgerblue", size = 2) +
  theme_minimal(base_size = 14) +
  labs(title    = "Estimated Mean Temporal Effects (uₜ)",
       subtitle = "Shows the shared temporal trend across sites",
       x        = "Time",
       y        = expression("Mean u"[t]))
print(plot_temporal_effects)


# --- 10.3. Visualize σ_s² vs. SD(s_i) using the helper function ---
# This plot helps understand the influence of the spatial variance parameter
# on the actual variability observed in the sampled spatial effects.
if (nrow(s_effects_posterior_samples) > 0 && length(sigma_s2_posterior_samples) > 0) {
  visualize_spatial_variance_effect(sigma_s2_posterior_samples, s_effects_posterior_samples)
}


# --- 10.4. Plot Spatial Correlation Decay ---
# Compare prior median correlation decay with posterior mean correlation decay

distance_sequence <- seq(0, max(distance_matrix_km, na.rm=TRUE) * 1.1, length.out = 300) # Extend to max distance

phi_s_prior_median <- qgamma(0.5, shape = 2, rate = 1/50) 
phi_s_posterior_mean <- mean(phi_s_posterior_samples, na.rm = TRUE)

correlation_decay_df <- rbind(
  data.frame(
    Distance_km = distance_sequence,
    Correlation = exp(-distance_sequence / phi_s_posterior_mean),
    Curve_Type  = paste0("Posterior Mean (φₛ ≈ ", round(phi_s_posterior_mean, 1), " km)")
  ),
  data.frame(
    Distance_km = distance_sequence,
    Correlation = exp(-distance_sequence / phi_s_prior_median),
    Curve_Type  = paste0("Prior Median (φₛ ≈ ", round(phi_s_prior_median, 1), " km)")
  )
)

plot_correlation_decay <- ggplot(correlation_decay_df, aes(x = Distance_km, y = Correlation, color = Curve_Type, linetype = Curve_Type)) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  scale_color_manual(values = c("steelblue", "orangered")) + # Adjust colors as preferred
  scale_linetype_manual(values = c("solid", "dashed")) +
  coord_cartesian(xlim = c(0, max(distance_sequence)), ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title    = "Spatial Correlation Decay Function",
    subtitle = "Exponential Correlation: Corr(d) = exp(-d/φₛ)",
    x        = "Distance (km)",
    y        = "Spatial Correlation",
    color    = "Effective Range (φₛ)",
    linetype = "Effective Range (φₛ)"
  ) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))
print(plot_correlation_decay)

# --- 11. Quantitative Assessment of Spatial Effect Magnitude ---

range_s_effects <- range(s_effects_posterior_mean)
range_Y_original <- range(df$Y, na.rm = TRUE) # Using original Y for context before centering

spread_s_effects <- diff(range_s_effects)
spread_Y_original <- diff(range_Y_original)

cat("\n--- Quantitative Summary of Effects ---\n")
cat("Range of estimated mean spatial effects (s_i): [", round(range_s_effects[1], 2), ",", round(range_s_effects[2], 2), "]\n")
cat("Spread of estimated mean spatial effects (s_i):", round(spread_s_effects, 2), "\n")
cat("Range of original response variable Y: [", round(range_Y_original[1], 2), ",", round(range_Y_original[2], 2), "]\n")
cat("Spread of original response variable Y:", round(spread_Y_original, 2), "\n")

if (spread_Y_original > 0) {
  proportion_spatial_variation <- (spread_s_effects / spread_Y_original) * 100
  cat("Proportion of spatial effect spread relative to original Y spread:", round(proportion_spatial_variation, 1), "%\n")
}

cat("\n--- Script Finished ---\n")
