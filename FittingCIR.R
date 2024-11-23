####################################################################################################################
####################################################################################################################
#---------------------------------------- Fitting CIR ----------------------------------------------------------
####################################################################################################################
####################################################################################################################

library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots
library("fHMM") ### Model fitting
library(Quandl) ### Data
library(dplyr) ### Data manipulations
library(lubridate) ### Time-variables
library(stats) ### ACF plots
library(matrixcalc) ### Matrix calculations
library("RColorBrewer") ### Colors
library(latex2exp) ### Text for plots
library(matrixStats) ### ColSds
set.seed(123)
#----------------------------------------------data check-----------------------------------------------------#
sum(yields_df$"3M" <= 0)  # Check for non-positive values
# Remove rows using logical indexing
yields_df <- yields_df[yields_df$"3M" > 0, ]

# Check the result
sum(yields_df$"3M" <= 0)  # Should return 0

#----------------------------------------------Log-likelihood--------------------------------------------------#
CIR_loglik <- function(params, r, dt = 1) {
  # Extract parameters
  kappa <- exp(params[1])  # Mean reversion speed
  theta <- exp(params[2])  # Long-term mean level
  sigma <- exp(params[3])  # Volatility coefficient
  
  # Check for valid parameters
  if (kappa <= 0 || theta <= 0 || sigma <= 0) {
    return(1e10)  # Large penalty
  }
  
  n <- length(r)
  ll <- 0  # Initialize log-likelihood
  
  # Loop through observations and compute log-likelihood
  for (i in 2:n) {
    # Define constants for the chi-square distribution
    n_t_T <- (4 * kappa * exp(-kappa * dt)) / (sigma^2 * (1 - exp(-kappa * dt)))
    scaling_factor <- exp(-kappa * dt) / n_t_T
    df <- (4 * kappa * theta) / sigma^2
    ncp <- n_t_T * r[i - 1]
    
    # Scaled observation
    scaled_x <- (1/scaling_factor) * r[i]
    
    # Compute log-density
    ll <- ll + dchisq(scaled_x, df = df, ncp = ncp, log = TRUE) - log(scaling_factor)
  }
  
  return(-ll)  # Return negative log-likelihood for minimization
}


# CIR initial parameters estimation using OLS
#----------------------------------------------Initial parameters-------------------------------------------------#
# Define function for OLS estimation of CIR parameters
estimate_initial_params <- function(data, dt = 1) {
  # Remove first and last observations for differencing
  x <- data[-length(data)]
  dx <- diff(data)
  
  # Transform the series
  dx_transformed <- dx / sqrt(x)
  regressors <- cbind(dt / sqrt(x), dt * sqrt(x))
  
  # Perform OLS regression
  drift <- solve(t(regressors) %*% regressors) %*% t(regressors) %*% dx_transformed
  
  # Extract parameters
  kappa <- -drift[2]
  theta <- -drift[1] / drift[2]
  residuals <- dx_transformed - regressors %*% drift
  sigma <- sqrt(var(residuals) / dt)
  
  # Return the estimated parameters as a named list
  list(kappa = kappa, theta = theta, sigma = sigma)
}


# Example usage with yields_df$"3M"
initial_params <- estimate_initial_params(as.numeric(yields_df$"3M"), dt = 1)

initial_params
# Print the estimated parameters
print(initial_params)

# Use these initial parameters for your optimization
init_params <- c(log(initial_params$kappa), log(initial_params$theta), log(initial_params$sigma))


#--------------------------------------------ACF for Kappa:-----------------------------------------------------#
# Function to estimate kappa using ACF
estimate_kappa_acf <- function(data, dt = 1) {
  # Compute ACF for lag 1
  acf_value <- acf(data, lag.max = 1, plot = FALSE)$acf[2]
  
  # Check if ACF is positive to avoid issues with log
  if (acf_value <= 0) {
    stop("ACF is non-positive, unable to estimate kappa.")
  }
  
  # Estimate kappa
  kappa <- -log(acf_value) / dt
  return(kappa)
}

# Example usage with yields_df$"3M"
kappa_acf <- estimate_kappa_acf(as.numeric(yields_df$"3M"), dt = 1)

# Print the estimated kappa
cat("Estimated kappa using ACF method:", kappa_acf, "\n")

# Insert it in the original OLS vector if wanted
init_params[1] <- log(kappa_acf)


#----------------------------------------------Random Guesses------------------------------------------------------------#
#init_params <- c(1, 1, 1)
#init_params <- c(0.005, 0.005, 0.005)
#--------------------------------------------Fitting-------------------------------------------------------------------#

# Continue with the CIR log-likelihood optimization as before
mod <- nlm(f = CIR_loglik, p = init_params, r = as.numeric(yields_df$"3M"), print.level = 2, gradtol = 1e-12, iterlim = 10000)
print(mod)

# Extract and print fitted parameters
kappa_hat <- exp(mod$estimate[1])
theta_hat <- exp(mod$estimate[2])
sigma_hat <- exp(mod$estimate[3])
loglikelihood <- mod$minimum
cat("Fitted parameters:\n")
cat("kappa:", kappa_hat, "\n") # 0.0005323919 
cat("theta:", theta_hat, "\n") # 0.020025 
cat("sigma:", sigma_hat, "\n") # 0.004474768 
cat("loglikelihood:", loglikelihood, "\n") #-59118.33 with OLS
#-59146.24 with ACF
#-59133.05 with 1. guess

