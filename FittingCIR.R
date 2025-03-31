####################################################################################################################
####################################################################################################################
#---------------------------------------- Fitting CIR ----------------------------------------------------------
####################################################################################################################
######################  ##############################################################################################

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

set.seed(112233)
#----------------------------------------------data check-----------------------------------------------------#
sum(yields_df$"3M" <= 0)  # Check for non-positive values
#yields_df <- yields_df[yields_df$"3M" > 0, ]

# Check the result
sum(yields_df$"3M" <= 0)  # Should return 0

#----------------------------------------------Log-likelihood--------------------------------------------------#
CIR_loglik <- function(params, r, dt = 1/252) {
  # Extract parameters
  kappa <- exp(params[1])  # Mean reversion speed
  theta <- exp(params[2])  # Long-term mean level
  sigma <- exp(params[3])  # Volatility coefficient
  
  n <- length(r)
  
  # Initialize log-likelihood
  ll <- 0
  # Calculate constants
  # Loop through observations and compute log-likelihood
  for (i in 2:n) {
    c <- (2 * kappa) / (sigma^2 * (1 - exp(-kappa * dt)))
    q <- (2 * kappa * theta) / sigma^2 - 1
    u <- c * r[i-1] * exp(-kappa * dt)
    v <- c * r[i]
    z <- 2 * sqrt(u * v)
    # Calculate log-likelihood using the provided formula 
    ll <- ll + log(c) - u - v + q / 2 * (log(v) - log(u)) + log(besselI(x = z, nu = q, expon.scaled = TRUE)) + abs(z)
  }
  
  return(-ll)  # Return negative log-likelihood for minimization
}


# CIR initial parameters estimation using OLS
#----------------------------------------------Initial parameters-------------------------------------------------#
# Define function for OLS estimation of CIR parameters
estimate_initial_params <- function(data, dt = 1/252) {
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
initial_params <- estimate_initial_params(as.numeric(yields_df$"3M"), dt = 1/252)

initial_params
# Print the estimated parameters
print(initial_params) 
### with dt = 1
#kappa 0.000820703
#theta 0.0247371
#sigma 0.004693142

### with dt = 1/252
#kappa 0.2068172
#theta 0.0247371
#sigma 0.07450132


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

# fast and close to same
#----------------------------------------------Random Guesses------------------------------------------------------------#
#init_params <- c(1, 1, 1) long and not correct
#init_params <- c(0.005, 0.005, 0.005) long and not correct
#--------------------------------------------Fitting---------------------------------------------------------------------#
# Continue with the CIR log-likelihood optimization as before
mod <- optim(
  par = init_params, 
  fn = CIR_loglik, 
  r = as.numeric(yields_df$"3M"), 
  dt = 1, 
  method = "Nelder-Mead", 
  control = list(maxit = 1000000, reltol = 1e-16)
)
print(mod)

# Extract and print fitted parameters
kappa_hat <- exp(mod$par[1])*252
theta_hat <- exp(mod$par[2])
sigma_hat <- exp(mod$par[3])*sqrt(252)
loglikelihood <- mod$value
cat("Fitted parameters:\n")                          #dt=1 but scaled correctly
cat("kappa:", format(kappa_hat, digits = 15), "\n")  #w/0 acf 0.13451935626599 
cat("theta:", format(theta_hat, digits = 15), "\n")  #w/0 acf  0.0199792346688647
cat("sigma:", format(sigma_hat, digits = 15), "\n")  #w/0 acf  0.071028365376124
cat("loglikelihood:", format(loglikelihood, digits = 15), "\n") # w/0 acf -59157.9352174842 


 