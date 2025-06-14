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
  x <- data[-length(data)]
  dx <- diff(data)
  
  dx_transformed <- dx / sqrt(x)
  regressors <- cbind(dt / sqrt(x), dt * sqrt(x))
  
  XtX_inv <- solve(t(regressors) %*% regressors)
  drift <- XtX_inv %*% t(regressors) %*% dx_transformed
  residuals <- dx_transformed - regressors %*% drift
  sigma2 <- var(residuals)
  
  cov_matrix <- XtX_inv * sigma2  # element-wise is okay here
  
  std_errors_drift <- sqrt(diag(cov_matrix))
  
  # Extract drift coefficients
  beta0 <- drift[1]
  beta1 <- drift[2]
  
  # Delta method for standard errors of transformed parameters
  # kappa = -beta1
  se_kappa <- sqrt(cov_matrix[2, 2])
  
  # theta = -beta0 / beta1
  dtheta_dbeta0 <- -1 / beta1
  dtheta_dbeta1 <- beta0 / (beta1^2)
  grad_theta <- c(dtheta_dbeta0, dtheta_dbeta1)
  var_theta <- t(grad_theta) %*% cov_matrix %*% grad_theta
  se_theta <- sqrt(var_theta)
  
  # CIR volatility
  se_sigma <- (1 / (2 * sqrt(dt * sigma2))) * sqrt(2 * sigma2^2 / (length(residuals) - 2))
  
  list(
    kappa = -beta1,
    theta = -beta0 / beta1,
    sigma = sqrt(sigma2 / dt),
    std_errors = c(se_kappa = se_kappa, se_theta = se_theta, se_sigma = se_sigma),
    drift_coeffs = drift,
    drift_std_errors = std_errors_drift,
    cov_matrix = cov_matrix
  )
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
estimate_kappa_acf <- function(data, dt = 1 / 252) {
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
kappa_acf <- estimate_kappa_acf(as.numeric(yields_df$"3M"), dt = 1/252)

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
  dt = 1 / 252, 
  method = "Nelder-Mead", 
  control = list(maxit = 1000000, reltol = 1e-16),
  hessian = TRUE
)
print(mod)

# Extract and print fitted parameters
kappa_hat <- exp(mod$par[1])
theta_hat <- exp(mod$par[2])
sigma_hat <- exp(mod$par[3])
loglikelihood <- mod$value
cat("Fitted parameters:\n")                          #dt=1 but scaled correctly
cat("kappa:", format(kappa_hat, digits = 15), "\n")  #w/0 acf 0.13451935626599 
cat("theta:", format(theta_hat, digits = 15), "\n")  #w/0 acf  0.0199792346688647
cat("sigma:", format(sigma_hat, digits = 15), "\n")  #w/0 acf  0.071028365376124
cat("loglikelihood:", format(loglikelihood, digits = 15), "\n") # w/0 acf -59157.9352174842 





cov_log <- solve(mod$hessian)  # Inverse of Hessian
se_log <- sqrt(diag(cov_log))  # Standard errors on log-scale

# Extract log estimates
log_estimates <- mod$par
# Compute SEs on original scale
se_original <- exp(log_estimates) * se_log

# Parameter estimates
kappa_hat <- exp(mod$par[1])
theta_hat <- exp(mod$par[2])
sigma_hat <- exp(mod$par[3])

# Standard errors
kappa_se <- se_original[1]
theta_se <- se_original[2]
sigma_se <- se_original[3]

# Print results
cat("kappa_hat:", kappa_hat, "SE:", kappa_se, "\n")
cat("theta_hat:", theta_hat, "SE:", theta_se, "\n")
cat("sigma_hat:", sigma_hat, "SE:", sigma_se, "\n")
















compute_CIR_residuals <- function(data, kappa, theta, sigma, dt = 1/252) {
  # Remove first and last observations for differencing
  x <- data[-length(data)]
  dx <- diff(data)
  
  # Transform the series using the CIR discretization
  dx_transformed <- dx / sqrt(x)
  regressors <- cbind(dt / sqrt(x), dt * sqrt(x))
  
  # Fitted drift component based on estimated parameters
  fitted_drift <- -(kappa * theta * dt / sqrt(x)) + (kappa * dt * sqrt(x))
  
  # Residuals
  residuals <- dx_transformed - fitted_drift
  residuals <- residuals / (sigma * sqrt(dt))
  
  return(residuals)
}

r_data <- as.numeric(yields_df$"3M")
residuals <- compute_CIR_residuals(r_data, kappa_hat, theta_hat, sigma_hat)


# Define theme color
theme_col <- "#901a1E"

# Remove first date to align with residuals
residuals_df <- data.frame(
  DateCont = yields_df$DateCont[-1],
  Residuals = residuals
)

ggplot(residuals_df, aes(x = DateCont, y = Residuals)) +
  geom_line(color = "#901a1E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  xlab(TeX("Year")) +
  ylab(TeX("Residual")) +
  scale_x_continuous(
    breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7)
  ) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
qqnorm(residuals, col = "#901a1E", main = "")
qqline(residuals, col = "black", lwd = 2)

# Compute 1st and 99th percentiles
q_low <- quantile(residuals, 0.01)
q_high <- quantile(residuals, 0.99)

# Filter residuals to keep only those within the 1%–99% range
residuals_trimmed <- residuals[residuals >= q_low & residuals <= q_high]

# Plot histogram 1 to 99 quantile
ggplot(data.frame(residuals_trimmed), aes(x = residuals_trimmed)) +
  geom_histogram(
    bins = 100,
    fill = "#901a1E",
    color = "white"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  xlab("Residuals") +
  ylab("Frequency") +
  ggtitle("Histogram of CIR Residuals (1st–99th Percentile)") +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13)
  )





library(ggplot2)
library(latex2exp)
library(patchwork)

# Define theme color
theme_col <- "#901a1E"

# Ensure residuals is a numeric vector
residuals <- as.numeric(residuals)

# Create residuals_df for time series plot
residuals_df <- data.frame(
  DateCont = yields_df$DateCont[-1],
  Residuals = residuals
)

# ------------------ 1. Time Series Plot ------------------
p1 <- ggplot(residuals_df, aes(x = DateCont, y = Residuals)) +
  geom_line(color = theme_col) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  xlab(TeX("Year")) +
  ylab(TeX("Residual")) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

# ------------------ 2. Q-Q Plot ------------------
p2 <- ggplot(data.frame(residuals = residuals), aes(sample = residuals)) +
  stat_qq(color = theme_col) +
  stat_qq_line(color = "black", linetype = "dashed") +
  theme_bw() +
  xlab("Theoretical Quantiles") +
  ylab("Sample Quantiles") +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13)
  )

# ------------------ 3. Histogram (1%–99% Quantiles) ------------------
q_low <- quantile(residuals, 0.01)
q_high <- quantile(residuals, 0.99)
#residuals_trimmed <- residuals[residuals >= q_low & residuals <= q_high]

p3 <- ggplot(data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(bins = 100, fill = theme_col, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  xlab("Residuals") +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13)
  )
p3
# ------------------ Combine all 3 plots ------------------
(p1 / (p2 | p3)) + plot_layout(heights = c(1, 1))



