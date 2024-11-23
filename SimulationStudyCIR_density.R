####################################################################################################################
####################################################################################################################
#---------------------------------------- CIR Simulation -------------------------------------------------------
####################################################################################################################
####################################################################################################################

library("GMCM")
library("ggplot2")
library("matrixStats")
library("gridExtra")

### Set parameters for the CIR process
n <- 1000
dt <- 1
kappa_true <- 0.02
theta_true <- 0.06
sigma_true <- 0.003
C <- 1.5
r0 <- theta_true

### Function to simulate CIR process
generateCIRPathQEDisc <- function(r0, kappa, theta, sigma, n, C) {
  r <- numeric(n)
  r[1] <- r0
  
  exponent <- exp(-kappa * dt)
  
  # Pre-generate random numbers
  Uv_array <- runif(n)
  Zv_array <- qnorm(Uv_array)
  
  for (i in 2:(n)) {
    m <- theta + (r[i - 1] - theta) * exponent
    s2 <- ((r[i - 1] * sigma^2 * exponent * (1 - exponent)) / kappa +
             (theta * sigma^2 * (1 - exponent)^2 / (2 * kappa)))
    
    psi <- s2 / m^2
    Uv <- Uv_array[i - 1]
    
    # Switching rule
    if (psi <= C) {
      Zv <- Zv_array[i - 1]
      b2 <- (2 / psi) - 1 + sqrt(2 / psi) * sqrt((2 / psi) - 1)
      a <- m / (1 + b2)
      r_next <- a * (sqrt(b2) + Zv)^2
    } else {  # psi > C
      p <- (psi - 1) / (psi + 1)
      beta <- (1 - p) / m
      if (0 <= Uv && Uv <= p) {
        r_next <- 0
      } else if (p < Uv && Uv <= 1) {
        r_next <- (1 / beta) * log((1 - p) / (1 - Uv))
      }
    }
    
    r[i] <- r_next
  }
  
  return(r)
}

CIR_loglik <- function(params, r, dt = 1) {
  # Extract parameters
  kappa <- exp(params[1])  # Mean reversion speed
  theta <- exp(params[2])  # Long-term mean level
  sigma <- exp(params[3])  # Volatility coefficient
  
  # Check for valid parameters
  if (kappa <= 0 || theta <= 0 || sigma <= 0) {
    return(1e6)  # Large penalty
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




### Function to fit the CIR model
fit_CIR <- function(r) {
  init_params <- c(log(kappa_true), log(theta_true), log(sigma_true))
  result <- nlm(f = CIR_loglik, p = init_params, r = r, print.level = 0, gradtol = 1e-12, iterlim = 10000)
  return(c(exp(result$estimate[1]), exp((result$estimate[2])), exp(result$estimate[3])))
}

### Number of simulations
num_simulations <- 100

### List to store the results
simulated_series_CIR <- vector("list", num_simulations)
fitted_params_CIR <- matrix(NA, nrow = num_simulations, ncol = 3)


### Simulate and fit the CIR model
for (i in 1:num_simulations) {
  set.seed(i)  ### For reproducibility
  
  ### Simulate the CIR process
  simulated_series_CIR[[i]] <- generateCIRPathQEDisc(r0, kappa_true, theta_true, sigma_true, n, C) 
  
  ### Fit the CIR model to the simulated series
  fitted_params_CIR[i, ] <- fit_CIR(simulated_series_CIR[[i]])
}


CIRFit <- numeric(3)
for (i in 1:3){
  CIRFit[i] <- mean(fitted_params_CIR[,i])
}
colMeans(fitted_params_CIR)






####################################################################################################################
#------------------------------------------ Load RData -------------------------------------------------------------
####################################################################################################################


#save(fitted_params_CIR, file = "fitted_params_CIR.RData")
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master Thesis Files")
#load("fitted_params_CIR.RData")








### Plot the fitted parameters
SimCIR <- as.data.frame(fitted_params_CIR)
colnames(SimCIR) <- c("kappa", "theta", "sigma")
SimCIR$Index <- c(1:num_simulations)











#------------------------------------------ kappa -------------------------------------------------------------


param_names <- c("kappa")
SimCIR_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimCIR$kappa)
)
true_vals <- data.frame(
  Group = c("kappa"),
  TrueValue = c(kappa_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(kappa))


p1 <- ggplot(SimCIR_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, show.legend = FALSE, size = 7) +
  scale_color_manual(values=c("#404080", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = kappa_true, col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression(""),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=10)) +
  scale_x_discrete(labels = latex_labels) + 
  scale_y_continuous(limits = c()) + 
  theme(plot.margin=unit(c(0,.3,-0.5,-.7),"cm")) ### c(top,right,bottom,left)


p1



#------------------------------------------ theta -------------------------------------------------------------


param_names <- c("theta")
SimCIR_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimCIR$theta)
)
true_vals <- data.frame(
  Group = c("theta"),
  TrueValue = c(theta_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(theta))


p2 <- ggplot(SimCIR_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, show.legend = FALSE, size = 7) +
  scale_color_manual(values=c("#FFBCD9", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = theta_true, col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression(""),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=10)) +
  scale_x_discrete(labels = latex_labels) +
  scale_y_continuous(limits = c()) + 
  theme(plot.margin=unit(c(.01,.3,-0.5,-.7),"cm")) ### c(top,right,bottom,left)




p2
#------------------------------------------ sigma -------------------------------------------------------------


param_names <- c("sigma")
SimCIR_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimCIR$sigma)
)
true_vals <- data.frame(
  Group = c("sigma"),
  TrueValue = c(sigma_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(sigma))


p3 <- ggplot(SimCIR_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 7, show.legend = FALSE) +
  scale_color_manual(values=c("#87A96B", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = sigma_true, col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression(""),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=10)) +
  scale_x_discrete(labels = latex_labels) +
  scale_y_continuous(limits = c()) + 
  theme(plot.margin=unit(c(.01,.3,-0.5,-.7),"cm")) ### c(top,right,bottom,left)


p3



grid.arrange(p1, p2, p3, ncol=3)
colMeans(fitted_params_CIR)
colSds(fitted_params_CIR)
