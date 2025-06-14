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
n <- 10000
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
    ll <- ll + log(c) - u - v + q / 2 * (log(v) - log(u)) + log(besselI(x = z, nu = q, expon.scaled = TRUE)) + z
  }
  
  return(-ll)  # Return negative log-likelihood for minimization
}




### Function to fit the CIR model
fit_CIR <- function(r) {
  init_params <- c(log(kappa_true), log(theta_true), log(sigma_true))
  result <- optim(
    par = init_params,
    fn = CIR_loglik,
    r = r,
    dt = 1,
    method = "Nelder-Mead",
    control = list(maxit = 1000000, reltol = 1e-12)
  )
  return(c(exp(result$par[1]), exp(result$par[2]), exp(result$par[3])))
}

### Number of simulations
num_simulations <- 100

### List to store the results
simulated_series_CIR <- vector("list", num_simulations)
fitted_params_CIR <- matrix(NA, nrow = num_simulations, ncol = 3)


### Simulate and fit the CIR model
for (i in 1:num_simulations) {
  set.seed(i)  ### For reproducibility
  print(i)
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






### Plot the fitted parameters
SimCIR <- as.data.frame(fitted_params_CIR)
colnames(SimCIR) <- c("kappa", "theta", "sigma")
SimCIR$Index <- c(1:num_simulations)



#save(fitted_params_CIR, file = "fitted_params_CIR.RData")
#setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master Thesis Files")
#load("fitted_params_CIR.RData")











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
range(SimCIR$kappa)
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
  scale_y_continuous(limits = c(range(SimCIR$kappa), range(SimCIR$kappa))) + 
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
  scale_y_continuous(limits = c(range(SimCIR$theta), range(SimCIR$theta))) + 
  theme(plot.margin=unit(c(.01,.3,-0.5,-.7),"cm")) ### c(top,right,bottom,left)


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
  scale_y_continuous(limits = c(range(SimCIR$sigma), range(SimCIR$sigma))) + 
  theme(plot.margin=unit(c(.01,.3,-0.5,-.7),"cm")) ### c(top,right,bottom,left)


p3



grid.arrange(p1, p2, p3, ncol=3) 
colMeans(fitted_params_CIR, na.rm = TRUE) #  0.39943136 0.01800398 0.04012070
summary(fitted_params_CIR)
colSds(fitted_params_CIR, na.rm = TRUE) # 0.0118862935 0.0001209848 0.0004018984


#
