
#---------------------------------------- N-state Simulation -------------------------------------------------------
library(ggplot2)
library(gridExtra)
library(GMCM)
library(matrixStats)
library(scales)
### Parameters
n <- 10000
dt <- 1
N <- 5
C <- 1.5
kappa_true <- c(0.01, 0.03, 0.06, 0.08, 0.1)
theta_true <- c(0.01, 0.03, 0.06, 0.08, 0.1)
sigma_true <- c(0.01, 0.03, 0.06, 0.08, 0.1)

# Observations:
### Nelder-Mead is the only optimizer that works
### All other optimizers give errors of allocation of large vectors
### Cant initialize if i use 0.00[something]?? 
# - a standard convention is 250 working days (dt = 1/250) but i cant do that
# - Does it matter? E.g. i want to use OLS initials but they are too low.

### Results
# Pretty decent; sometimes probabilities explode even though they should be constrained to equal 1?
# Gammas are tough as transition probabilities are low and we have 1000 observations pr. sim.
Gamma_true <- diag(N)
Gamma_true[!Gamma_true] <- seq(0.01, 0.05, length = N * (N - 1))
Gamma_true <- Gamma_true / rowSums(Gamma_true)
delta_true <- solve(t(diag(N) - Gamma_true + 1), rep(1, N))
r0 <- theta_true[1]

Gamma_true
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

### Function to simulate the regime-switching CIR process
simulate_HMM <- function(n, N, kappa, theta, sigma, Gamma, delta, r0) {
  r <- numeric(n)
  S <- numeric(n)
  S[1] <- sample(x = 1:N, size = 1, prob = delta)
  r[1] <- r0
  
  for (i in 2:n) {
    S[i] <- sample(x = 1:N, size = 1, prob = Gamma[S[i - 1], ])
    r[i] <- generateCIRPathQEDisc(r[i - 1], kappa[S[i]], theta[S[i]], sigma[S[i]], 2, C = 1.5)[2]
  }
  
  return(list(r = r, S = S))
}


mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  # Gamma <- Gamma_true
  
  ### Compute the stationary distribution delta
  delta <- rep(1 / N, times = N)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[j] * theta[j] / sigma[j]^2 - 1
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j]))*sigma[j]^2)
      u <- c * x[ind[i-1]] * exp(-kappa[j])
      v <- c * x[ind[i]]
      
      # Log transition density of CIR process
      log_prob <- log(c) - u - v + q/2 * (log(v) - log(u)) + 
        log(besselI(x = 2*sqrt(u*v), nu = q, expon.scaled = TRUE)) + 
        2 * sqrt(u * v)
      
      allprobs[ind[i], j] <- exp(log_prob)
    }
  }
  
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)) {
    if (t > nrow(allprobs)) {
      stop(paste("Subscript out of bounds: t =", t, "exceeds nrow(allprobs) =", nrow(allprobs)))
    }
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  
  return(-l)
}



fit_HMM <- function(r, N, theta.star) {
  # Define the objective function for Nelder-Mead optimization
  objective_function <- function(params) {
    return(mllk(params, x = r, N = N))
  }
  
  # Use optim with Nelder-Mead
  result <- optim(
    par = theta.star, 
    fn = objective_function, 
    method = "Nelder-Mead", 
    control = list(maxit = 100000)
  )
  
  # Extract and return fitted parameters
  return(c(exp(result$par[1:N]),                  ### kappa
           exp(result$par[(N + 1):(2 * N)]),     ### theta
           exp(result$par[(2 * N + 1):(3 * N)]), ### sigma
           exp(result$par[(3 * N + 1):length(theta.star)]))) ### Gamma
}



### Number of simulations
num_simulations <- 10

### Initial parameters for fitting
theta.star <- c(log(kappa_true),
                log(theta_true),
                log(sigma_true),
                log(Gamma_true[row(Gamma_true) != col(Gamma_true)]))

### List to store the results
fitted_params_HMM <- matrix(NA, nrow = num_simulations, ncol = length(theta.star))



for (i in 1:num_simulations) {
  set.seed(i)  ### For reproducibility
  print(i)
  
  ### Simulate the CIR process
  simulation <- simulate_HMM(n, N, kappa_true, theta_true, sigma_true, Gamma_true, delta_true, r0)
  r <- simulation$r
  
  ### Fit the CIR model to the simulated series
  fitted_params_HMM[i, ] <- fit_HMM(r, N, theta.star)
}

### Print the fitted parameters
print(fitted_params_HMM)
#save(fitted_params_HMM, file = "fitted_params_HMM.RData")
load("fitted_params_HMM.RData")



SimHMM <- as.data.frame(fitted_params_HMM)
colnames(SimHMM) <- c("kappa1", "kappa2", "kappa3", "kappa4", "kappa5", 
                      "theta1", "theta2", "theta3", "theta4", "theta5", 
                      "sigma1","sigma2","sigma3","sigma4","sigma5",
                      "P(2->1)", "P(3->1)", "P(4->1)", "P(5->1)",
                      "P(1->2)", "P(3->2)", "P(4->2)", "P(5->2)",
                      "P(1->3)", "P(2->3)", "P(4->3)", "P(5->3)",
                      "P(1->4)", "P(2->4)", "P(3->4)", "P(5->4)",
                      "P(1->5)", "P(2->5)", "P(3->5)", "P(4->5)")
SimHMM$Index <- c(1:num_simulations)

#------------------------------------------ kappa -------------------------------------------------------------



param_names <- c("kappa1", "kappa2", "kappa3", "kappa4", "kappa5")
SimHMM_kappa <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(SimHMM$kappa1,
            SimHMM$kappa2,
            SimHMM$kappa3,
            SimHMM$kappa4,
            SimHMM$kappa5)
)
true_kappa <- data.frame(
  Group = c("kappa1", "kappa2", "kappa3", "kappa4", "kappa5"),
  TrueValue = kappa_true
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(kappa[1]), expression(kappa[2]), expression(kappa[3]), expression(kappa[4]), expression(kappa[5]))

ggplot(SimHMM_kappa, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = kappa_true[1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_kappa, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated speed of mean reversion, " * kappa),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=30, hjust=0)) +
  theme(axis.title = element_text(size=30)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  theme(plot.margin=unit(c(.01,.3,-0.7,-.7),"cm")) ### c(top,right,bottom,left)






#------------------------------------------ sigma -------------------------------------------------------------



param_names <- c("sigma1", "sigma2", "sigma3", "sigma4", "sigma5")
SimHMM_sigma <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(SimHMM$sigma1,
            SimHMM$sigma2,
            SimHMM$sigma3,
            SimHMM$sigma4,
            SimHMM$sigma5)
)
true_sigma <- data.frame(
  Group = c("sigma1", "sigma2", "sigma3", "sigma4", "sigma5"),
  TrueValue = sigma_true
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4]), expression(sigma[5]))

ggplot(SimHMM_sigma, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = sigma_true[1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_sigma, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated volatility, " * sigma),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=30, hjust=0)) +
  theme(axis.title = element_text(size=30)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) + 
  theme(plot.margin=unit(c(.01,.3,-0.7,-.7),"cm")) ### c(top,right,bottom,left)





#------------------------------------------ theta -------------------------------------------------------------



param_names <- c("theta1", "theta2", "theta3", "theta4", "theta5")
SimHMM_theta <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(SimHMM$theta1,
            SimHMM$theta2,
            SimHMM$theta3,
            SimHMM$theta4,
            SimHMM$theta5)
)
true_theta <- data.frame(
  Group = c("theta1", "theta2", "theta3", "theta4", "theta5"),
  TrueValue = theta_true
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]))

ggplot(SimHMM_theta, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = theta_true[1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_theta, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated long-run mean, " * theta),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=30, hjust=0)) +
  theme(axis.title = element_text(size=30)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  theme(plot.margin=unit(c(.01,.3,-0.7,-.7),"cm")) ### c(top,right,bottom,left)






#------------------------------------------ Gamma -------------------------------------------------------------

### Structure
### "P(1->2)", "P(1->3)", "P(1->4)", "P(1->5)"
### "P(2->1)", "P(2->3)", "P(2->4)", "P(2->5)"
### "P(3->1)", "P(3->2)", "P(3->4)", "P(3->5)"
### "P(4->1)", "P(4->2)", "P(4->3)", "P(4->5)"
### "P(5->1)", "P(5->2)", "P(5->3)", "P(5->4)"



### Gamma1

param_names <- c("P(1->1)", "P(1->2)", "P(1->3)", "P(1->4)", "P(1->5)")
SimHMM_Gamma1 <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(rep(NA, 10),
            as.numeric(SimHMM$"P(1->2)"),
            as.numeric(SimHMM$"P(1->3)"),
            as.numeric(SimHMM$"P(1->4)"),
            as.numeric(SimHMM$"P(1->5)"))
)
true_Gamma1 <- data.frame(
  Group = c("P(1->1)", "P(1->2)", "P(1->3)", "P(1->4)", "P(1->5)"),
  TrueValue = c(NA, Gamma_true[1,2], Gamma_true[1,3], Gamma_true[1,4], Gamma_true[1,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[11]), expression(gamma[12]), expression(gamma[13]), expression(gamma[14]), expression(gamma[15]))

HMM_Gamma1 <- ggplot(SimHMM_Gamma1, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[1,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[1,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[1,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[1,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma1, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated transition probability matrix, " * Gamma),
       x = "",
       y = expression(gamma["1j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=30, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(labels = label_number(scale = 1, accuracy = 0.0001))  + 
  theme(plot.margin=unit(c(0,.3,-1.3,.3),"cm"))


### Gamma2

param_names <- c("P(2->1)", "P(2->2)", "P(2->3)", "P(2->4)", "P(2->5)")
SimHMM_Gamma2 <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(as.numeric(SimHMM$"P(2->1)"),
            rep(NA, 10),
            as.numeric(SimHMM$"P(2->3)"),
            as.numeric(SimHMM$"P(2->4)"),
            as.numeric(SimHMM$"P(2->5)"))
)
true_Gamma2 <- data.frame(
  Group = c("P(2->1)", "P(2->2)", "P(2->3)", "P(2->4)", "P(2->5)"),
  TrueValue = c(Gamma_true[2,1], NA, Gamma_true[2,3], Gamma_true[2,4], Gamma_true[2,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[21]), expression(gamma[22]), expression(gamma[23]), expression(gamma[24]), expression(gamma[25]))

HMM_Gamma2 <- ggplot(SimHMM_Gamma2, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[2,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[2,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[2,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[2,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma2, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["2j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL)  + 
  theme(plot.margin=unit(c(-1.5,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)







### Gamma3

param_names <- c("P(3->1)", "P(3->2)", "P(3->3)", "P(3->4)", "P(3->5)")
SimHMM_Gamma3 <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(as.numeric(SimHMM$"P(3->1)"),
            as.numeric(SimHMM$"P(3->2)"),
            rep(NA, 10),
            as.numeric(SimHMM$"P(3->4)"),
            as.numeric(SimHMM$"P(3->5)"))
)
true_Gamma3 <- data.frame(
  Group = c("P(3->1)", "P(3->2)", "P(3->3)", "P(3->4)", "P(3->5)"),
  TrueValue = c(Gamma_true[3,1], Gamma_true[3,2], NA, Gamma_true[3,4], Gamma_true[3,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[31]), expression(gamma[32]), expression(gamma[33]), expression(gamma[34]), expression(gamma[35]))

HMM_Gamma3 <- ggplot(SimHMM_Gamma3, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[3,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[3,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[3,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[3,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma3, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["3j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL)  + 
  theme(plot.margin=unit(c(-1.5,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)








### Gamma4

param_names <- c("P(4->1)", "P(4->2)", "P(4->3)", "P(4->4)", "P(4->5)")
SimHMM_Gamma4 <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(as.numeric(SimHMM$"P(4->1)"),
            as.numeric(SimHMM$"P(4->2)"),
            as.numeric(SimHMM$"P(4->3)"),
            rep(NA, 10),
            as.numeric(SimHMM$"P(4->5)"))
)
true_Gamma4 <- data.frame(
  Group = c("P(4->1)", "P(4->2)", "P(4->3)", "P(4->4)", "P(4->5)"),
  TrueValue = c(Gamma_true[4,1], Gamma_true[4,2], Gamma_true[4,3], NA, Gamma_true[4,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[41]), expression(gamma[42]), expression(gamma[43]), expression(gamma[44]), expression(gamma[45]))

HMM_Gamma4 <- ggplot(SimHMM_Gamma4, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[4,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[4,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[4,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[4,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma4, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["4j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL)  + 
  theme(plot.margin=unit(c(-1.5,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)



### Gamma5

param_names <- c("P(5->1)", "P(5->2)", "P(5->3)", "P(5->4)", "P(5->5)")
SimHMM_Gamma5 <- data.frame(
  Parameter = rep(param_names, each = 10),
  Value = c(as.numeric(SimHMM$"P(5->1)"),
            as.numeric(SimHMM$"P(5->2)"),
            as.numeric(SimHMM$"P(5->3)"),
            as.numeric(SimHMM$"P(5->4)"),
            rep(NA, 10))
)
true_Gamma5 <- data.frame(
  Group = c("P(5->1)", "P(5->2)", "P(5->3)", "P(5->4)", "P(5->5)"),
  TrueValue = c(Gamma_true[5,1], Gamma_true[5,2], Gamma_true[5,3], Gamma_true[5,4], NA)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[51]), expression(gamma[52]), expression(gamma[53]), expression(gamma[54]), expression(gamma[55]))

HMM_Gamma5 <- ggplot(SimHMM_Gamma5, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[5,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[5,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[5,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[5,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma5, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["5j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = c(
    expression(j == 1), 
    expression(j == 2), 
    expression(j == 3), 
    expression(j == 4), 
    expression(j == 5)
  )) +
  coord_cartesian(ylim = NULL) + 
  theme(plot.margin=unit(c(-1.5,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right) 


HMM_Gamma1
HMM_Gamma2
HMM_Gamma3
HMM_Gamma4
HMM_Gamma5
colMeans(fitted_params_HMM)
colSds(fitted_params_HMM)

