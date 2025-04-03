####################################################################################################################
####################################################################################################################
#---------------------------------------- Simulating CIR -------------------------------------------------------
####################################################################################################################
####################################################################################################################

library("GMCM")
library("ggplot2")
library(latex2exp) ### Text for plots
library("matrixStats")
set.seed(112233)

### Prediction of the short rate
n <- 10000
T <- 10
dt <- T/n
kappa <- 0.6 
theta <- 0.03 
sigma <- 0.002 
r0 <- 0.027
r <- numeric(n)
r[1] <- r0
C <- 1.5
### Pre-generate random numbers
Uv_array <- runif(n)
Zv_array <- qnorm(Uv_array)


alpha <- 0.95
z <- qnorm((1 + alpha) / 2)
tau <- numeric(n)
for (i in 1:n){
  tau[i] <- i * T/n
}
CI_upper <- numeric(n)
CI_lower <- numeric(n)
CI_upper[1] <- (r0 * exp(-kappa * tau[1]) + theta * (1 - exp(-kappa * tau[1]))) + z * sqrt(sigma^2 * (1 - exp(-2 * kappa * tau[1])) / (2 * kappa))
CI_lower[1] <- (r0 * exp(-kappa * tau[1]) + theta * (1 - exp(-kappa * tau[1]))) - z * sqrt(sigma^2 * (1 - exp(-2 * kappa * tau[1])) / (2 * kappa))



for (i in 2:n){
  # Expected value and variance for the next step
  exponent <- exp(-kappa * dt)
  m <- theta + (r[i - 1] - theta) * exponent
  s2 <- ((r[i - 1] * sigma^2 * exponent * (1 - exponent)) / kappa +
           (theta * sigma^2 * (1 - exponent)^2 / (2 * kappa)))
  
  psi <- s2 / m^2
  Uv <- Uv_array[i - 1]
  
  # Switching rule for random variable generation
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
  expected_r <- theta + (r[i - 1] - theta) * exp(-kappa * dt)
  var_r <- ((r[i - 1] * sigma^2 * exp(-kappa * dt)) / kappa) * (1 - exp(-kappa * dt)) + ((theta * sigma^2) / (2 * kappa)) * (1 - exp(-kappa * dt))^2
  CI_upper[i] <- (theta + (r0 - theta) * exp(-kappa * tau[i])) + z *
    sqrt(((r[i - 1] * sigma^2 * exp(-kappa * tau[i])) / kappa) * (1 - exp(-kappa * tau[i])) + ((theta * sigma^2) / (2 * kappa)) * (1 - exp(-kappa * tau[i]))^2)
  CI_lower[i] <- (theta + (r0 - theta) * exp(-kappa * tau[i])) - z *
    sqrt(((r[i - 1] * sigma^2 * exp(-kappa * tau[i])) / kappa) * (1 - exp(-kappa * tau[i])) + ((theta * sigma^2) / (2 * kappa)) * (1 - exp(-kappa * tau[i]))^2) 
}

CI_upper_stat <- theta + z * sqrt(sigma^2  * theta / (2 * kappa))
CI_lower_stat <- theta - z * sqrt(sigma^2  * theta / (2 * kappa))



### Plot simulated short rate


data <- data.frame(
  Time = tau,
  SR = r,
  Lower = CI_lower,
  Upper = CI_upper,
  Lower_stat = CI_lower_stat,
  Upper_stat = CI_upper_stat
)
ggplot(data, aes(x = Time)) +
  geom_line(aes(y = SR), color = "#87A96B", size = 1) +
  geom_ribbon(aes(ymin = Lower_stat, ymax = Upper_stat), alpha = 0.2, fill = "#96C8A2", col = "#96C8A2") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#801818", col = "#801818", linetype = "longdash") +
  theme_minimal(base_family = "serif") +
  labs(x = TeX("Time in Years"),
       y = TeX("Short Rate")) +
  scale_y_continuous(limits = c(0.026, 0.033)) +
  theme(
    text = element_text(family = "serif", size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0)
  ) +
  annotate(geom = "text", x = 8, y = 0.0325, label = "CI - cond. dist.",
           color = "#801818", size = 9, family = "serif") +
  annotate(geom = "text", x = 7.85, y = 0.0315, label = "CI - asymp. dist.",
           color = "#96C8A2", size = 9, family = "serif")
