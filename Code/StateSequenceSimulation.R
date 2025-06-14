set.seed(123) # For reproducibility

library(ggplot2)
library(gridExtra)
library(GMCM)
library(matrixStats)
library(latex2exp)
library(scales)
### Parameters
n <- 10000 #27.4 years
dt <- 1
N <- 5
C <- 1.5
# Define the state vectors for the CIR model
kappa_true <- c(0.04, 0.020, 0.025, 0.0010, 0.005)  # Fast mean reversion
theta_true <- c(0.001, 0.01, 0.04, 0.05, 0.09)  # Short-term yields
sigma_true <- c(0.004, 0.003, 0.001, 0.003, 0.004)  # Moderate volatility
Gamma_true <- diag(N)
Gamma_true[!Gamma_true] <- exp(-5)
Gamma_true <- Gamma_true / rowSums(Gamma_true)
delta_true <- solve(t(diag(N) - Gamma_true + 1), rep(1, N))
r0 <- 0.04


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
simulate_HMM <- function(n, N, kappa, theta, sigma, Gamma, delta, r0, C) {
  r <- numeric(n)
  S <- numeric(n)
  S[1] <- sample(x = 1:N, size = 1, prob = delta)
  r[1] <- r0
  
  for (i in 2:n) {
    S[i] <- sample(x = 1:N, size = 1, prob = Gamma[S[i - 1], ])
    r[i] <- generateCIRPathQEDisc(r[i - 1], kappa[S[i]], theta[S[i]], sigma[S[i]], 2, C)[2]
  }
  
  return(list(r = r, S = S))
}

set.seed(123) # For reproducibility

# Simulate 1 path of the regime-switching CIR process
simulation <- simulate_HMM(
  n = n, 
  N = N, 
  kappa = kappa_true, 
  theta = theta_true, 
  sigma = sigma_true, 
  Gamma = Gamma_true, 
  delta = delta_true, 
  r0 = r0, 
  C = C
)

# Extract short rates
r <- simulation$r

# Create a dataframe for plotting
data <- data.frame(
  Time = (1:n)/365.25,
  ShortRate = r
)

# Plot the short rate path with centered title and LaTeX-style text
p <- ggplot(data, aes(x = Time, y = ShortRate)) +
  geom_line(color = "#901a1E") +
  labs(
    y = TeX("Short Rate"),
    x = TeX("Time (years)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "serif"),
    axis.title = element_text(family = "serif")
  )


p_density <- ggplot(data, aes(x = ShortRate)) +
  geom_density(aes(x=r),colour="#901a1E") +
  labs(
    x = TeX("Short Rate"),
    y = "Density (%)"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +  # Normalize to percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "serif"),
    axis.title = element_text(family = "serif")
  )

# Arrange both plots vertically
grid.arrange(p_density, p, ncol = 1, heights = c(1, 2))
