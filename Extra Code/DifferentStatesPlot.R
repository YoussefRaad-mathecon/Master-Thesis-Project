library(ggplot2)
library(latex2exp)

# Set parameters
total_steps <- 10000  # Number of time steps
dt <- 1  # Time step
r0 <- 0.05  # Initial interest rate
years <- total_steps / 252  # Convert steps to years (assuming 252 trading days per year)

# Define economic regimes
regimes <- list(
  "Bad Recession" = list(kappa = 0.01, theta = 0.02, sigma = 0.001, C = 1.5),
  "Recession" = list(kappa = 0.015, theta = 0.03, sigma = 0.002, C = 1.5),
  "Recovery" = list(kappa = 0.02, theta = 0.05, sigma = 0.003, C = 1.5),
  "Expansion" = list(kappa = 0.03, theta = 0.07, sigma = 0.004, C = 1.5),
  "Good Expansion" = list(kappa = 0.04, theta = 0.1, sigma = 0.005, C = 1.5)
)

# Function to simulate CIR process
generateCIRPath <- function(r0, kappa, theta, sigma, n, C) {
  r <- numeric(n)
  r[1] <- r0
  exponent <- exp(-kappa * dt)
  
  Uv_array <- runif(n)
  Zv_array <- rnorm(n)
  
  for (i in 2:n) {
    m <- theta + (r[i - 1] - theta) * exponent
    s2 <- ((r[i - 1] * sigma^2 * exponent * (1 - exponent)) / kappa +
             (theta * sigma^2 * (1 - exponent)^2 / (2 * kappa)))
    psi <- s2 / m^2
    Uv <- Uv_array[i - 1]
    
    if (psi <= C) {
      Zv <- Zv_array[i - 1]
      b2 <- (2 / psi) - 1 + sqrt(2 / psi) * sqrt((2 / psi) - 1)
      a <- m / (1 + b2)
      r_next <- a * (sqrt(b2) + Zv)^2
    } else {
      p <- (psi - 1) / (psi + 1)
      beta <- (1 - p) / m
      if (Uv <= p) {
        r_next <- 0
      } else {
        r_next <- (1 / beta) * log((1 - p) / (1 - Uv))
      }
    }
    
    r[i] <- r_next
  }
  return(r)
}

# Simulate CIR paths for each regime
set.seed(42)
simulated_paths <- lapply(regimes, function(params) {
  generateCIRPath(r0, params$kappa, params$theta, params$sigma, total_steps, params$C)
})

# Convert to data frame for plotting
df <- data.frame(Year = rep(seq(0, years, length.out = total_steps), times = length(simulated_paths)),
                 ShortRate = unlist(simulated_paths),
                 Regime = rep(names(simulated_paths), each = total_steps))

# Plot results
ggplot(df, aes(x = Year, y = ShortRate, color = Regime)) +
  geom_line() +
  labs(title = "Simulated Interest Rate Paths for Different Economic Regimes",
       x = TeX("Year"), y = TeX("Short\\ Rate")) +
  theme_minimal(base_family = "serif") +
  theme(
    text = element_text(family = "serif", size = 14),
    plot.title = element_text(hjust = 0.5)
  )
