# CIR Log-likelihood function (from your post)
CIR_loglik <- function(params, r, dt = 1/252) {
  kappa <- exp(params[1])  # Mean reversion speed
  theta <- exp(params[2])  # Long-term mean level
  sigma <- exp(params[3])  # Volatility coefficient
  
  n <- length(r)
  ll <- 0
  for (i in 2:n) {
    c <- (2 * kappa) / (sigma^2 * (1 - exp(-kappa * dt)))
    q <- (2 * kappa * theta) / sigma^2 - 1
    u <- c * r[i-1] * exp(-kappa * dt)
    v <- c * r[i]
    z <- 2 * sqrt(u * v)
    ll <- ll + log(c) - u - v + q / 2 * (log(v) - log(u)) +
      log(besselI(x = z, nu = q, expon.scaled = TRUE)) + abs(z)
  }
  return(-ll)  # Negative log-likelihood for minimization
}

# Wrapper for fixed-parameter likelihood surface
CIR_loglik_profile <- function(kappa, theta, sigma, r, dt = 1/252) {
  params <- log(c(kappa, theta, sigma))
  -CIR_loglik(params, r, dt) # Return log-likelihood (not negative)
}

# MLEs from your question
kappa_hat <- 0.206817
theta_hat <- 0.024737
sigma_hat <- 0.0745013

# Grids for each parameter
kappa_seq <- seq(0.2, 0.7, length.out = 50)
theta_seq <- seq(0.02, 0.1, length.out = 50)
sigma_seq <- seq(0.05, 0.25, length.out = 50)

# Provide your observed rates, e.g.
r <- as.numeric(yields_df$"3M")

# 1. Log-likelihood surface: Kappa vs Theta (fix Sigma)
loglik_mat1 <- matrix(NA, length(kappa_seq), length(theta_seq))
for (i in seq_along(kappa_seq)) {
  for (j in seq_along(theta_seq)) {
    loglik_mat1[i, j] <- CIR_loglik_profile(
      kappa = kappa_seq[i],
      theta = theta_seq[j],
      sigma = sigma_hat,
      r = r
    )
  }
}

# 2. Log-likelihood surface: Kappa vs Sigma (fix Theta)
loglik_mat2 <- matrix(NA, length(kappa_seq), length(sigma_seq))
for (i in seq_along(kappa_seq)) {
  for (j in seq_along(sigma_seq)) {
    loglik_mat2[i, j] <- CIR_loglik_profile(
      kappa = kappa_seq[i],
      theta = theta_hat,
      sigma = sigma_seq[j],
      r = r
    )
  }
}

# 3. Log-likelihood surface: Theta vs Sigma (fix Kappa)
loglik_mat3 <- matrix(NA, length(theta_seq), length(sigma_seq))
for (i in seq_along(theta_seq)) {
  for (j in seq_along(sigma_seq)) {
    loglik_mat3[i, j] <- CIR_loglik_profile(
      kappa = kappa_hat,
      theta = theta_seq[i],
      sigma = sigma_seq[j],
      r = r
    )
  }
}



par(mfrow = c(1, 3))
persp(
  x = kappa_seq, y = theta_seq, z = loglik_mat1,
  xlab = expression(kappa), ylab = expression(theta), zlab = expression(ell),
  theta = 40, phi = 30, col = "lightblue", ticktype = "detailed",
  main = expression("Log-Likelihood: " * kappa * " vs " * theta * "  (" * sigma * " fixed)")
)
persp(
  x = kappa_seq, y = sigma_seq, z = loglik_mat2,
  xlab = expression(kappa), ylab = expression(sigma), zlab = expression(ell),
  theta = 40, phi = 30, col = "lightgreen", ticktype = "detailed",
  main = expression("Log-Likelihood: " * kappa * " vs " * sigma * "  (" * theta * " fixed)")
)
persp(
  x = theta_seq, y = sigma_seq, z = loglik_mat3,
  xlab = expression(theta), ylab = expression(sigma), zlab = expression(ell),
  theta = 40, phi = 30, col = "salmon", ticktype = "detailed",
  main = expression("Log-Likelihood: " * theta * " vs " * sigma * "  (" * kappa * " fixed)")
)
par(mfrow = c(1,1))


if (requireNamespace("plotly", quietly = TRUE)) {
  library(plotly)
  
  bigfont <- list(size = 22)  # You can increase/decrease as needed
  
  # 1. Kappa vs Theta
  fig1 <- plot_ly(x = kappa_seq, y = theta_seq, z = ~loglik_mat1) %>%
    add_surface(colorbar = list(title = list(text = "\u2113", font = bigfont))) %>%
    layout(
      title = "Log-Likelihood Surface: \u03BA vs. \u03B8 (\u03C3 fixed)",
      scene = list(
        xaxis = list(title = "\u03BA", titlefont = bigfont),
        yaxis = list(title = "\u03B8", titlefont = bigfont),
        zaxis = list(title = "\u2113", titlefont = bigfont)
      )
    )
  
  # 2. Kappa vs Sigma
  fig2 <- plot_ly(x = kappa_seq, y = sigma_seq, z = ~loglik_mat2) %>%
    add_surface(colorbar = list(title = list(text = "\u2113", font = bigfont))) %>%
    layout(
      title = "Log-Likelihood Surface: \u03BA vs. \u03C3 (\u03B8 fixed)",
      scene = list(
        xaxis = list(title = "\u03BA", titlefont = bigfont),
        yaxis = list(title = "\u03C3", titlefont = bigfont),
        zaxis = list(title = "\u2113", titlefont = bigfont)
      )
    )
  
  # 3. Theta vs Sigma
  fig3 <- plot_ly(x = theta_seq, y = sigma_seq, z = ~loglik_mat3) %>%
    add_surface(colorbar = list(title = list(text = "\u2113", font = bigfont))) %>%
    layout(
      title = "Log-Likelihood Surface: \u03B8 vs. \u03C3 (\u03BA fixed)",
      scene = list(
        xaxis = list(title = "\u03B8", titlefont = bigfont),
        yaxis = list(title = "\u03C3", titlefont = bigfont),
        zaxis = list(title = "\u2113", titlefont = bigfont)
      )
    )
  
  fig1
  fig2
  fig3
}
print(fig1)
print(fig2)
print(fig3)
