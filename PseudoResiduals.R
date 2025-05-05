####################################################################################################################
####################################################################################################################
#----------------------------------------- Pseudo Residuals --------------------------------------------------------
####################################################################################################################
####################################################################################################################
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis")
load("fitted_models_HMM.RData")


#----------------------------------------- N-state HMM --------------------------------------------------------


### Extract CIR-ARHMM parameters
extract_parameters_CIR <- function(theta.star, N) {
  kappa <- exp(theta.star[1:N])
  theta <- exp(theta.star[(N + 1):(2 * N)])
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)]) 
  
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[(3 * N + 1):length(theta.star)])
  Gamma <- Gamma / rowSums(Gamma)
  
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  
  list(kappa = kappa, theta = theta, sigma = sigma, Gamma = Gamma, delta = delta)
}


### Log-forward algorithm for CIR-ARHMM
lForward_CIR <- function(y, mod, N, dt = 1) {
  params <- extract_parameters_CIR(mod$par, N)
  T <- length(y)
  
  kappa <- params$kappa
  theta <- params$theta
  sigma <- params$sigma
  delta <- params$delta
  Gamma <- params$Gamma
  
  q <- 2 * kappa * theta / sigma^2 - 1
  c <- 2 * kappa / ((1 - exp(-kappa * dt)) * sigma^2)
  
  lalpha <- matrix(NA, N, T)
  lalpha[, 1] <- log(delta)  # log-initialize with stationary distribution
  
  lscale <- 0  # running log-scale tracker for stability (optional)
  
  for (i in 2:T) {
    u <- c * y[i - 1] * exp(-kappa * dt)
    df <- 2 * q + 2
    ncp <- 2 * u
    P <- dchisq(y[i] * (2 * c), df = df, ncp = ncp) * (2 * c)
    
    # Log-sum-exp for numerical stability
    cmax <- max(lalpha[, i - 1])
    foo <- exp(lalpha[, i - 1] - cmax)
    foo <- foo / sum(foo)  # normalize
    foo <- foo %*% Gamma * P
    
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}


### Compute pseudo-residuals for CIR-ARHMM
PseudoResiduals <- function(y, mod, N, dt = 1) {
  params <- extract_parameters_CIR(mod$par, N)
  kappa <- params$kappa
  theta <- params$theta
  sigma <- params$sigma
  delta <- params$delta
  Gamma <- params$Gamma
  
  T <- length(y)
  q <- 2 * kappa * theta / sigma^2 - 1
  c <- 2 * kappa / ((1 - exp(-kappa * dt)) * sigma^2)
  
  la <- t(lForward_CIR(y = y, mod = mod, N = N, dt = dt))
  Res <- rep(NA, T - 1)
  pMat <- matrix(NA, nrow = T - 1, ncol = N)
  
  # Residual at time t = 2 (based on y[1])
  u <- c * y[1] * exp(-kappa * dt)
  df <- 2 * q + 2
  ncp <- 2 * u
  pMat[1, ] <- pchisq(y[2] * (2 * c), df = df, ncp = ncp)
  Res[1] <- qnorm(drop(params$delta %*% pMat[1, ]))
  
  # Residuals for t = 3 to T
  # Residuals for t = 3 to T
  for (i in 3:T) {
    if (i %% 100 == 0) {
      message("Computing residual ", i - 1, " of ", T - 1)
    }
    
    u <- c * y[i - 1] * exp(-kappa * dt)
    df <- 2 * q + 2
    ncp <- 2 * u
    pMat[i - 1, ] <- pchisq(y[i] * (2 * c), df = df, ncp = ncp)
    
    cmax <- max(la[i - 2, ])
    a <- exp(la[i - 2, ] - cmax)
    weights <- a / sum(a)
    predictive_probs <- drop(weights %*% Gamma)
    
    Res[i - 1] <- qnorm(predictive_probs %*% pMat[i - 1, ])
  }
  
  
  return(list(Res = Res, pMat = pMat))
}


#----------------------------------------- PseudoResiduals --------------------------------------------------------



### Missing Values
### theta
pseudo_res_2_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_theta, N = 2)
pseudo_res_3_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_theta, N = 3)
pseudo_res_4_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_theta, N = 4) # 4 errors
pseudo_res_5_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_theta, N = 5) # 2 errors
sum(!is.finite(pseudo_res_2_theta$Res)) # 0
sum(!is.finite(pseudo_res_3_theta$Res)) # 1
sum(!is.finite(pseudo_res_4_theta$Res)) # 0
sum(!is.finite(pseudo_res_5_theta$Res)) # 0

### kappa
pseudo_res_2_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa, N = 2)
pseudo_res_3_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa, N = 3)
pseudo_res_4_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa, N = 4)
pseudo_res_5_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa, N = 5)
sum(!is.finite(pseudo_res_2_kappa$Res)) # 0
sum(!is.finite(pseudo_res_3_kappa$Res)) # 0 
sum(!is.finite(pseudo_res_4_kappa$Res)) # 0
sum(!is.finite(pseudo_res_5_kappa$Res)) # 0

### sigma
pseudo_res_2_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_sigma, N = 2)
pseudo_res_3_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_sigma, N = 3)
pseudo_res_4_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_sigma, N = 4)
pseudo_res_5_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_sigma, N = 5)
sum(!is.finite(pseudo_res_2_sigma$Res)) # 3435
sum(!is.finite(pseudo_res_3_sigma$Res)) # 3443
sum(!is.finite(pseudo_res_4_sigma$Res)) # 0
sum(!is.finite(pseudo_res_5_sigma$Res)) # 3437


### sigma theta
pseudo_res_2_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_sigma_theta, N = 2)
pseudo_res_3_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_sigma_theta, N = 3)
pseudo_res_4_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_sigma_theta, N = 4)
pseudo_res_5_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_sigma_theta, N = 5)
sum(!is.finite(pseudo_res_2_sigma_theta$Res)) # 3437
sum(!is.finite(pseudo_res_3_sigma_theta$Res)) # 0
sum(!is.finite(pseudo_res_4_sigma_theta$Res)) # 0
sum(!is.finite(pseudo_res_5_sigma_theta$Res)) # 0




### sigma kappa
pseudo_res_2_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa_sigma, N = 2)
pseudo_res_3_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa_sigma, N = 3)
pseudo_res_4_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa_sigma, N = 4)
pseudo_res_5_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa_sigma, N = 5)
sum(!is.finite(pseudo_res_2_kappa_sigma$Res)) # 3437
sum(!is.finite(pseudo_res_3_kappa_sigma$Res)) # 0
sum(!is.finite(pseudo_res_4_kappa_sigma$Res)) # 0
sum(!is.finite(pseudo_res_5_kappa_sigma$Res)) # 0




### kappa theta
pseudo_res_2_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa_theta, N = 2)
pseudo_res_3_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa_theta, N = 3)
pseudo_res_4_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa_theta, N = 4) # 3 errors
pseudo_res_5_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa_theta, N = 5)
sum(!is.finite(pseudo_res_2_kappa_theta$Res)) # 9
sum(!is.finite(pseudo_res_3_kappa_theta$Res)) # 0
sum(!is.finite(pseudo_res_4_kappa_theta$Res)) # 3439
sum(!is.finite(pseudo_res_5_kappa_theta$Res)) # 0






### sigma theta kappa
pseudo_res_2_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2, N = 2)
pseudo_res_3_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3, N = 3)
pseudo_res_4_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4, N = 4)
pseudo_res_5_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5, N = 5)
sum(!is.finite(pseudo_res_2_sigma_kappa_theta$Res)) # 0
sum(!is.finite(pseudo_res_3_sigma_kappa_theta$Res)) # 0
sum(!is.finite(pseudo_res_4_sigma_kappa_theta$Res)) # 0
sum(!is.finite(pseudo_res_5_sigma_kappa_theta$Res)) # 0


# Load
load("pseudoresiduals.RData")

# Collect all pseudo-residuals into a list
pseudo_residuals_list <- list(
  pseudo_res_2_theta = pseudo_res_2_theta,
  pseudo_res_3_theta = pseudo_res_3_theta,
  pseudo_res_4_theta = pseudo_res_4_theta,
  pseudo_res_5_theta = pseudo_res_5_theta,
  
  pseudo_res_2_kappa = pseudo_res_2_kappa,
  pseudo_res_3_kappa = pseudo_res_3_kappa,
  pseudo_res_4_kappa = pseudo_res_4_kappa,
  pseudo_res_5_kappa = pseudo_res_5_kappa,
  
  pseudo_res_2_sigma = pseudo_res_2_sigma,
  pseudo_res_3_sigma = pseudo_res_3_sigma,
  pseudo_res_4_sigma = pseudo_res_4_sigma,
  pseudo_res_5_sigma = pseudo_res_5_sigma,
  
  pseudo_res_2_sigma_theta = pseudo_res_2_sigma_theta,
  pseudo_res_3_sigma_theta = pseudo_res_3_sigma_theta,
  pseudo_res_4_sigma_theta = pseudo_res_4_sigma_theta,
  pseudo_res_5_sigma_theta = pseudo_res_5_sigma_theta,
  
  pseudo_res_2_kappa_sigma = pseudo_res_2_kappa_sigma,
  pseudo_res_3_kappa_sigma = pseudo_res_3_kappa_sigma,
  pseudo_res_4_kappa_sigma = pseudo_res_4_kappa_sigma,
  pseudo_res_5_kappa_sigma = pseudo_res_5_kappa_sigma,
  
  pseudo_res_2_kappa_theta = pseudo_res_2_kappa_theta,
  pseudo_res_3_kappa_theta = pseudo_res_3_kappa_theta,
  pseudo_res_4_kappa_theta = pseudo_res_4_kappa_theta,
  pseudo_res_5_kappa_theta = pseudo_res_5_kappa_theta,
  
  pseudo_res_2_sigma_kappa_theta = pseudo_res_2_sigma_kappa_theta,
  pseudo_res_3_sigma_kappa_theta = pseudo_res_3_sigma_kappa_theta,
  pseudo_res_4_sigma_kappa_theta = pseudo_res_4_sigma_kappa_theta,
  pseudo_res_5_sigma_kappa_theta = pseudo_res_5_sigma_kappa_theta
)


#----------------------------------------- QQ-plot --------------------------------------------------------


png("big_plot.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(7, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))



### theta
qqnorm(pseudo_res_2_theta$Res[is.finite(pseudo_res_2_theta$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_theta$Res[is.finite(pseudo_res_2_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_theta$Res[is.finite(pseudo_res_3_theta$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_theta$Res[is.finite(pseudo_res_3_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_theta$Res[is.finite(pseudo_res_4_theta$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_theta$Res[is.finite(pseudo_res_4_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_theta$Res[is.finite(pseudo_res_5_theta$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_theta$Res[is.finite(pseudo_res_5_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

### kappa
qqnorm(pseudo_res_2_kappa$Res[is.finite(pseudo_res_2_kappa$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_kappa$Res[is.finite(pseudo_res_2_kappa$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 1", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_kappa$Res[is.finite(pseudo_res_3_kappa$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_kappa$Res[is.finite(pseudo_res_3_kappa$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_kappa$Res[is.finite(pseudo_res_4_kappa$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_kappa$Res[is.finite(pseudo_res_4_kappa$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_kappa$Res[is.finite(pseudo_res_5_kappa$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_kappa$Res[is.finite(pseudo_res_5_kappa$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

### sigma
qqnorm(pseudo_res_2_sigma$Res[is.finite(pseudo_res_2_sigma$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_sigma$Res[is.finite(pseudo_res_2_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 3436", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_sigma$Res[is.finite(pseudo_res_3_sigma$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_sigma$Res[is.finite(pseudo_res_3_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 3437", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_sigma$Res[is.finite(pseudo_res_4_sigma$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_sigma$Res[is.finite(pseudo_res_4_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_sigma$Res[is.finite(pseudo_res_5_sigma$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_sigma$Res[is.finite(pseudo_res_5_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 3437", bty = "n", cex = 0.8)

### theta sigma
qqnorm(pseudo_res_2_sigma_theta$Res[is.finite(pseudo_res_2_sigma_theta$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_sigma_theta$Res[is.finite(pseudo_res_2_sigma_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 3437", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_sigma_theta$Res[is.finite(pseudo_res_3_sigma_theta$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_sigma_theta$Res[is.finite(pseudo_res_3_sigma_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_sigma_theta$Res[is.finite(pseudo_res_4_sigma_theta$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_sigma_theta$Res[is.finite(pseudo_res_4_sigma_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_sigma_theta$Res[is.finite(pseudo_res_5_sigma_theta$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_sigma_theta$Res[is.finite(pseudo_res_5_sigma_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

### kappa theta
qqnorm(pseudo_res_2_kappa_theta$Res[is.finite(pseudo_res_2_kappa_theta$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_kappa_theta$Res[is.finite(pseudo_res_2_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_kappa_theta$Res[is.finite(pseudo_res_3_kappa_theta$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_kappa_theta$Res[is.finite(pseudo_res_3_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_kappa_theta$Res[is.finite(pseudo_res_4_kappa_theta$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_kappa_theta$Res[is.finite(pseudo_res_4_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 3437", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_kappa_theta$Res[is.finite(pseudo_res_5_kappa_theta$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_kappa_theta$Res[is.finite(pseudo_res_5_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

### sigma kappa
qqnorm(pseudo_res_2_kappa_sigma$Res[is.finite(pseudo_res_2_kappa_sigma$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_kappa_sigma$Res[is.finite(pseudo_res_2_kappa_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 3437", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_kappa_sigma$Res[is.finite(pseudo_res_3_kappa_sigma$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_kappa_sigma$Res[is.finite(pseudo_res_3_kappa_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_kappa_sigma$Res[is.finite(pseudo_res_4_kappa_sigma$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_kappa_sigma$Res[is.finite(pseudo_res_4_kappa_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 2", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_kappa_sigma$Res[is.finite(pseudo_res_5_kappa_sigma$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_kappa_sigma$Res[is.finite(pseudo_res_5_kappa_sigma$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

### theta sigma kappa
qqnorm(pseudo_res_2_sigma_kappa_theta$Res[is.finite(pseudo_res_2_sigma_kappa_theta$Res)], main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_sigma_kappa_theta$Res[is.finite(pseudo_res_2_sigma_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 2", bty = "n", cex = 0.8)

qqnorm(pseudo_res_3_sigma_kappa_theta$Res[is.finite(pseudo_res_3_sigma_kappa_theta$Res)], main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_sigma_kappa_theta$Res[is.finite(pseudo_res_3_sigma_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 1", bty = "n", cex = 0.8)

qqnorm(pseudo_res_4_sigma_kappa_theta$Res[is.finite(pseudo_res_4_sigma_kappa_theta$Res)], main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_sigma_kappa_theta$Res[is.finite(pseudo_res_4_sigma_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 0", bty = "n", cex = 0.8)

qqnorm(pseudo_res_5_sigma_kappa_theta$Res[is.finite(pseudo_res_5_sigma_kappa_theta$Res)], main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_sigma_kappa_theta$Res[is.finite(pseudo_res_5_sigma_kappa_theta$Res)], col = "steelblue", lwd = 3)
legend("bottomright", legend = "Missing values: 1", bty = "n", cex = 0.8)



sum(pseudo_res_5_sigma_kappa_theta$Res >= -2.5 & pseudo_res_5_sigma_kappa_theta$Res <= 2.5, na.rm = TRUE)

9479/9616

mtext("Sample Quantiles", side = 2, line = 1.2, outer = TRUE, cex = 1.4)


# Properly centered row labels
at_vals <- rev(((1:7) - 0.5) / 7)

mtext(expression(theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[1])
mtext(expression(kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[2])
mtext(expression(sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[3])
mtext(expression(theta * ", " * sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[4])
mtext(expression(kappa * ", " * theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[5])
mtext(expression(sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[6])
mtext(expression(theta * ", " * sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[7])
dev.off()


#----------------------------------------------Time Series plot------------------------------------------#

# Set up PNG output
png("TimeSeries.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(7, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))

# Aligned dates
date_aligned <- yields_df$DateCont[-1]

# Loop through residuals
for (name in names(pseudo_residuals_list)) {
  res <- pseudo_residuals_list[[name]]$Res
  
  # Skip if length doesn't match
  if (length(res) != length(date_aligned)) {
    warning(paste("Skipping", name, "due to length mismatch"))
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # blank plot
    next
  }
  
  # Remove NA values
  complete_idx <- complete.cases(res)
  dates_clean <- date_aligned[complete_idx]
  res_clean <- res[complete_idx]
  
  # Plot residual time series (with axis ticks, no title)
  plot(dates_clean, res_clean, type = "l", col = "#901a1E", lwd = 1,
       xlab = "", ylab = "")
}

# Global y-axis label
mtext("Residuals", side = 2, line = 1.2, outer = TRUE, cex = 1.4)

# Add labels for each row group on the right side
at_vals <- rev(((1:7) - 0.5) / 7)
mtext(expression(theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[1])
mtext(expression(kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[2])
mtext(expression(sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[3])
mtext(expression(theta * ", " * sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[4])
mtext(expression(kappa * ", " * theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[5])
mtext(expression(sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[6])
mtext(expression(theta * ", " * sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[7])

# Close device
dev.off()

#----------------------------------------------ACF plot--------------------------------------------------#
# Set up PNG for ACF plots
png("ACFSeries.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(7, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))

# Custom red color
acf_col <- "#901a1E"

for (name in names(pseudo_residuals_list)) {
  res <- pseudo_residuals_list[[name]]$Res
  
  # Check for valid residuals: at least 2 finite values
  if (is.null(res) || sum(is.finite(res)) < 2) {
    warning(paste("Skipping", name, "due to insufficient finite values"))
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # blank plot
    next
  }
  
  # Use only finite values
  res_clean <- res[is.finite(res)]
  
  # Compute ACF but don't plot yet
  acf_obj <- acf(res_clean, plot = FALSE)
  
  # Now plot manually with custom color
  plot(acf_obj, main = "", xlab = "", ylab = "", col = acf_col, lwd = 2)
}

# Global y-axis label
mtext("Autocorrelation", side = 2, line = 1.2, outer = TRUE, cex = 1.4)

# Right-side row labels
at_vals <- rev(((1:7) - 0.5) / 7)
mtext(expression(theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[1])
mtext(expression(kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[2])
mtext(expression(sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[3])
mtext(expression(theta * ", " * sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[4])
mtext(expression(kappa * ", " * theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[5])
mtext(expression(sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[6])
mtext(expression(theta * ", " * sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[7])

# Close PNG device
dev.off()


#------------------------------------------Histogram Plot--------------------------------------------------#
# Set up PNG for histograms
png("HistogramSeries.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(7, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))

# Custom red color for histogram bars
hist_col <- "#901a1E"

for (name in names(pseudo_residuals_list)) {
  res <- pseudo_residuals_list[[name]]$Res
  
  # Check for valid residuals: at least 2 finite values
  if (is.null(res) || sum(is.finite(res)) < 2) {
    warning(paste("Skipping", name, "due to insufficient finite values"))
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # blank plot
    next
  }
  
  # Use only finite values
  res_clean <- res[is.finite(res)]
  
  # Trim to 1st–99th percentile range
  q_low <- quantile(res_clean, 0.01, na.rm = TRUE)
  q_high <- quantile(res_clean, 0.99, na.rm = TRUE)
  res_trimmed <- res_clean[res_clean >= q_low & res_clean <= q_high]
  
  # Plot histogram
  hist(res_trimmed,
       breaks = 20,
       col = hist_col,
       border = "white",
       main = "",
       xlab = "",
       ylab = "")
}

# Global y-axis label
mtext("Frequency", side = 2, line = 1.2, outer = TRUE, cex = 1.4)

# Right-side row labels
at_vals <- rev(((1:7) - 0.5) / 7)
mtext(expression(theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[1])
mtext(expression(kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[2])
mtext(expression(sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[3])
mtext(expression(theta * ", " * sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[4])
mtext(expression(kappa * ", " * theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[5])
mtext(expression(sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[6])
mtext(expression(theta * ", " * sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[7])

# Close PNG device
dev.off()


#------------------------------------------Histogram Plot (All Residuals)------------------------------------------#
# Set up PNG for histograms
png("HistogramSeries_AllResiduals.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(7, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))

# Custom red color for histogram bars
hist_col <- "#901a1E"

for (name in names(pseudo_residuals_list)) {
  res <- pseudo_residuals_list[[name]]$Res
  
  # Check for valid residuals: at least 2 finite values
  if (is.null(res) || sum(is.finite(res)) < 2) {
    warning(paste("Skipping", name, "due to insufficient finite values"))
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # blank plot
    next
  }
  
  # Use all finite residuals
  res_clean <- res[is.finite(res)]
  
  # Plot histogram
  hist(res_clean,
       breaks = 20,
       col = hist_col,
       border = "white",
       main = "",
       xlab = "",
       ylab = "")
}

# Global y-axis label
mtext("Frequency", side = 2, line = 1.2, outer = TRUE, cex = 1.4)

# Right-side row labels
at_vals <- rev(((1:7) - 0.5) / 7)
mtext(expression(theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[1])
mtext(expression(kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[2])
mtext(expression(sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[3])
mtext(expression(theta * ", " * sigma), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[4])
mtext(expression(kappa * ", " * theta), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[5])
mtext(expression(sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[6])
mtext(expression(theta * ", " * sigma * ", " * kappa), side = 4, line = 1, outer = TRUE, cex = 2, at = at_vals[7])

# Close PNG device
dev.off()

