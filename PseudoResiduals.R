####################################################################################################################
####################################################################################################################
#----------------------------------------- Pseudo Residuals --------------------------------------------------------
####################################################################################################################
####################################################################################################################
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis")
load("fitted_models_HMM.RData")


#----------------------------------------- N-state HMM --------------------------------------------------------

extract_parameters_CIR <- function(theta.star, N) {
  kappa <- exp(theta.star[1:N]) * 252  # Scale kappa
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Theta (not scaled)
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)]) * sqrt(252)  # Scale sigma
  
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[(3 * N + 1):length(theta.star)])
  Gamma <- Gamma / rowSums(Gamma)
  
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  
  list(kappa = kappa, theta = theta, sigma = sigma, Gamma = Gamma, delta = delta)
}

lForward_CIR <- function(y, mod, N) {
  params <- extract_parameters_CIR(mod$par, N)
  T <- length(y)
  mus <- matrix(params$theta, nrow = T, ncol = N, byrow = TRUE)
  sigmas <- matrix(sqrt(params$sigma^2 / (2 * params$kappa)), nrow = T, ncol = N, byrow = TRUE)
  
  lalpha <- matrix(NA, N, T)
  P <- dnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  foo <- params$delta * P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  for (i in 2:T) {
    P <- dnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    foo <- foo %*% params$Gamma * P
    sumfoo <- sum(foo) 
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}

PseudoResiduals <- function(y, mod, N) {
  params <- extract_parameters_CIR(mod$par, N)
  mus <- matrix(params$theta, nrow = length(y), ncol = N, byrow = TRUE)
  sigmas <- matrix(sqrt(params$sigma^2 / (2 * params$kappa)), nrow = length(y), ncol = N, byrow = TRUE)
  
  la <- t(lForward_CIR(y = y, mod = mod, N = N))
  n <- length(y)
  Res <- rep(NA, n)
  pMat <- matrix(NA, nrow = n, ncol = N)
  
  pMat[1, ] <- pnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  Res[1] <- qnorm(params$delta %*% pMat[1, ])
  
  for (i in 2:n) {
    pMat[i, ] <- pnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    c <- max(la[i - 1, ])
    a <- exp(la[i - 1, ] - c)
    weighted_Gamma <- params$Gamma / sum(a)
    Res[i] <- qnorm(a %*% weighted_Gamma %*% pMat[i, ])
  }
  
  return(list(Res = Res))
}
#----------------------------------------- PseudoResiduals --------------------------------------------------------


### theta
pseudo_res_2_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_theta, N = 2)


pseudo_res_3_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_theta, N = 3)
pseudo_res_4_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_theta, N = 4)
pseudo_res_5_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_theta, N = 5)




### kappa
pseudo_res_2_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa, N = 2)
pseudo_res_3_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa, N = 3)
pseudo_res_4_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa, N = 4)
pseudo_res_5_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa, N = 5)





### sigma
pseudo_res_2_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_sigma, N = 2)
pseudo_res_3_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_sigma, N = 3)
pseudo_res_4_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_sigma, N = 4)
pseudo_res_5_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_sigma, N = 5)







### sigma theta
pseudo_res_2_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_sigma_theta, N = 2)
pseudo_res_3_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_sigma_theta, N = 3)
pseudo_res_4_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_sigma_theta, N = 4)
pseudo_res_5_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_sigma_theta, N = 5)





### sigma kappa
pseudo_res_2_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa_sigma, N = 2)
pseudo_res_3_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa_sigma, N = 3)
pseudo_res_4_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa_sigma, N = 4)
pseudo_res_5_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa_sigma, N = 5)





### kappa theta
pseudo_res_2_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa_theta, N = 2)
pseudo_res_3_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa_theta, N = 3)
pseudo_res_4_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa_theta, N = 4)
pseudo_res_5_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa_theta, N = 5)







### sigma theta kappa
pseudo_res_2_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2, N = 2)
pseudo_res_3_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3, N = 3)
pseudo_res_4_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4, N = 4)
pseudo_res_5_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5, N = 5)







#----------------------------------------- SMALL PLOTS --------------------------------------------------------


par(oma=c(3,3,0,1),mar=c(3,3,2,3),mfrow=c(3,4))

### theta
qqnorm(pseudo_res_2_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_theta$Res, col = "steelblue", lwd = 3)



### kappa
qqnorm(pseudo_res_2_kappa$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa$Res, col = "steelblue", lwd = 3)


### sigma
qqnorm(pseudo_res_2_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma$Res, col = "steelblue", lwd = 3)




mtext(text="Theoretical Quantiles",side=1,line=1,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2.7, adj = 0.5)
mtext(text=expression(theta),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.85)
mtext(text=expression(kappa),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.5)
mtext(text=expression(sigma),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.15)













par(oma=c(3,3,0,1),mar=c(3,3,2,3),mfrow=c(3,4))

### theta sigma
qqnorm(pseudo_res_2_sigma_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma_theta$Res, col = "steelblue", lwd = 3)



### kappa theta
qqnorm(pseudo_res_2_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa_theta$Res, col = "steelblue", lwd = 3)


### sigma kappa
qqnorm(pseudo_res_2_kappa_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa_sigma$Res, col = "steelblue", lwd = 3)




mtext(text="Theoretical Quantiles",side=1,line=1,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2.7, adj = 0.5)
mtext(text=expression(theta * ", " * sigma),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.85)
mtext(text=expression(theta * ", " * kappa),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.5)
mtext(text=expression(kappa * ", " * sigma),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.15)










par(oma=c(3,3,0,1),mar=c(3,3,2,3),mfrow=c(3,4))

### theta sigma kappa
qqnorm(pseudo_res_2_sigma_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

mtext(text="Theoretical Quantiles",side=1,line=1,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2.7, adj = 0.5)



#----------------------------------------- BIG PLOT --------------------------------------------------------


png("big_plot.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(7, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))



### theta
qqnorm(pseudo_res_2_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_theta$Res, col = "steelblue", lwd = 3)



### kappa
qqnorm(pseudo_res_2_kappa$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_kappa$Res, col = "steelblue", lwd = 3)


### sigma
qqnorm(pseudo_res_2_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_sigma$Res, col = "steelblue", lwd = 3)



### theta sigma
qqnorm(pseudo_res_2_sigma_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_sigma_theta$Res, col = "steelblue", lwd = 3)




### kappa theta
qqnorm(pseudo_res_2_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_kappa_theta$Res, col = "steelblue", lwd = 3)


### sigma kappa
qqnorm(pseudo_res_2_kappa_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_kappa_sigma$Res, col = "steelblue", lwd = 3)





### theta sigma kappa
qqnorm(pseudo_res_2_sigma_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_2_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_3_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_4_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 1.5)
qqline(pseudo_res_5_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)





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

