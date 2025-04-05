####################################################################################################################
####################################################################################################################
#------------------------------------------- AIC & BIC -------------------------------------------------------------
####################################################################################################################
####################################################################################################################


#-------------------------------- Finite State Space HMM -------------------------------------------------------------



### Recall AIC
### AIC = 2k − 2ln(L)
### Where k = # estimated parameters and L = maximized value of the likelihood function

### Recall BIC
### BIC = k ln(n) − 2ln(L)
### Where k = # estimated parameters and L = maximized value of the likelihood function and n = # observations


### Extract the negative log-likelihood value at the minimum
neg_loglik_2 <- mod2$minimum
neg_loglik_3 <- mod3$minimum
neg_loglik_4 <- mod4$minimum
neg_loglik_5 <- mod5$minimum

neg_loglik_2_theta <- mod2_theta$minimum
neg_loglik_3_theta <- mod3_theta$minimum
neg_loglik_4_theta <- mod4_theta$minimum
neg_loglik_5_theta <- mod5_theta$minimum

neg_loglik_2_kappa <- mod2_kappa$minimum
neg_loglik_3_kappa <- mod3_kappa$minimum
neg_loglik_4_kappa <- mod4_kappa$minimum
neg_loglik_5_kappa <- mod5_kappa$minimum

neg_loglik_2_sigma <- mod2_sigma$minimum
neg_loglik_3_sigma <- mod3_sigma$minimum
neg_loglik_4_sigma <- mod4_sigma$minimum
neg_loglik_5_sigma <- mod5_sigma$minimum

neg_loglik_2_sigma_theta <- mod2_sigma_theta$minimum
neg_loglik_3_sigma_theta <- mod3_sigma_theta$minimum
neg_loglik_4_sigma_theta <- mod4_sigma_theta$minimum
neg_loglik_5_sigma_theta <- mod5_sigma_theta$minimum

neg_loglik_2_kappa_theta <- mod2_kappa_theta$minimum
neg_loglik_3_kappa_theta <- mod3_kappa_theta$minimum
neg_loglik_4_kappa_theta <- mod4_kappa_theta$minimum
neg_loglik_5_kappa_theta <- mod5_kappa_theta$minimum

neg_loglik_2_kappa_sigma <- mod2_kappa_sigma$minimum
neg_loglik_3_kappa_sigma <- mod3_kappa_sigma$minimum
neg_loglik_4_kappa_sigma <- mod4_kappa_sigma$minimum
neg_loglik_5_kappa_sigma <- mod5_kappa_sigma$minimum


### Number of parameters estimated
num_params_2 <- length(mod2$par)
num_params_3 <- length(mod3$par)
num_params_4 <- length(mod4$par)
num_params_5 <- length(mod5$par)

num_params_2_theta <- length(mod2_theta$par) - 1 - 1
num_params_3_theta <- length(mod3_theta$par) - 2 - 2
num_params_4_theta <- length(mod4_theta$par) - 3 - 3
num_params_5_theta <- length(mod5_theta$par) - 4 - 4

num_params_2_kappa <- length(mod2_kappa$par) - 1 - 1
num_params_3_kappa <- length(mod3_kappa$par) - 2 - 2
num_params_4_kappa <- length(mod4_kappa$par) - 3 - 3
num_params_5_kappa <- length(mod5_kappa$par) - 4 - 4

num_params_2_sigma <- length(mod2_sigma$par) - 1 - 1
num_params_3_sigma <- length(mod3_sigma$par) - 2 - 2
num_params_4_sigma <- length(mod4_sigma$par) - 3 - 3
num_params_5_sigma <- length(mod5_sigma$par) - 4 - 4

num_params_2_sigma_theta <- length(mod2_sigma_theta$par) - 1
num_params_3_sigma_theta <- length(mod3_sigma_theta$par) - 2 
num_params_4_sigma_theta <- length(mod4_sigma_theta$par) - 3
num_params_5_sigma_theta <- length(mod5_sigma_theta$par) - 4 

num_params_2_kappa_theta <- length(mod2_kappa_theta$par) - 1 
num_params_3_kappa_theta <- length(mod3_kappa_theta$par) - 2 
num_params_4_kappa_theta <- length(mod4_kappa_theta$par) - 3 
num_params_5_kappa_theta <- length(mod5_kappa_theta$par) - 4

num_params_2_kappa_sigma <- length(mod2_kappa_sigma$par) - 1 
num_params_3_kappa_sigma <- length(mod3_kappa_sigma$par) - 2 
num_params_4_kappa_sigma <- length(mod4_kappa_sigma$par) - 3 
num_params_5_kappa_sigma <- length(mod5_kappa_sigma$par) - 4 


### Number of observations
num_obs <- length(yields_df$"3M")

### Calculate AIC
aic_2 <- 2 * num_params_2 + 2 * neg_loglik_2
aic_3 <- 2 * num_params_3 + 2 * neg_loglik_3
aic_4 <- 2 * num_params_4 + 2 * neg_loglik_4
aic_5 <- 2 * num_params_5 + 2 * neg_loglik_5

aic_2_theta <- 2 * num_params_2_theta + 2 * neg_loglik_2_theta
aic_3_theta <- 2 * num_params_3_theta + 2 * neg_loglik_3_theta
aic_4_theta <- 2 * num_params_4_theta + 2 * neg_loglik_4_theta
aic_5_theta <- 2 * num_params_5_theta + 2 * neg_loglik_5_theta

aic_2_kappa <- 2 * num_params_2_kappa + 2 * neg_loglik_2_kappa
aic_3_kappa <- 2 * num_params_3_kappa + 2 * neg_loglik_3_kappa
aic_4_kappa <- 2 * num_params_4_kappa + 2 * neg_loglik_4_kappa
aic_5_kappa <- 2 * num_params_5_kappa + 2 * neg_loglik_5_kappa

aic_2_sigma <- 2 * num_params_2_sigma + 2 * neg_loglik_2_sigma
aic_3_sigma <- 2 * num_params_3_sigma + 2 * neg_loglik_3_sigma
aic_4_sigma <- 2 * num_params_4_sigma + 2 * neg_loglik_4_sigma
aic_5_sigma <- 2 * num_params_5_sigma + 2 * neg_loglik_5_sigma

aic_2_sigma_theta <- 2 * num_params_2_sigma_theta + 2 * neg_loglik_2_sigma_theta
aic_3_sigma_theta <- 2 * num_params_3_sigma_theta + 2 * neg_loglik_3_sigma_theta
aic_4_sigma_theta <- 2 * num_params_4_sigma_theta + 2 * neg_loglik_4_sigma_theta
aic_5_sigma_theta <- 2 * num_params_5_sigma_theta + 2 * neg_loglik_5_sigma_theta

aic_2_kappa_theta <- 2 * num_params_2_kappa_theta + 2 * neg_loglik_2_kappa_theta
aic_3_kappa_theta <- 2 * num_params_3_kappa_theta + 2 * neg_loglik_3_kappa_theta
aic_4_kappa_theta <- 2 * num_params_4_kappa_theta + 2 * neg_loglik_4_kappa_theta
aic_5_kappa_theta <- 2 * num_params_5_kappa_theta + 2 * neg_loglik_5_kappa_theta

aic_2_kappa_sigma <- 2 * num_params_2_kappa_sigma + 2 * neg_loglik_2_kappa_sigma
aic_3_kappa_sigma <- 2 * num_params_3_kappa_sigma + 2 * neg_loglik_3_kappa_sigma
aic_4_kappa_sigma <- 2 * num_params_4_kappa_sigma + 2 * neg_loglik_4_kappa_sigma
aic_5_kappa_sigma <- 2 * num_params_5_kappa_sigma + 2 * neg_loglik_5_kappa_sigma


cbind(aic_2, aic_3, aic_4, aic_5)
cbind(aic_2_theta, aic_3_theta, aic_4_theta, aic_5_theta)
cbind(aic_2_kappa, aic_3_kappa, aic_4_kappa, aic_5_kappa)
cbind(aic_2_sigma, aic_3_sigma, aic_4_sigma, aic_5_sigma)
cbind(aic_2_sigma_theta, aic_3_sigma_theta, aic_4_sigma_theta, aic_5_sigma_theta)
cbind(aic_2_kappa_theta, aic_3_kappa_theta, aic_4_kappa_theta, aic_5_kappa_theta)
cbind(aic_2_kappa_sigma, aic_3_kappa_sigma, aic_4_kappa_sigma, aic_5_kappa_sigma)

### Calculate BIC
bic_2 <- num_params_2 * log(num_obs) + 2 * neg_loglik_2
bic_3 <- num_params_3 * log(num_obs) + 2 * neg_loglik_3
bic_4 <- num_params_4 * log(num_obs) + 2 * neg_loglik_4
bic_5 <- num_params_5 * log(num_obs) + 2 * neg_loglik_5

bic_2_theta <- num_params_2_theta * log(num_obs) + 2 * neg_loglik_2_theta
bic_3_theta <- num_params_3_theta * log(num_obs) + 2 * neg_loglik_3_theta
bic_4_theta <- num_params_4_theta * log(num_obs) + 2 * neg_loglik_4_theta
bic_5_theta <- num_params_5_theta * log(num_obs) + 2 * neg_loglik_5_theta

bic_2_kappa <- num_params_2_kappa * log(num_obs) + 2 * neg_loglik_2_kappa
bic_3_kappa <- num_params_3_kappa * log(num_obs) + 2 * neg_loglik_3_kappa
bic_4_kappa <- num_params_4_kappa * log(num_obs) + 2 * neg_loglik_4_kappa
bic_5_kappa <- num_params_5_kappa * log(num_obs) + 2 * neg_loglik_5_kappa

bic_2_sigma <- num_params_2_sigma * log(num_obs) + 2 * neg_loglik_2_sigma
bic_3_sigma <- num_params_3_sigma * log(num_obs) + 2 * neg_loglik_3_sigma
bic_4_sigma <- num_params_4_sigma * log(num_obs) + 2 * neg_loglik_4_sigma
bic_5_sigma <- num_params_5_sigma * log(num_obs) + 2 * neg_loglik_5_sigma

bic_2_sigma_theta <- num_params_2_sigma_theta * log(num_obs) + 2 * neg_loglik_2_sigma_theta
bic_3_sigma_theta <- num_params_3_sigma_theta * log(num_obs) + 2 * neg_loglik_3_sigma_theta
bic_4_sigma_theta <- num_params_4_sigma_theta * log(num_obs) + 2 * neg_loglik_4_sigma_theta
bic_5_sigma_theta <- num_params_5_sigma_theta * log(num_obs) + 2 * neg_loglik_5_sigma_theta

bic_2_kappa_theta <- num_params_2_kappa_theta * log(num_obs) + 2 * neg_loglik_2_kappa_theta
bic_3_kappa_theta <- num_params_3_kappa_theta * log(num_obs) + 2 * neg_loglik_3_kappa_theta
bic_4_kappa_theta <- num_params_4_kappa_theta * log(num_obs) + 2 * neg_loglik_4_kappa_theta
bic_5_kappa_theta <- num_params_5_kappa_theta * log(num_obs) + 2 * neg_loglik_5_kappa_theta

bic_2_kappa_sigma <- num_params_2_kappa_sigma * log(num_obs) + 2 * neg_loglik_2_kappa_sigma
bic_3_kappa_sigma <- num_params_3_kappa_sigma * log(num_obs) + 2 * neg_loglik_3_kappa_sigma
bic_4_kappa_sigma <- num_params_4_kappa_sigma * log(num_obs) + 2 * neg_loglik_4_kappa_sigma
bic_5_kappa_sigma <- num_params_5_kappa_sigma * log(num_obs) + 2 * neg_loglik_5_kappa_sigma

cbind(bic_2, bic_3, bic_4, bic_5)
cbind(bic_2_theta, bic_3_theta, bic_4_theta, bic_5_theta)
cbind(bic_2_kappa, bic_3_kappa, bic_4_kappa, bic_5_kappa)
cbind(bic_2_sigma, bic_3_sigma, bic_4_sigma, bic_5_sigma)
cbind(bic_2_sigma_theta, bic_3_sigma_theta, bic_4_sigma_theta, bic_5_sigma_theta)
cbind(bic_2_kappa_theta, bic_3_kappa_theta, bic_4_kappa_theta, bic_5_kappa_theta)
cbind(bic_2_kappa_sigma, bic_3_kappa_sigma, bic_4_kappa_sigma, bic_5_kappa_sigma)
