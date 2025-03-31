####################################################################################################
####################################################################################################
#------------------------------------------ seed and parameters -----------------------------------#
####################################################################################################
####################################################################################################

set.seed(112233)


####################################################################################################
#------------------------------------------   theta  ----------------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N, dt = 1) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[1] * theta[j] / sigma[1]^2 - 1
      c <- 2 * kappa[1] / ((1 - exp(-kappa[1] * dt))*sigma[1]^2)
      u <- c * x[ind[i-1]] * exp(-kappa[1] * dt)
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

### Generate initial values function that works dynamically for N
generate_theta_star <- function(N) {
  c(log(seq(0.001, 0.2, length = N)),     # kappa
    log(seq(0.001, 0.08, length = N)),    # theta
    log(seq(0.004, 0.01, length = N)),    # sigma
    rep(-5, (N - 1) * N))                 # Gamma off-diagonal elements
}

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N)
mod2_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)


#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#----------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)



####################################################################################################
#------------------------------------------   kappa  ----------------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[j] * theta[1] / sigma[1]^2 - 1
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j]))*sigma[1]^2)
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

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N) 
mod2_kappa <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3_kappa <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4_kappa <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#----------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5_kappa <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)







####################################################################################################
#------------------------------------------   sigma  ----------------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[1] * theta[1] / sigma[j]^2 - 1
      c <- 2 * kappa[1] / ((1 - exp(-kappa[1]))*sigma[j]^2)
      u <- c * x[ind[i-1]] * exp(-kappa[1])
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

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N) 
mod2_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#----------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)





####################################################################################################
#------------------------------------------ sigma & theta -----------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[1] * theta[j] / sigma[j]^2 - 1
      c <- 2 * kappa[1] / ((1 - exp(-kappa[1]))*sigma[j]^2)
      u <- c * x[ind[i-1]] * exp(-kappa[1])
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

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N) 
mod2_sigma_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3_sigma_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4_sigma_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#----------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5_sigma_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)



####################################################################################################
#------------------------------------------ kappa & theta -----------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[j] * theta[j] / sigma[1]^2 - 1
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j]))*sigma[1]^2)
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

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N) 
mod2_kappa_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3_kappa_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4_kappa_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#----------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5_kappa_theta <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)




####################################################################################################
#------------------------------------------ sigma & kappa -----------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    # CIR density calculation for each state
    for (i in 2:length(ind)) {  # Start from the second observation
      q <- 2 * kappa[j] * theta[1] / sigma[j]^2 - 1
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

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N) 
mod2_kappa_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3_kappa_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4_kappa_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)
#----------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5_kappa_sigma <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)





####################################################################################################
#------------------------------------------ sigma & kappa & theta ---------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1):(2 * N)])  # Long-term mean level
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  
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

#------------------------------------2 state----------------------------------------------------#
### 2-state HMM
N <- 2
theta.star <- generate_theta_star(N) 
mod2 <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)


#------------------------------------3 state----------------------------------------------------#
### 3-state HMM
N <- 3
theta.star <- generate_theta_star(N) 
mod3 <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)



#------------------------------------4 state----------------------------------------------------#
### 4-state HMM
N <- 4
theta.star <- generate_theta_star(N) 
mod4 <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)


#------------------------------------5 state----------------------------------------------------#
### 5-state HMM
N <- 5
theta.star <- generate_theta_star(N) 
mod5 <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)

#19067 function evaluations used
#lots of In besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE) :
#bessel_i(565.684,nu=1439): præcision tabt i resultat
Gamma5_hat <- diag(N)
Gamma5_hat[!Gamma5_hat] <- exp(mod5$par[(3 * N + 1):(length(theta.star))])
Gamma5_hat <- Gamma5_hat / rowSums(Gamma5_hat)
Gamma5_hat #[1,] 0.910089013 0.011195376 0.003562239 0.021939568 0.053213804
# [2,] 0.196296521 0.750352177 0.009416788 0.038750992 0.005183522
# [3,] 0.002462499 0.007568432 0.985645473 0.001989003 0.002334593
# [4,] 0.011177746 0.003365324 0.002830302 0.971760997 0.010865631
# [5,] 0.016234105 0.002116563 0.007805607 0.029587920 0.944255805

delta5_hat <- solve(t(diag(N) - Gamma5_hat + 1), rep(1, N)) ### 0.1775503 0.8224497
delta5_hat #0.13822867 0.02051568 0.24165131 0.38129939 0.21830494

kappa5_hat <- exp(mod5$par[1:N])
theta5_hat <- exp(mod5$par[(N + 1):(2 * N)]) 
sigma5_hat <- exp(mod5$par[(2 * N + 1):(3 * N)])
kappa5_hat ### 0.0054995529 0.0026673149 0.0009752358 0.0001473490 0.0010062147
theta5_hat ### 0.00255053 0.12278648 0.06509860 0.02471914 0.01251813
sigma5_hat ### 0.006278071 0.027900250 0.002108080 0.001602708 0.003243255
##loglikelihood: ### -64246.96
mod5$value




models_list <- list(
  mod2 = mod2_theta,
  mod3 = mod3_theta,
  mod4 = mod4_theta,
  mod5 = mod5_theta,
  
  mod2_kappa = mod2_kappa,
  mod3_kappa = mod3_kappa,
  mod4_kappa = mod4_kappa,
  mod5_kappa = mod5_kappa,
  
  mod2_sigma = mod2_sigma,
  mod3_sigma = mod3_sigma,
  mod4_sigma = mod4_sigma,
  mod5_sigma = mod5_sigma,
  
  mod2_kappa_sigma = mod2_kappa_sigma,
  mod3_kappa_sigma = mod3_kappa_sigma,
  mod4_kappa_sigma = mod4_kappa_sigma,
  mod5_kappa_sigma = mod5_kappa_sigma,
  
  mod2_sigma_theta = mod2_sigma_theta,
  mod3_sigma_theta = mod3_sigma_theta,
  mod4_sigma_theta = mod4_sigma_theta,
  mod5_sigma_theta = mod5_sigma_theta,
  
  mod2_kappa_theta = mod2_kappa_theta,
  mod3_kappa_theta = mod3_kappa_theta,
  mod4_kappa_theta = mod4_kappa_theta,
  mod5_kappa_theta = mod5_kappa_theta
)


save(models_list, file = "fitted_models.RData")
#load("fitted_models.RData")

