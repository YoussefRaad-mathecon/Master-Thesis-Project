####################################################################################################
####################################################################################################
#------------------------------------------ seed and parameters -----------------------------------#
####################################################################################################
####################################################################################################

set.seed(112233)


####################################################################################################
#------------------------------------------   theta  ----------------------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N, dt = 1 / 252) {
  # Parameter extraction (explicit)
  kappa <- exp(theta.star[1])                                # scalar
  theta <- exp(theta.star[2:(1 + N)])                        # state-dependent
  sigma <- exp(theta.star[2 + N])                            # scalar
  off_diag_elements <- exp(theta.star[(3 + N):length(theta.star)])  # remaining parameters
  
  # Transition matrix Gamma
  Gamma <- diag(N)
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  # Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    for (i in 2:length(ind)) {
      q <- 2 * kappa * theta[j] / sigma^2 - 1
      c <- 2 * kappa / ((1 - exp(-kappa * dt)) * sigma^2)
      u <- c * x[ind[i - 1]] * exp(-kappa * dt)
      v <- c * x[ind[i]]
      
      log_prob <- log(c) - u - v + q / 2 * (log(v) - log(u)) +
        log(besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE)) +
        2 * sqrt(u * v)
      
      allprobs[ind[i], j] <- exp(log_prob)
    }
  }
  
  # Forward algorithm
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)) {
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
mllk <- function(theta.star, x, N, dt = 1/252) {
  kappa <- exp(theta.star[1:N])  # Mean reversion speed
  theta <- exp(theta.star[(N + 1)])  # Long-term mean level
  sigma <- exp(theta.star[N + 2])  # Volatility coefficient
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(N + 3):(length(theta.star))])
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
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j] * dt))*sigma[1]^2)
      u <- c * x[ind[i-1]] * exp(-kappa[j] * dt)
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
mllk <- function(theta.star, x, N, dt = 1/252) {
  kappa <- exp(theta.star[1])                      # scalar
  theta <- exp(theta.star[2])                      # scalar
  sigma <- exp(theta.star[3:(2 + N)])              # state-dependent volatility
  off_diag_elements <- exp(theta.star[(3 + N):length(theta.star)])
  
  # Transition matrix
  Gamma <- diag(N)
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  # Likelihood
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    for (i in 2:length(ind)) {
      q <- 2 * kappa * theta / sigma[j]^2 - 1
      c <- 2 * kappa / ((1 - exp(-kappa * dt)) * sigma[j]^2)
      u <- c * x[ind[i - 1]] * exp(-kappa * dt)
      v <- c * x[ind[i]]
      
      log_prob <- log(c) - u - v + q / 2 * (log(v) - log(u)) +
        log(besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE)) +
        2 * sqrt(u * v)
      
      allprobs[ind[i], j] <- exp(log_prob)
    }
  }
  
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)) {
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
#------------------------------------- sigma & theta state-dependent ------------------------------#
####################################################################################################
mllk <- function(theta.star, x, N, dt = 1/252) {
  # Parameter extraction
  kappa <- exp(theta.star[1])                                # scalar
  theta <- exp(theta.star[2:(1 + N)])                        # state-dependent theta
  sigma <- exp(theta.star[(2 + N):(1 + 2 * N)])              # state-dependent sigma
  off_diag_elements <- exp(theta.star[(2 + 2 * N):length(theta.star)])  # Gamma off-diagonal
  
  # Transition matrix Gamma
  Gamma <- diag(N)
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  # Likelihood setup
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    for (i in 2:length(ind)) {
      q <- 2 * kappa * theta[j] / sigma[j]^2 - 1
      c <- 2 * kappa / ((1 - exp(-kappa * dt)) * sigma[j]^2)
      u <- c * x[ind[i - 1]] * exp(-kappa * dt)
      v <- c * x[ind[i]]
      
      log_prob <- log(c) - u - v + q / 2 * (log(v) - log(u)) +
        log(besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE)) +
        2 * sqrt(u * v)
      
      allprobs[ind[i], j] <- exp(log_prob)
    }
  }
  
  # Forward algorithm
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)) {
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
mllk <- function(theta.star, x, N, dt = 1/252) {
  # Parameter extraction
  kappa <- exp(theta.star[1:N])                             # state-dependent
  theta <- exp(theta.star[(N + 1):(2 * N)])                 # state-dependent
  sigma <- exp(theta.star[2 * N + 1])                       # scalar
  off_diag_elements <- exp(theta.star[(2 * N + 2):length(theta.star)])
  
  # Transition matrix Gamma
  Gamma <- diag(N)
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  # Likelihood computation
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    for (i in 2:length(ind)) {
      q <- 2 * kappa[j] * theta[j] / sigma^2 - 1
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j] * dt)) * sigma^2)
      u <- c * x[ind[i - 1]] * exp(-kappa[j] * dt)
      v <- c * x[ind[i]]
      
      log_prob <- log(c) - u - v + q / 2 * (log(v) - log(u)) +
        log(besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE)) +
        2 * sqrt(u * v)
      
      allprobs[ind[i], j] <- exp(log_prob)
    }
  }
  
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)) {
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
mllk <- function(theta.star, x, N, dt = 1/252) {
  # Parameter extraction
  kappa <- exp(theta.star[1:N])                              # state-dependent
  theta <- exp(theta.star[N + 1])                            # scalar
  sigma <- exp(theta.star[(N + 2):(2 * N + 1)])              # state-dependent
  off_diag_elements <- exp(theta.star[(2 * N + 2):length(theta.star)])
  
  # Transition matrix Gamma
  Gamma <- diag(N)
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  # Likelihood computation
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    for (i in 2:length(ind)) {
      q <- 2 * kappa[j] * theta / sigma[j]^2 - 1
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j] * dt)) * sigma[j]^2)
      u <- c * x[ind[i - 1]] * exp(-kappa[j] * dt)
      v <- c * x[ind[i]]
      
      log_prob <- log(c) - u - v + q / 2 * (log(v) - log(u)) +
        log(besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE)) +
        2 * sqrt(u * v)
      
      allprobs[ind[i], j] <- exp(log_prob)
    }
  }
  
  # Forward algorithm
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)) {
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
mllk <- function(theta.star, x, N, dt = 1 /252) {
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
      c <- 2 * kappa[j] / ((1 - exp(-kappa[j] * dt ))*sigma[j]^2)
      u <- c * x[ind[i-1]] * exp(-kappa[j]* dt)
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
#theta.star <- generate_theta_star(N)
theta.star <- c(log(c(0.055, 0.027, 0.001, 0.001, 0.010)),     # kappa
                log(c(0.026, 0.1228, 0.651, 0.247, 0.125)),     # theta
                log(c(0.063, 0.079, 0.021, 0.016, 0.032)),     # sigma
                rep(-5, (5 - 1) * 5))


mod5 <- optim(
  par = theta.star, 
  fn = mllk, 
  method = "Nelder-Mead", 
  x = as.numeric(yields_df$"3M"), 
  N = N, 
  control = list(maxit = 100000, trace = 2)
)


N <- 5

#19067 function evaluations used
#lots of In besselI(x = 2 * sqrt(u * v), nu = q, expon.scaled = TRUE) :
#bessel_i(565.684,nu=1439): præcision tabt i resultat
Gamma5_hat <- diag(N)
Gamma5_hat[!Gamma5_hat] <- exp(mod5$par[(3 * N + 1):(length(theta.star))])
Gamma5_hat <- Gamma5_hat / rowSums(Gamma5_hat)
Gamma5_hat 

#[,1]        [,2]        [,3]        [,4]        [,5]
#[1,] 0.968723346 0.006821676 0.006821676 0.006821676 0.010811627
#[2,] 0.006849003 0.972603989 0.006849003 0.006849003 0.006849003
#[3,] 0.006849003 0.006849003 0.972603989 0.006849003 0.006849003
#[4,] 0.006878083 0.002632234 0.006878083 0.976733518 0.006878083
#[5,] 0.006849003 0.006849003 0.006849003 0.006849003 0.972603989

delta5_hat <- solve(t(diag(N) - Gamma5_hat + 1), rep(1, N)) ### 0.1775503 0.8224497s
delta5_hat #0.13822867 0.02051568 0.24165131 0.38129939 0.21830494

kappa5_hat <- exp(mod5$par[1:N])
theta5_hat <- exp(mod5$par[(N + 1):(2 * N)]) 
sigma5_hat <- exp(mod5$par[(2 * N + 1):(3 * N)])
kappa5_hat ### 0.0057481275 0.0028218081 0.0001045114 0.0001045114 0.0010451141
theta5_hat ### 0.002717297 0.128340011 0.068036928 0.025814318 0.013063926
sigma5_hat ### 0.006372595 0.028564435 0.002122917 0.001592439 0.003238651
##loglikelihood: ### -64246.96
mod5$value

load("fitted_models_HMM.RData")

