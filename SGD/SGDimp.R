
#####################################################
### PACKAGES
#####################################################

library(tidyverse)
library(CSwR)
library(reshape2)
library(tibble)
library(LaplacesDemon)
library(see)
library(knitr)
library(magrittr)
library(ggplot2)
library(latex2exp)
library(splines)
library(ggthemr)
library(patchwork)
library(bench)

#####################################################
### UTILITY FUNCTIONS
#####################################################

inverse_logit <- \(x) exp(x) / (1 + exp(x))

# Omega
pen_mat <- function(inner_knots) {
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots))
  d <- diff(inner_knots)  # The vector of knot differences; b - a
  g_ab <- splineDesign(knots, inner_knots, derivs = 2)
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) +
      4 * crossprod(d * g_ab_mid, g_ab_mid) +
      crossprod(d * g_b, g_b)) / 6
}  

decay_scheduler <- function(gamma0 = 1, a = 1, K = 1, gamma1, n1) {
  force(a)
  if (!missing(gamma1) && !missing(n1))
    K <- n1^a * gamma1 / (gamma0 - gamma1)
  b <- gamma0 * K
  function(n) b / (K + n^a)
}

horses <- read_csv("P:/KU/CompStat/Data/Horses.csv")
horses <- na.omit(horses)


get_knots <- function(inner_knots) {
  sort(c(rep(range(inner_knots), 3), inner_knots))
}

Phi <- function(x_vec, knots)
  splineDesign(knots, x_vec)

grad_two_norm_squared <- function(beta_vec, inner_knots)
  pen_mat(inner_knots) %*% beta_vec

p <- function(beta, x_vec, knots) {
  inverse_logit(Phi(x_vec, knots) %*% beta)
}

p_fast <- function(beta, phi_mat) {
  inverse_logit(phi_mat %*% beta)
}

grad_L <- function(beta, lambda, x_vec, y_vec, knots, inner_knots) {
  -crossprod(Phi(x_vec, knots),
             y_vec - p(beta, x_vec, knots)) / length(x_vec) +
    lambda * 2 * grad_two_norm_squared(beta, inner_knots)
}

get_objective_H <- function(lambda, x_vec, y_vec, knots, inner_knots) {
  force(x_vec); force(y_vec); force(knots);force(inner_knots); force(lambda)
  N <- length(x_vec)
  function(beta) {
    -1/N * (crossprod(y_vec, log(p(beta, x_vec, knots))) + 
            crossprod(1 - y_vec, log(1 - p(beta, x_vec, knots)))) + 
    lambda * crossprod(beta, grad_two_norm_squared(beta, inner_knots))
  }
}

get_grad_L_fast <- function(lambda, x_vec, y_vec, knots, inner_knots) {
  force(x_vec); force(y_vec); force(knots);force(inner_knots); force(lambda)
  Phi_mat <- Phi(x_vec, knots)
  function(beta, i) {
    Phi_i <- matrix(Phi_mat[i,], nrow = length(i))
    -crossprod(Phi_i,
             y_vec[i] - p_fast(beta, Phi_mat[i,])) / length(x_vec[i]) +
    lambda * 2*grad_two_norm_squared(beta, inner_knots)
  }
}


get_grad_L <- function(lambda, x_vec, y_vec, knots, inner_knots) {
  force(x_vec); force(y_vec); force(knots);force(inner_knots); force(lambda)
  function(beta, i) {
    -crossprod(Phi(x_vec[i], knots),
             y_vec[i] - p(beta, x_vec[i], knots)) / length(x_vec[i]) +
    lambda * grad_two_norm_squared(beta, inner_knots)
  }
}


#####################################################
### STOCHASTIC GRADIENT DESCENT
#####################################################

sgd <- function(n_iter, beta_init, lambda, knots, inner_knots,x, y,
                decay_schedule, batch_size = 1, cb = NULL) {
  beta_new <- beta_init
  n <- length(x)
  for (i in 1:n_iter) { 
    beta_old <- beta_new
    indexes <- sample(n, batch_size, replace = T)
    x_vec <- x[indexes]
    y_vec <- y[indexes]
    beta_new <- beta_old - 
      decay_schedule(i) * grad_L(beta_old, lambda, x_vec, y_vec, knots, inner_knots)
    beta_v <- as.vector(beta_new)
    if (!is.null(cb)) cb()
  }
  beta_new
}


#####################################################
### STOCHASTIC GRADIENT DESCENT - NRH VERSION
#####################################################


SG <- function(
  par, 
  grad,              # Function of parameter and observation index
  N,                 # Sample size
  gamma,             # Decay schedule or a fixed learning rate
  maxiter = 100,     # Max epoch iterations
  sampler = sample,  # How data is resampled. Default is a random permutation
  cb = NULL, 
  ...
) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter) 
  for(k in 1:maxiter) {
    samp <- sampler(N)   
    if(!is.null(cb)) cb()
    for(j in 1:N) {
      i <- samp[j]
      par <- par - gamma[k] * grad(par, i, ...)
    }
  }
  par
}



################################
### BATCH
###############################

batch <- function(
    par, 
    samp,
    gamma,  
    grad,              # Function of parameter and observation index
    m = 50,            # Mini-batch size
    ...
  ) {
    M <- floor(length(samp) / m)
    for(j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      par <- par - gamma * grad(par, i, ...)
    }
    par
}

################################
### SG-full
###############################

SGFULL <- function(
  par, 
  N,                 # Sample size
  gamma,             # Decay schedule or a fixed learning rate
  epoch = batch,     # Epoch update function
  ...,               # Other arguments passed to epoch updates
  maxiter = 100,     # Max epoch iterations
  sampler = sample,  # How data is resampled. Default is a random permutation
  cb = NULL
) {
  gamma1 <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter) 
  for(k in 1:maxiter) {
    if(!is.null(cb)) cb()
    samp <- sampler(N)
    par <- epoch(par, samp, gamma1[k], ...)
  }
  par
}


################################
### MOMENTUM
###############################

momentum <- function() {
  rho <- 0
  function(
    par, 
    samp,
    gamma,  
    grad,
    m = 50,             # Mini-batch size
    beta = 0.95,        # Momentum memory 
    ...
  ) {
    M <- floor(length(samp) / m)
    for(j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      # Using '<<-' assigns the value to rho in the enclosing environment
      rho <<- beta * rho + (1 - beta) * grad(par, i, ...)  
      par <- par - gamma * rho
    }
    par
  }
}

#####################################################
### NEWTON RAPHSON
#####################################################

newton_update <- function(num_param, x, y, lambda) {
  force(num_param); force(x); force(y); force(lambda)
  NUM         <- num_param
  n_sim       <- length(x)
  inner_knots <- seq(min(x), max(x), length.out = NUM - 2)
  knots       <- get_knots(inner_knots)
  Phix        <- Phi(x, knots)
  Omegax      <- pen_mat(inner_knots)
  
  W <- function(beta) {
    p1 <- as.vector(p(beta, x, knots))
    diag(p1*(1 - p1))
  }
  
  Winv <- function(beta) {
    p1 <- as.vector(p(beta, x, knots))
    diag(1/(p1*(1 - p1)))
  }
  
  z <- \(beta) Phix %*% beta + Winv(beta) %*% (y - p(beta, x, knots))
  updateStep <-
    \(beta) solve(t(Phix) %*% W(beta) %*% Phix + n_sim * lambda * Omegax) %*% t(Phix) %*% W(beta) %*% z(beta)
  
  updateStep
}


NEWTON <- 
  function(x,y, max_iter = 100, init_value,epsilon = 1e-5,num_param, lambda = 0, cb = NULL) {
  beta_old <- init_value
  newton_step <- newton_update(num_param, x, y, lambda)
  for (i in 1:max_iter) {
    beta_new <- newton_step(beta_old)
    if(!is.null(cb)) cb()
    if(sum((beta_new - beta_old)^2) <= epsilon * (sum(beta_new^2) + epsilon)) 
      break
    beta_old = beta_new
  }
  beta_new
}






