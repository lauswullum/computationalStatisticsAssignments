
#####################################################
### PACKAGES 
#####################################################

library(CSwR)
library(caret)
library(tidyverse)
library(reshape2)
library(tibble)
library(LaplacesDemon)
library(see)
library(ggplot2)
library(latex2exp)
library(ggthemr)
library(patchwork)
library(bench)

#####################################################
### FACTORY
#####################################################



EM_factory <- function(n, nu, mu, sigma2) {
  force(n);force(nu);force(mu);force(sigma2)
  ######################
  # SIMULATE DATA 
  ######################
  w <- rchisq(n, df = nu)
  x <- rnorm(n, mean = mu, sd = sqrt(nu*sigma2 / w))
  
  ######################
  # Q(THETA | THETA')
  ######################
  Q <- function(par, par_prime) {
      - n * log(par[2]) / 2 -
    sum((x - par[1])^2 / (1 + (x - par_prime[1])^2 / (nu * par_prime[2]))) *
    (1 + nu) / (2 * nu * par[2])
  }
  
  ######################
  # W-step
  ######################
  W_step <- function(par) {
    # Generate beta_list and alpha for the conditional mean
    beta_list   <- 0.5 * (1 + ((x - par[1])^2 /(par[2] * nu))  )
    alpha       <- (nu + 1) * 0.5
    alpha/beta_list
  }

  ######################
  # EM STEP
  ######################
  EM_step <- function(par) {

    w_mean <- W_step(par)

    # update mu_hat and sigmasq_hat
    mu_hat      <- sum(x*w_mean) / sum(w_mean)
    sigmasq_hat <- (1 / (nu * n)) * sum(w_mean * (x - par[1])^2)
    par         <- c(mu_hat, sigmasq_hat)
  }
  
  ######################
  # Q GRADIENT
  ######################

  Q_grad <- function(par, par_prime) {
    w_mean <- W_step(par_prime)

    # compute gradient
    mu_grad <- sum(w_mean*(x - par[1])) / (nu * par[2]) 
    sigma2_grad <- - n/(2*par[2]) + sum(w_mean * (x - par[1])^2) / (2 * nu * par[2]^2) 
    c(mu_grad, sigma2_grad)
  }
  
  ######################
  # FISHER INFORMATION
  ######################
  EM_fisher <- function(par, center = T) {
      
    w_mean <- W_step(par)
    
    # Gradient lists
    mu_grads <- w_mean * ((x - par[1]) / (nu * par[2]) )
    sigmasq_grads <- -1/(2*par[2]) + (w_mean * (x - par[1])^2) / (2 * nu * par[2]^2 ) 
    grad_list <- cbind(mu = mu_grads, sigma2 = sigmasq_grads)
    
    if (center) { 
      #t(grad_list) %*% grad_list
      crossprod(grad_list, grad_list)
    }
    else {
      summen <- 0
      mean_mu <- mean(mu_grads)
      mean_sigmasq <- mean(sigmasq_grads)
      grad_loglik_sum <- c(mean_mu, mean_sigmasq)
      print(grad_loglik_sum)
      for (i in 1:n) {
        summen <- summen + (grad_list[i,] - grad_loglik_sum) %*% t(grad_list[i,] - grad_loglik_sum)
      }
      summen
    }
  }
   
  ######################
  # COMPLETE OPT
  ######################
  complete_loglik_opt <- function() {
    mu <- sum(w*x) / sum(w)
    K <- sum(w * (x - mu)^2)
    sigma_sq <- 1/(n*nu)*K
    c(mu = mu, sigma_sq = sigma_sq)
  }

  #########################
  # MARGINAL LOGLIKELIHOOD 
  #########################

  marg_loglik <- function(par) {
    -n*log(gamma(0.5*(nu + 1))) + n*0.5*log(nu *par[2]*pi) + n*log(gamma(0.5*nu)) + (nu + 1)*0.5*sum( log(1 + (x - par[1])^2 / (nu * par[2])))
  }

  #########################
  # MARGINAL GRADIENTs
  #########################
  marg_grad_i <- function(par, i) {
    mu_grad <- 
      -sum(((nu + 1)*(x[i] - par[1])) / ( (x[i] - par[1])^2 + nu * par[2]))
    sigma_grad <- 
      n/(2*par[2]) - sum( (nu + 1)/(2*par[2]) * ((x[i] - par[1])^2 / (nu*par[2] + (x[i] - par[1])^2)) )
    
    c(mu_grad, sigma_grad)
  }
  
  marg_grad <- function(par) {
    
    mu_grad <- 
      -sum(((nu + 1)*(x - par[1])) / ( (x - par[1])^2 + nu * par[2]))
    sigma_grad <- 
      n/(2*par[2]) - sum( (nu + 1)/(2*par[2]) * ((x - par[1])^2 / (nu*par[2] + (x - par[1])^2)) )
    
    c(mu_grad, sigma_grad)
  }
  
  
  list(
    w = w,
    x = x,
    Q = Q,
    EM_step = EM_step,
    EM_fisher = EM_fisher,
    Q_grad = Q_grad, 
    complete_loglik_opt = complete_loglik_opt, 
    marg_grad = marg_grad, 
    marg_grad_i = marg_grad_i, 
    marg_loglik = marg_loglik
  )
}



#####################################################
### EM - ALGORTIHM
#####################################################

EM <- function(par, EM_step, max_iter = 1000, epsilon = 1e-15, trace = NULL, ...) {
  
  for (i in 1:max_iter) {
    par0  <- par
    # Update mu_hat and sigmasq_hat
    par   <- EM_step(par0, ...)
    
    # Tracer
    if(!is.null(trace)) trace()
    
    # Convergence Criterion
    if(sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon)) 
      break
  }
  par
}


#####################################################
### OPTIMIZATION OF THE MARGINAL LIKELIHOOD
#####################################################

GD <- function(
  par, 
  H,
  gr,
  d = 0.8, 
  c = 0.1, 
  gamma0 = 0.01, 
  epsilon = 1e-4, 
  maxiter = 1000,
  cb = NULL
) {
  for(i in 1:maxiter) {
    value <- H(par)
    grad <- gr(par)
    h_prime <- sum(grad^2)
    # Convergence criterion based on gradient norm
    #if(h_prime <= epsilon) break
    gamma <- gamma0
    # Proposed descent step
    par1 <- par - gamma * grad
    # Backtracking while descent is insufficient
    while(H(par1) > value - c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par - gamma * grad
    }
    if(sum((par - par1)^2) <= epsilon * (sum(par^2) + epsilon)) 
      break
    if(!is.null(cb)) cb()
    par <- par1
  }
  if(i == maxiter)
    warning("Maximal number, ", maxiter, ", of iterations reached")
  par
}

SG <- function(
  par, 
  N,                 # Sample size
  gamma,             # Decay schedule or a fixed learning rate
  epoch = batch,     # Epoch update function
  ...,               # Other arguments passed to epoch updates
  maxiter = 100,     # Max epoch iterations
  sampler = sample,  # How data is resampled. Default is a random permutation
  cb = NULL
) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter) 
  for(k in 1:maxiter) {
    if(!is.null(cb)) cb()
    samp <- sampler(N)
    par <- epoch(par, samp, gamma[k], ...)
  }
  par
}

######################################
### BATCH
######################################

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

decay_scheduler <- function(gamma0 = 1, a = 1, K = 1, gamma1, n1) {
  force(a)
  if (!missing(gamma1) && !missing(n1))
    K <- n1^a * gamma1 / (gamma0 - gamma1)
  b <- gamma0 * K
  function(n) b / (K + n^a)
}


  
