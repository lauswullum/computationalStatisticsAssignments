
#####################################################
### PACKAGES 
#####################################################
library(Rcpp)
library(posterior)
library(cmdstanr)
library(bayesplot)
library(tidyverse)
library(reshape2)
library(tibble)
library(knitr)
library(ggplot2)
library(latex2exp)
library(patchwork)
library(bench)

#####################################################
### UTILITY FUNCTIONS
#####################################################

# Read in data
pois_data <- read_csv("P:/KU/CompStat/Ass/simulation/Poisson.csv")
x_vec <- pois_data$x
z_vec <- pois_data$z

# Implement density function
get_f_dens_vec <- function(predictor, response) {
  f_dens        <- \(y) prod( exp(y * response * predictor - exp(y*predictor) ) )
  f_dens_vec    <- \(y) map_dbl(y, f_dens)
  f_dens_vec
}

f_dens        <- \(y) prod( exp(y * z_vec * x_vec - exp(y*x_vec) ) )
f_dens_vec    <- \(y) map_dbl(y, f_dens)

# Log-space target density
flog_dens     <- \(y)  (y >= 0) * sum(y*x_vec*z_vec - exp(y*x_vec))
flog_dens_vec <- \(y)  map_dbl(y, flog_dens)


gauss_log_proposal <- \(y) - 0.5 * y^2
gauss_proposal <- \(y) exp(-0.5 * y^2)
gauss_proposal_mean <- function(mean) \(y) exp(-0.5 * (y - mean)^2)

# Get the alpha prime 
alpha_prime_fun_log <- \(y) gauss_log_proposal(y) - flog_dens(y)
alpha_prime_log <- optimize(alpha_prime_fun_log, interval = c(0,3))$objective
alpha_prime <- exp(alpha_prime_log)

#####################################################
### ALPHA PRIMER
#####################################################

get_starting_values <- function(predictor, response) {
  mod     <- glm(response ~ predictor - 1, family = "poisson")
  conf    <- suppressMessages(confint(mod, level = 0.99))
  confLwr <- conf[1] 
  confUpr <- conf[2] 
  list(coef = mod$coefficients, confLwr, confUpr)
}

simulate_data <- function(num, beta,...) {
  predictor <- rnorm(num, ...)
  response <- rpois(num, exp(beta * predictor) )
  list(predictor = predictor, response = response )
}


alpha_prime_finder <- function(
 log_target, 
 predictor, 
 response,
 interval = c(0,40)
 ) {
  mode = coef(glm(response ~ predictor - 1, family = "poisson"))
  gauss_log_proposal <- \(y) - 0.5 * (y - mode)^2
  alpha_prime_fun_log <- \(x) gauss_log_proposal(x) - log_target(x)
  optimize(alpha_prime_fun_log, interval = interval)
}

alpha_prime_ny <- alpha_prime_finder(flog_dens, x_vec, z_vec)


#####################################################
### REJECTION SAMPLING 
#####################################################

post_rejec <- function(n, trace = F) {
  y <- numeric(n)
  count <- 0 
  for(i in 1:n) {
    reject <- TRUE
    while (reject){
      count  <- count + 1
      y0     <- rnorm(1)
      u      <- runif(1)
      reject <- u > ((f_dens(y0) * alpha_prime) / gauss_proposal(y0)) 
    }
    y[i] <- y0
  }
  if(trace)
    cat("Rejection freq:", (count - n)/ count, "\n")  ## Rejection frequency
  y
}


#####################################################
### C++ IMPLEMENTATION 
#####################################################
sourceCpp(file = "P:/KU/CompStat/Ass/simulation/rejecCpp.cpp" )


#####################################################
### VECTORIZED IMPLEMENTATION
#####################################################

post_rejec_vec <- function(n) {
  fact <- 1
  j    <- 1
  l    <- 0  ## The number of accepted samples
  y    <- list()
  while(l < n) {
    m      <- floor(fact * (n - l))  ## equals n the first time
    y0     <- rnorm(m)
    u      <- runif(m)
    accept <- (u <= (f_dens_vec(y0) * alpha_prime) / gauss_proposal(y0))
    l      <- l + sum(accept)
    y[[j]] <- y0[accept]
    j      <- j + 1
    if(fact == 1) fact <- n / l
  }
  unlist(y)[1:n]
}


#####################################################
### ADAPTIVE ENVELOPING
#####################################################

post_adap <- function(n, x1, x2, response, predictor, trace = F ) {
  num <- length(response)
  prod1 <- response * predictor
  lf        <- \(x) map_dbl(x, \(x) sum( x * prod1 - exp(x * predictor) ))
  lf_deriv  <- \(x) map_dbl(x, \(x) sum( prod1 - predictor * exp(x * predictor) ))
  a1 <- lf_deriv(x1)
  a2 <- lf_deriv(x2)
  b1 <- lf(x1) - a1 * x1
  b2 <- lf(x2) - a2 * x2
  z0 <- 0
  z2 <- 1
  z1 <- (b2 - b1) / (a1 - a2)
  Q1 <- exp(b1) * (exp(a1 * z1) - exp(a1*z0)) / a1 
  c  <- Q1 + exp(b2) * (exp(a2 * z2) - exp(a2 * z1)) / a2
  y  <- numeric(n)
  count <- 0
  for(i in 1:n) {
    reject <- TRUE
    while (reject) {
      count <- count + 1
      uy    <- runif(1)
      u     <- runif(1)
      u0    <- c * uy
      if (u0 < Q1) {
        z <- log(a1 * exp(-b1) * u0 + 1) / a1
        reject <- u > exp(lf(z) - a1 * z - b1)
      } else {
        z <- log(a2 * exp(-b2) * (u0 - Q1) + exp(a2 * z1)) / a2
        reject <- u > exp(lf(z) - a2 * z - b2)
      }
    }
    y[i] <- z
  }
  if(trace)
    cat("Rejection freq:", (count - n)/ count, "\n")  ## Rejection frequency
  y
}


#####################################################
### TRULY ADAPTIVE REJECTION SAMPLING - OLD
#####################################################

make_grid <- function(len, response, predictor) {
  start_vals <- get_starting_values(predictor, response)
  seq(from = start_vals[[2]], to = start_vals[[3]], length.out = len )
}


#####################################################
### ADAPTIVE ENVELOPE
#####################################################

make_envelope_utils <- function(response, predictor, len, endpoint) {
  grid  <- make_grid(len, response, predictor)
  #grid <- c(0.2, 0.28)
  prod1 <- response * predictor
  
  # lf's
  lf        <- \(x) map_dbl(x, \(x) sum( x * prod1 - exp(x * predictor) ))
  lf_deriv  <- \(x) map_dbl(x, \(x) sum( prod1 - predictor * exp(x * predictor) ))
  
  # a and b list
  a_list <- lf_deriv(grid)
  b_list <- lf(grid) - a_list * grid
  
  # cut points
  z_list        <- numeric(len + 1)
  z_list[1]     <- 0 # Startpoint
  z_list[len+1] <- endpoint # EndPoint
  z_list[2:len] <- ( b_list[2:len] - b_list[1:len - 1] ) / (a_list[1:(len - 1)] - a_list[2:len])
  
  Q_list    <- numeric(len)
  Q_list[1] <- exp(b_list[1]) * ( exp(a_list[1] * z_list[2]) - exp(a_list[1] * z_list[1]) ) / a_list[1]
  for (i in 2:len) {
    Q_list[i] <- Q_list[i - 1] + exp(b_list[i]) * ( exp(a_list[i] * z_list[i + 1]) - exp(a_list[i] * z_list[i]) ) / a_list[i]
  }
  return(list(
    Q_list = Q_list,
    z_list = z_list,
    a_list = a_list,
    b_list = b_list
  ))
}


post_adap_more <- function(n, response, predictor, len, trace = F, endpoint = 1) {
  envelope_utils <- make_envelope_utils(response, predictor, len, endpoint)
  # Recomputing??
  prod1 <- response * predictor
  lf        <- \(x) map_dbl(x, \(x) sum( x * prod1 - exp(x * predictor) ))
  lf_deriv  <- \(x) map_dbl(x, \(x) sum( prod1 - predictor * exp(x * predictor) ))
  
  Q_list  <- envelope_utils$Q_list
  c       <- Q_list[len] # fine
  z_list  <- envelope_utils$z_list
  a_list  <- envelope_utils$a_list
  b_list  <- envelope_utils$b_list
  
  y <- numeric(n)
  count <- 0 
  for(i in 1:n) {
    reject <- TRUE
    while(reject) {
      count   <- count + 1
      u       <- runif(1)
      uy      <- runif(1) 
      u0      <- c * uy
      u0_idx  <- which.max(u0 <= Q_list)
      if (u0 < Q_list[1]) {
        z      <- log(a_list[1] * exp(-b_list[1]) * u0 + exp(a_list[1] * z_list[1]) ) / a_list[1]
        reject <- u > exp(lf(z) - a_list[1] * z - b_list[1])
      }
      else {
        z      <- log( a_list[u0_idx] * exp(-b_list[u0_idx]) * (u0 - Q_list[u0_idx-1]) + exp(a_list[u0_idx] * z_list[u0_idx]) ) / a_list[u0_idx]
        reject <- u > exp( lf(z) - a_list[u0_idx] * z - b_list[u0_idx])
      }
    }
    y[i] <- z
  }
  if(trace)
    cat("Rejection freq:", (count - n)/ count, "\n")  ## Rejection frequency
  y
}


#####################################################
### STAN COMPARISON
#####################################################

daten = list(NUM = 100, z = z_vec, x = x_vec)

file <- file.path("P:/KU/CompStat/Ass/simulation/stanSampling.stan")
mod_pois <- cmdstan_model(file)


fitCMDSTAN <- mod_pois$sample(
  data            = daten,
  chains          = 2,
  parallel_chains = 2,
  refresh         = 500,
  iter_sampling   = 4500
  #seed           = 124,
)
draws_array <- fitCMDSTAN$draws()
draws_df    <- as_draws_df(draws_array)
p_mcmc      <- mcmc_hist(fitCMDSTAN$draws("y"))


  
