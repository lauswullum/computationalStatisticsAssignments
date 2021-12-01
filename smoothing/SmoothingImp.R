#####################################################
### PACKAGES 
#####################################################
library(Rcpp)
library(caret)
library(tidyverse)
library(reshape2)
library(tibble)
library(LaplacesDemon)
library(see)
library(knitr)
library(ggplot2)
library(latex2exp)
library(ggthemr)
library(patchwork)
library(bench)

#####################################################
### UTILITY FUNCTIONS
#####################################################
# Likelihood loss 
loglik_loss <- \(x) sum(log(x))

# Unbiased cross validation loss
UCV_loss <- function(x, n) 3/(5*n*h) - 2/n * sum(x)

# make CV sets
make_cv_fast <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# The epanenchnikoc kernel
epa_ker <- \(x) 3 / 4 * (1 - x ^ 2) * (x > -1) * (x < 1)

# The second derivative of the gaussian kernel
gauss_sec_deriv <- \(x) - sqrt(2 / pi) * (exp(-x ^ 2) - 2 * x ^ 2 * exp(-x^2))

# Robust estimate of variance by Silverman
sigma_tilde <- \(x) min(var(x), IQR(x) / 1.34)

# Silvermans rule of thumb
silvermans_rt <- function(x) {
  n <- length(x)
  0.9 * sigma_tilde(x) * n^(-0.2)
}

silvermans_rt_epa <- function(x) {
  n <- length(x)
  (8*sqrt(pi)*5)^0.2 * sigma_tilde(x) * n^(-0.2)
}

# Fourth derivative of the gaussian kernel with normalization
#   Loader, Clive R. “Bandwidth Selection: Classical or Plug-In?” 
#   The Annals of Statistics 27, no. 2 (1999): 415–38.
gauss_fourth_deriv <- \(x) sqrt(1/(2*pi))*exp(-x^2 / 2)*(x^4 - 6*x^2 + 3)

# Estimate of where to search for h in optimization
find_range_of_h <- function(x) {
  max_dist <- max(x) - min(x)
  n <- length(x) 
  len_mat <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      len_mat[i, j] <- abs(x[i] - x[j])
    }
  }
  len_vec <- as.vector(len_mat) |> na.omit()
  ret <- len_vec[len_vec != 0]
  quantile(ret, c(0.25, 0.75))
}

#####################################################
### KERNEL ESTIMATION PROCEDURES
#####################################################

# Epanenchnikov kernel estimation
epa_estim <- function(x, h, m = 512) {
  xx <- seq(min(x) - 3*h, max(x) + 3*h, length.out = m) # grid points
  y <- numeric(m)
  for (i in seq_along(xx)) {
    y[i] <- sum(3/4 * (1 - (xx[i] - x)^2 / h^2) * ((xx[i] - x) >= -h) * ((xx[i] - x) <= h))
  } 
  list(x = xx, y = y / (length(x)*h))
}

epa_estim_test <- function(x, h, m = 512, xx = seq(min(x) - 3*h, max(x) + 3*h, length.out = m)) {
  y <- numeric(m)
  for (i in seq_along(xx)) {
    y[i] <- sum(3/4 * (1 - (xx[i] - x)^2 / h^2) * ((xx[i] - x) >= -h) * ((xx[i] - x) <= h))
  } 
  list(x = xx, y = y / (length(x)*h))
}

# Kernel estimation at a single observation.
epa_estim_single <- function(x, x_obs, h) {
  y <- sum(3/4 * (1 - (x_obs - x)^2 / h^2) * ((x_obs - x) > -h) *((x_obs - x) < h) )
  y /  (h * length(x))
}


cppFunction('double epa_estim_single_cpp(NumericVector x, double x_obs ,double h, int m = 512) {
  int n = x.size();
  NumericVector lwr(n);
  NumericVector upr(n);
  NumericVector proddet(n); 
  for(int i = 0; i < n; ++i) {
    lwr[i] = ((x_obs - x[i]) > -h ? 1 : 0);
    upr[i] = ((x_obs - x[i]) < h ? 1 : 0);
  }
  return sum(0.75 * (1 - (x_obs - x)*(x_obs - x) / (h*h)) * lwr * upr) / (h * n);
}')

kern_bin <- function(x, lo, hi, m) {
  w <- numeric(m)
  delta <- (hi - lo) / (m - 1)
  for(j in seq_along(x)) {
    i <- floor((x[j] - lo) / delta + 0.5) + 1
    w[i] <- w[i] + 1
  }
  w / sum(w)
}

epa_estim_binning <- function(x, h, m = 512) {
  n <- length(x)
  rg <- range(x) + c(- 3 * h, 3 * h) 
  xx <- seq(rg[1], rg[2], length.out = m)
  weights <- kern_bin(x, rg[1], rg[2], m)
  delt = xx[2] - xx[1]
  w <- get_weights(x, delt, minx = rg[1], n, m)
  kerneval <- 3/4 * (1 - (xx - xx[1])^2 / h^2) * ((xx - xx[1]) > -h) *((xx - xx[1]) < h) / h
  kerndif <- toeplitz(kerneval)
  y <- colSums(w * kerndif)
  list(x = xx, y = y, h = h)
}
# Binning in c++
#sourceCpp(file = "epa_estimRcpp.cpp")
#sourceCpp(file = "epa_estim_binning_cpp.cpp")
sourceCpp(file = "P:/KU/CompStat/Ass/smoothing/epa_estim_binning_cpp.cpp")
sourceCpp(file = "P:/KU/CompStat/Ass/smoothing/epa_estimRcpp.cpp")

#####################################################
### CROSS VALIDATION ###
#####################################################

#LOOCV
compute_loglik <- function(x, h) {
  fihat_vec <- numeric(length(x))
  for (i in seq_along(x)) {
    loo_set <- x[-i]
    x_eval <- x[i]
    fihat_vec[i] <- epa_estim_single(loo_set, x_eval, h)
  }
  # loss function
  sum(log(fihat_vec))
}

# GENERAL CV with k number of folds
CV_fast <- function(x, h, k, loss = loglik_loss, ...) {
  n <- length(x)
  # Create fold indices via the caret package
  fold_idx <- createFolds(x, k = k) %>% unname
  fihat_vec <- numeric(n)
  
  for(i in 1:k) {
    # Leave one fold out
    #other_fold_idx <- .Internal(unlist(fold_idx[-i], FALSE, FALSE))
    other_fold_idx <- unlist(fold_idx[-i])
    # Fold to evaluate
    eval_fold_idx <- fold_idx[[i]] 
    # Evaluations
    eval_set <- x[eval_fold_idx]
    eval_other <- x[other_fold_idx]
    # Compute f_i hat in all left out
    for(j in eval_fold_idx) {
      fihat_vec[j] <- epa_estim_single_cpp(eval_other, x[j], h)
    }
  }
  loss(fihat_vec)
}

# GENERAL CV WITH K number of folds
CV_slow <- function(x, h, k, loss = loglik_loss, ...) {
  n <- length(x)
  # Create fold indices via the caret package
  folds <- make_cv_fast(x, n = k)
  fihat_vec <- numeric(n)
  count <- 1 
  for(i in 1:k) {
    # other folds
    other_fold <- unlist(folds[-i])
    #other_fold <- .Internal(unlist(folds[-i], FALSE, FALSE))
    # Leave one fold out
    eval_fold <- folds[[i]]
    # Compute f_i hat in all left out
    for (j in seq_along(eval_fold)) {
      fihat_vec[count] <- epa_estim_single(other_fold, eval_fold[j], h)
      count <- count + 1
    }
  }
  loss(fihat_vec)
}

#####################################################
### PLUG IN ESTIMATION ###
#####################################################

# ||f''_h||_2^2 - curvature estimate using pilot curve
pilot_curve_estimate <- function(x,h) {
  n <- length(x)
  normalizer <- n^2*(sqrt(2)*h)^5
  ret_sum <- 0   
  for (xi in x) {
    ret_sum <- ret_sum + sum(gauss_fourth_deriv((xi - x)/(sqrt(2)*h)))
  }
  ret_sum / normalizer
}

# Plug in estimation
plug_in_estimate <- function(x) {
  # pilot bandwidth
  #pilot_bandwidth <- silvermans_rt(x)
  pilot_bandwidth <- silvermans_rt_epa(x)
  # pilot curvature estimate
  fbarbar <- pilot_curve_estimate(x, pilot_bandwidth)
  # Oracle bandwidth
  (15 / (fbarbar * length(x)))^(1 / 5) 
}


#####################################################
### UNBIASED CROSS-VALIDATION CRITERION ###
#####################################################

# generate LOO sets and compute cross validation
compute_UCV <- function(x, h) {
  n <- length(x) 
  fihat_vec <- numeric(n)
  for (i in seq_along(x)) {
    loo_set <- x[-i]
    x_eval <- x[i]
    fihat_vec[i] <- epa_estim_single(loo_set, x_eval, h)
  }
  # loss function
  3/(5*n*h) - 2/n * sum(fihat_vec)
}
# Cross validation using other loss: UCV

#####################################################
### epa_estim
#####################################################

estim_gen <- function(x, len, m, h) {
  ifelse(len < 2*m, epa_estim_cpp(x, h), epa_estim_binning_cpp(x, h))
}
own_density <- function(x, m = 512, method = "default", interval, ...) {
  len <- length(x)
  method_list = c("default", "loocv", "kcv", "plugin")
  method_list_p = c(" default ", " loocv ", " kcv ", " plugin ")
  if (is.na(match(method, method_list))) {
    stop(method,
         " must be one of ",
         method_list_p,
         " or none if you want to use default method. ")
  }
  if (method == "default") {
    h <- silvermans_rt_epa(x)
    estim_gen(x, len, m, h)
  }
  else if (method == "loocv") { 
    opt_loocv <- optimize(\(h) compute_loglik(sett, h), interval, maximum = T)
    estim_gen(x, len, m, h = opt_loocv$maximum)
  }
  else if (method == "kcv") { 
    opt_kcv <- optimize(\(h) CV_fast(x, h,...), interval, maximum = T)
    estim_gen(x,len, m, opt_kcv$maximum)
  }
  else { 
    h <- plug_in_estimate(x)
    estim_gen(x,len, m, h)
  }
}
