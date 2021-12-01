
data {
  int<lower=0> NUM;               // Number of Observations
  int<lower=0> z[NUM];            // Observed response
  vector[NUM] x;                  // Covariate
}

parameters {
  real y;                     // Theta parameter
}
model {
  //y ~ uniform(0, 10);
  z ~ poisson_log(y * x);
}
