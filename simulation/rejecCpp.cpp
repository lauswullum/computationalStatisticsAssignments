#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector post_rejec_cpp(int n, NumericVector x_vec, NumericVector z_vec, double alpha_prime, bool trace = false) {
  NumericVector y(n);
  double y0;
  int data_length = x_vec.size();
  bool reject;
  for(int i = 0; i < n; ++i) {
    do {
      y0 = R::rnorm(0, 1);
      double f_dens = 1; 
      bool geq = (y0 >= 0);
      for(int i = 0; i < data_length; i++) {
        f_dens = f_dens *  exp(y0 * z_vec[i] * x_vec[i] - exp(y0 * x_vec[i]));
      }
      double gauss_dens = exp(-0.5*y0);
      reject = R::runif(0, 1) > (geq * f_dens * alpha_prime / gauss_dens ) ;
    } while(reject);
    y[i] = y0;
  }
  return y;
}
