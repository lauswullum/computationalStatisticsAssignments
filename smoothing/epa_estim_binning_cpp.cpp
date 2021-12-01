#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;


//   http://gallery.rcpp.org/

// [[Rcpp::export]]
NumericVector rotator(NumericVector w, int num) {
  //Rotate wieghts to the left (need right)
  int m = w.size();
  for(int j = 0; j < num; ++j) {
    int temp = w.at(m - 1);
    for (int i = (m-1); i > 0 ; i--)
    {
      w.at(i) = w.at(i - 1);
    }
    w.at(0) = temp;
  }
  return w;
}


// [[Rcpp::export]]
NumericVector get_weights(NumericVector x, double delt, double minx, int n, int m) {
  NumericVector weights(m);
  for(int i = 0; i < n; ++i) {
    int s = std::floor( ((x[i] - minx) / delt) + 0.5);
    weights[s] = weights[s] + 1;
  }
  weights = weights / sum(weights);
  return weights;
}

 
// [[Rcpp::export]]
List epa_estim_binning_cpp(NumericVector x, double h, int m = 512) {
  int n = x.size();
  NumericVector xx(m);
  NumericVector y(m);
  NumericVector lwr(n);
  NumericVector upr(n);
  NumericVector bins(m);
  NumericVector weights(m);
  double minx = min(x); 
  double lwr_rg = minx - 3*h; 
  double maxx = max(x);
  double upr_rg = maxx + 3*h; 
  double lent = upr_rg - lwr_rg;
  double delt = lent / (m - 1);
  
  // compute bins and weights
  weights = get_weights(x, delt, minx = lwr_rg, n, m);
  
  // Make xx (values to evaluate kernel in) 
  // xx = seq(min(x) - 3*h, max(x) + 3*h, length.out = m)
  for(int i = 0; i < m; ++i) {
    xx[i] = lwr_rg + i*delt;
  }
  // Estimate y (kernel density estimate)
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < m; ++j) {
      lwr[j] = ((xx[i] - xx[j]) > -h ? 1 : 0);
      upr[j] = ((xx[i] - xx[j]) < h ? 1 : 0);
    }
    y[i] = sum(0.75 * (1 - (xx[i] - xx)*(xx[i] - xx) / (h*h)) * lwr * upr * weights) / h;
    
  }
  List L = List::create(Named("y") = y , _["x"] = xx);
  return L;
}
