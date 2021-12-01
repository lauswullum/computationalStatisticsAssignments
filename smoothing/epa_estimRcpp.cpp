#include <Rcpp.h>
using namespace Rcpp;

//   http://gallery.rcpp.org/

  
// [[Rcpp::export]]
List epa_estim_cpp(NumericVector x, double h, int m = 512) {
  int n = x.size();
  NumericVector xx(m);
  NumericVector y(m);
  NumericVector lwr(n);
  NumericVector upr(n);
  double minx = min(x); 
  double maxx = max(x);
  double lent = maxx + 3*h - minx + 3*h;
  double delt = lent / (m - 1);
  // Make xx (values to evaluate kernel in) 
  // xx = seq(min(x) - 3*h, max(x) + 3*h, length.out = m)
  for(int i = 0; i < m; ++i) {
    xx[i] = minx - 3*h + i*delt;
  }
  // Estimate y (kernel density estimate)
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      lwr[j] = ((xx[i] - x[j]) > -h ? 1 : 0);
      upr[j] = ((xx[i] - x[j]) < h ? 1 : 0);
    }
    y[i] = sum(0.75 * (1 - (xx[i] - x)*(xx[i] - x) / (h*h)) * lwr * upr) / (h * n);
    
  }
  List L = List::create(Named("y") = y , _["x"] = xx);
  return L;
}


