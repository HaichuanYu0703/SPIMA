// Fast distance computations (Rcpp accelerated)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double dist_euclidean(NumericVector x, NumericVector y) {
  int n = x.size();
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    double d = x[i] - y[i];
    sum += d * d;
  }
  return std::sqrt(sum);
}

// [[Rcpp::export]]
double dist_weighted_euclidean(NumericVector x, NumericVector y,
                                NumericVector weights) {
  int n = x.size();
  double sum = 0.0;
  double w_sum = 0.0;
  for (int i = 0; i < n; i++) {
    double d  = x[i] - y[i];
    double w  = weights[i];
    sum += w * d * d;
    w_sum += w;
  }
  return std::sqrt(sum / (w_sum + 1e-10));
}
