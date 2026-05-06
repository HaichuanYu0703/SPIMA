// Simulate binary pseudo-individual data (Rcpp accelerated)
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix simulate_binary_studies(NumericVector n_ctrl,
                                       NumericVector n_trt,
                                       NumericVector p_ctrl,
                                       NumericVector theta_i) {
  // Generates binary (Bernoulli) pseudo-IPD for both arms.
  // For each study:
  //   Control:  y ~ Bernoulli(p_ctrl[i])
  //   Treatment: logit(p) = logit(p_ctrl[i]) + theta_i[i]
  //
  // Returns: matrix with 3 columns [study_id, group, y]
  //   group: 0 = control, 1 = treatment
  int n_studies = n_ctrl.size();
  int total_n = sum(n_ctrl) + sum(n_trt);

  NumericMatrix out(total_n, 3);
  int row = 0;

  for (int i = 0; i < n_studies; i++) {
    // Control arm (group = 0)
    int nc = n_ctrl[i];
    for (int j = 0; j < nc; j++) {
      out(row, 0) = i + 1;       // study_id (1-based)
      out(row, 1) = 0;           // control
      out(row, 2) = R::rbinom(1, p_ctrl[i]);
      row++;
    }

    // Treatment arm (group = 1)
    double logit_pc = std::log(p_ctrl[i] / (1.0 - p_ctrl[i] + 1e-10));
    double p_trt = 1.0 / (1.0 + std::exp(-(logit_pc + theta_i[i])));
    int nt = n_trt[i];
    for (int j = 0; j < nt; j++) {
      out(row, 0) = i + 1;
      out(row, 1) = 1;           // treatment
      out(row, 2) = R::rbinom(1, p_trt);
      row++;
    }
  }
  return out;
}
