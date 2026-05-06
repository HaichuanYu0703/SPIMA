// Simulate continuous pseudo-individual data (Rcpp accelerated)
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix simulate_cont_studies(NumericVector n_ctrl,
                                     NumericVector n_trt,
                                     NumericVector mean_ctrl,
                                     NumericVector sd_ctrl,
                                     NumericVector sd_trt,
                                     NumericVector theta_i) {
  // Generates normal pseudo-IPD for both arms.
  // For each study:
  //   Control:  y ~ N(mean_ctrl[i], sd_ctrl[i])
  //   Treatment: y ~ N(mean_ctrl[i] + theta_i[i], sd_trt[i])
  //
  // Returns: matrix with 3 columns [study_id, group, y]
  //   group: 0 = control, 1 = treatment
  int n_studies = n_ctrl.size();
  int total_n = sum(n_ctrl) + sum(n_trt);

  NumericMatrix out(total_n, 3);
  int row = 0;

  for (int i = 0; i < n_studies; i++) {
    double m_c = mean_ctrl[i];
    double s_c = sd_ctrl[i];
    double s_t = sd_trt[i];
    double mu_t = m_c + theta_i[i];

    // Control arm (group = 0)
    for (int j = 0; j < n_ctrl[i]; j++) {
      out(row, 0) = i + 1;
      out(row, 1) = 0;          // control
      out(row, 2) = R::rnorm(m_c, s_c);
      row++;
    }

    // Treatment arm (group = 1)
    for (int j = 0; j < n_trt[i]; j++) {
      out(row, 0) = i + 1;
      out(row, 1) = 1;          // treatment
      out(row, 2) = R::rnorm(mu_t, s_t);
      row++;
    }
  }
  return out;
}
