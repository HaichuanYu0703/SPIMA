# ---- utils_skewness.R ----
# Helper functions for SE estimation from quantiles and skewness diagnosis.
# Used internally by the elastic-input quantile path in module_cont.R.

#' Approximate Standard Error of the Median from IQR
#'
#' For a normal distribution, IQR ~ 1.35 * sigma and
#' Var(median) ~ (pi/2) * sigma^2 / n.
#'
#' @param iqr Interquartile range (Q3 - Q1).
#' @param n Sample size.
#' @return Standard error of the median (scalar).
#' @keywords internal
se_from_iqr <- function(iqr, n) {
  if (n < 2 || iqr <= 0) return(NA_real_)
  sigma <- iqr / 1.35
  sqrt(1.570796 * sigma^2 / n)  # 1.570796 = pi/2
}


#' Approximate Standard Error of the Median from Range
#'
#' Uses the Wan et al. (2014) formula to estimate sigma from the range,
#' then computes SE(median) = sqrt(pi/2 * sigma^2 / n).
#'
#' Reference: Wan et al., BMC Medical Research Methodology 2014, 14:135.
#'
#' @param range_val Observed range (max - min).
#' @param n Sample size.
#' @return Standard error of the median (scalar), or NA if n < 5.
#' @keywords internal
se_from_range <- function(range_val, n) {
  if (n < 5 || range_val <= 0) return(NA_real_)
  # Clamp the qnorm argument to avoid extreme values
  p <- max(0.001, min(0.999, (n - 0.375) / (n + 0.25)))
  sigma <- range_val / (2 * stats::qnorm(p))
  sqrt(1.570796 * sigma^2 / n)
}


#' Diagnose Skewness from a Five-Number Summary
#'
#' Compares the right-side spread (Q3-median, max-median) to the left-side
#' spread (median-Q1, median-min) to infer the skew direction.
#'
#' This function is for **informational purposes only** and does not affect
#' the data-generation behaviour of the package.
#'
#' @param q1 First quartile.
#' @param median Median.
#' @param q3 Third quartile.
#' @param min Minimum value.
#' @param max Maximum value.
#' @return Character string: \code{"symmetric"}, \code{"right_skewed"}, or
#'   \code{"left_skewed"}.
#' @keywords internal
diagnose_skewness <- function(q1, median, q3, min, max) {
  left_iqr  <- median - q1
  right_iqr <- q3 - median
  left_rng  <- median - min
  right_rng <- max - median

  # Guard against zero divisions
  if (left_iqr <= 0 || right_iqr <= 0 || left_rng <= 0 || right_rng <= 0)
    return("symmetric")

  r_iqr   <- right_iqr / left_iqr
  r_range <- right_rng / left_rng
  r <- mean(c(r_iqr, r_range))

  if (r > 1.3) return("right_skewed")
  if (r < 0.75) return("left_skewed")
  "symmetric"
}
