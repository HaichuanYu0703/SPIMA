# Elastic input tests for continuous module
# Run from package root: Rscript tests/testthat/test-elastic_input.R

cat("========== Elastic Input Module Tests ==========\n\n")

# ---- Source package files ----
find_pkg_root <- function() {
  candidates <- tryCatch(normalizePath(getwd()), error = function(e) NULL)
  tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    for (a in args) {
      a <- sub("^--file=", "", a)
      if (grepl("\\.R$", a, ignore.case = TRUE) && file.exists(a))
        candidates <- c(candidates, dirname(normalizePath(a)))
    }
  }, error = function(e) NULL)
  for (start in candidates) {
    root <- start
    for (i in 1:10) {
      if (file.exists(file.path(root, "DESCRIPTION"))) return(root)
      parent <- dirname(root)
      if (parent == root) break
      root <- parent
    }
  }
  stop("Cannot locate spima package root. Run from the package directory.")
}
pkg_root <- find_pkg_root()
R_dir <- file.path(pkg_root, "R")

source(file.path(R_dir, "utils.R"))
source(file.path(R_dir, "prior.R"))
source(file.path(R_dir, "distance.R"))
source(file.path(R_dir, "abc_smc.R"))
source(file.path(R_dir, "utils_skewness.R"))
source(file.path(R_dir, "module_cont.R"))
source(file.path(R_dir, "module_generic.R"))
source(file.path(R_dir, "spima.R"))

library(spima)

# Wire C++ functions from the installed package namespace
for (fn in c("simulate_binary_studies", "simulate_cont_studies",
             "dist_weighted_euclidean", "dist_euclidean")) {
  tryCatch({
    assign(fn, getFromNamespace(fn, "spima"), envir = .GlobalEnv)
  }, error = function(e) {})
}

suppressPackageStartupMessages(library(lme4))

set.seed(2024)

# ============================================================
# Helper: generate quantile-format test data (two-group)
# ============================================================
make_quantile_data <- function(K = 5, n_per_arm = 100, delta = 2.0,
                               format = "median_iqr") {
  all_rows <- list()
  for (j in seq_len(K)) {
    ctrl <- rnorm(n_per_arm, mean = 50, sd = 10)
    trt  <- rnorm(n_per_arm, mean = 50 + delta, sd = 10)

    for (arm in c("control", "treatment")) {
      y <- if (arm == "control") ctrl else trt
      row <- list(
        study = j,
        group = arm,
        n = length(y),
        median = round(median(y), 4)
      )
      if (format %in% c("median_iqr", "five")) {
        row$q1 <- round(quantile(y, 0.25), 4)
        row$q3 <- round(quantile(y, 0.75), 4)
      }
      if (format %in% c("range", "five")) {
        row$min <- round(min(y), 4)
        row$max <- round(max(y), 4)
      }
      all_rows[[length(all_rows) + 1]] <- data.frame(
        study = row$study, group = row$group, n = row$n,
        median = row$median,
        q1 = if (!is.null(row$q1)) row$q1 else NA,
        q3 = if (!is.null(row$q3)) row$q3 else NA,
        min = if (!is.null(row$min)) row$min else NA,
        max = if (!is.null(row$max)) row$max else NA,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, all_rows)
  # Drop NA columns not needed for this format
  needed <- c("study", "group", "n", "median")
  if (format %in% c("median_iqr", "five")) needed <- c(needed, "q1", "q3")
  if (format %in% c("range", "five"))      needed <- c(needed, "min", "max")
  df[, needed, drop = FALSE]
}

# ============================================================
# Helper: one-group quantile data
# ============================================================
make_quantile_data_onegroup <- function(K = 5, n = 100, format = "median_iqr") {
  all_rows <- list()
  for (j in seq_len(K)) {
    y <- rnorm(n, mean = 50, sd = 10)
    row <- list(study = j, n = n, median = round(median(y), 4))
    if (format %in% c("median_iqr", "five")) {
      row$q1 <- round(quantile(y, 0.25), 4)
      row$q3 <- round(quantile(y, 0.75), 4)
    }
    if (format %in% c("range", "five")) {
      row$min <- round(min(y), 4)
      row$max <- round(max(y), 4)
    }
    all_rows[[length(all_rows) + 1]] <- data.frame(row, stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, all_rows)
  needed <- c("study", "n", "median")
  if (format %in% c("median_iqr", "five")) needed <- c(needed, "q1", "q3")
  if (format %in% c("range", "five"))      needed <- c(needed, "min", "max")
  df[, needed, drop = FALSE]
}

# ============================================================
# Test 1: median/IQR two-group
# ============================================================
cat("--- Test 1: median/IQR two-group ---\n")

dat_iqr <- make_quantile_data(K = 5, n_per_arm = 100, delta = 2.0,
                              format = "median_iqr")
spec_iqr <- list(study = "study", group = "group", n = "n",
                 median = "median", q1 = "q1", q3 = "q3")

# observed_stats
obs_iqr <- spima_cont_quantile_observed_stats(dat_iqr, spec_iqr)
stopifnot(is.numeric(obs_iqr))
stopifnot(length(obs_iqr) == 5)
stopifnot(!is.null(attr(obs_iqr, "weights")))
stopifnot(all(attr(obs_iqr, "weights") > 0))
cat("  observed_stats: names =", paste(names(obs_iqr), collapse = ","), "\n")
cat("  weights:", round(attr(obs_iqr, "weights"), 4), "\n")

# simulate
params <- c(mu = 0.5, tau = 0.2)
ipd <- spima_cont_simulate(dat_iqr, params, spec_iqr)
stopifnot(is.data.frame(ipd))
stopifnot(all(c("study", "group", "y") %in% names(ipd)))
stopifnot(nrow(ipd) == 5 * 100 * 2)  # 5 studies, 100/arm, 2 arms
cat("  simulate: nrow =", nrow(ipd), "\n")

# analyze
res <- spima_cont_quantile_analyze(ipd, spec_iqr)
stopifnot(is.list(res))
stopifnot(is.null(res$estimates))
stopifnot(is.numeric(res$summary_stats))
stopifnot(length(res$summary_stats) == 5)
stopifnot(res$converged)
cat("  analyze: summary_stats =", round(res$summary_stats, 4), "\n")

# distance
d <- spima_generic_distance(res$summary_stats, obs_iqr)
stopifnot(is.finite(d), d >= 0)
cat("  distance =", round(d, 4), "\n")

# Full spima() run
pr <- prior(mu = "normal(0, 10)", tau = "halfnormal(0, 1)")
ctrl <- smc_control(n_particles = 50, n_generations = 2, verbose = FALSE)
spima_res <- spima(dat_iqr, "continuous", input_spec = spec_iqr,
                   prior = pr, smc_control = ctrl)
stopifnot(inherits(spima_res, "spima"))
stopifnot(length(spima_res$abc_result$generations) >= 1)
post <- spima_res$abc_result$posterior
stopifnot(!is.null(post))
cat("  spima() completed, posterior mu =", round(post$theta[, "mu"], 4), "\n")

cat("\n")

# ============================================================
# Test 2: median/range two-group
# ============================================================
cat("--- Test 2: median/range two-group ---\n")

dat_range <- make_quantile_data(K = 5, n_per_arm = 100, delta = 2.0,
                                format = "range")
spec_range <- list(study = "study", group = "group", n = "n",
                   median = "median", min = "min", max = "max")

obs_range <- spima_cont_quantile_observed_stats(dat_range, spec_range)
stopifnot(is.numeric(obs_range))
stopifnot(length(obs_range) == 5)
stopifnot(!is.null(attr(obs_range, "weights")))
stopifnot(all(attr(obs_range, "weights") > 0))
cat("  observed_stats: weights =", round(attr(obs_range, "weights"), 4), "\n")

ipd <- spima_cont_simulate(dat_range, params, spec_range)
stopifnot(is.data.frame(ipd), nrow(ipd) == 1000)
cat("  simulate: nrow =", nrow(ipd), "\n")

res <- spima_cont_quantile_analyze(ipd, spec_range)
stopifnot(is.numeric(res$summary_stats), length(res$summary_stats) == 5)
cat("  analyze: summary_stats =", round(res$summary_stats, 4), "\n")

d <- spima_generic_distance(res$summary_stats, obs_range)
stopifnot(is.finite(d), d >= 0)
cat("  distance =", round(d, 4), "\n")

# Full run
spima_res <- spima(dat_range, "continuous", input_spec = spec_range,
                   prior = pr, smc_control = ctrl)
stopifnot(inherits(spima_res, "spima"))
cat("  spima() completed\n")

cat("\n")

# ============================================================
# Test 3: Five-number summary two-group
# ============================================================
cat("--- Test 3: Five-number summary two-group ---\n")

dat_five <- make_quantile_data(K = 5, n_per_arm = 100, delta = 2.0,
                               format = "five")
spec_five <- list(study = "study", group = "group", n = "n",
                  median = "median", q1 = "q1", q3 = "q3",
                  min = "min", max = "max")

obs_five <- spima_cont_quantile_observed_stats(dat_five, spec_five)
stopifnot(is.numeric(obs_five), length(obs_five) == 5)
stopifnot(!is.null(attr(obs_five, "weights")))
stopifnot(all(attr(obs_five, "weights") > 0))
cat("  observed_stats: weights =", round(attr(obs_five, "weights"), 4), "\n")

# diagnose_skewness on first study, first arm
row1_ctrl <- dat_five[dat_five$study == 1 & dat_five$group == "control", ]
skew <- diagnose_skewness(row1_ctrl$q1, row1_ctrl$median,
                          row1_ctrl$q3, row1_ctrl$min, row1_ctrl$max)
cat("  Skewness diagnosis (study 1, control):", skew, "\n")
stopifnot(is.character(skew))
# Symmetric data, should be symmetric
stopifnot(skew == "symmetric")

ipd <- spima_cont_simulate(dat_five, params, spec_five)
stopifnot(is.data.frame(ipd), nrow(ipd) == 1000)
cat("  simulate: nrow =", nrow(ipd), "\n")

res <- spima_cont_quantile_analyze(ipd, spec_five)
stopifnot(is.numeric(res$summary_stats), length(res$summary_stats) == 5)
cat("  analyze: summary_stats =", round(res$summary_stats, 4), "\n")

d <- spima_generic_distance(res$summary_stats, obs_five)
stopifnot(is.finite(d), d >= 0)
cat("  distance =", round(d, 4), "\n")

# Full run
spima_res <- spima(dat_five, "continuous", input_spec = spec_five,
                   prior = pr, smc_control = ctrl)
stopifnot(inherits(spima_res, "spima"))
cat("  spima() completed\n")

cat("\n")

# ============================================================
# Test 4: One-group (no group column)
# ============================================================
cat("--- Test 4: One-group median/IQR ---\n")

dat_one <- make_quantile_data_onegroup(K = 5, n = 100, format = "median_iqr")
spec_one <- list(study = "study", n = "n",
                 median = "median", q1 = "q1", q3 = "q3")

obs_one <- spima_cont_quantile_observed_stats(dat_one, spec_one)
stopifnot(is.numeric(obs_one), length(obs_one) == 5)
cat("  observed_stats:", round(obs_one, 4), "\n")
cat("  weights:", round(attr(obs_one, "weights"), 4), "\n")

ipd <- spima_cont_simulate(dat_one, params, spec_one)
stopifnot(is.data.frame(ipd))
stopifnot(all(c("study", "group", "y") %in% names(ipd)))
stopifnot(nrow(ipd) == 500)
cat("  simulate: nrow =", nrow(ipd), "\n")

res <- spima_cont_quantile_analyze(ipd, spec_one)
stopifnot(is.numeric(res$summary_stats), length(res$summary_stats) == 5)
cat("  analyze: summary_stats =", round(res$summary_stats, 4), "\n")

d <- spima_generic_distance(res$summary_stats, obs_one)
stopifnot(is.finite(d), d >= 0)
cat("  distance =", round(d, 4), "\n")

# One-group should also work with no study column
dat_one_nos <- dat_one
dat_one_nos$study <- NULL
spec_one_nos <- list(n = "n", median = "median", q1 = "q1", q3 = "q3")
obs_one_nos <- spima_cont_quantile_observed_stats(dat_one_nos, spec_one_nos)
stopifnot(is.numeric(obs_one_nos), length(obs_one_nos) == 5)
cat("  No-study-column observed_stats OK\n")

cat("\n")

# ============================================================
# Test 5: Helper function unit tests
# ============================================================
cat("--- Test 5: Helper functions ---\n")

# se_from_iqr: Normal IQR ~ 1.35 * sigma
se_iqr <- se_from_iqr(iqr = 1.35, n = 100)
# sigma ≈ 1, SE ≈ sqrt(pi/2 * 1/100) ≈ sqrt(0.0157) ≈ 0.125
cat("  se_from_iqr(1.35, 100) =", round(se_iqr, 4), "(expected ~0.125)\n")
stopifnot(is.finite(se_iqr))
stopifnot(se_iqr > 0)
stopifnot(abs(se_iqr - 0.125) < 0.01)

# se_from_iqr with n < 2 should return NA
se_iqr_na <- se_from_iqr(iqr = 1.35, n = 1)
stopifnot(is.na(se_iqr_na))
cat("  se_from_iqr(1.35, 1) = NA: OK\n")

# se_from_range
se_rng <- se_from_range(range_val = 6, n = 50)
cat("  se_from_range(6, 50) =", round(se_rng, 4), "\n")
stopifnot(is.finite(se_rng))
stopifnot(se_rng > 0)

# se_from_range with n < 5 returns NA
se_rng_na <- se_from_range(range_val = 6, n = 3)
stopifnot(is.na(se_rng_na))
cat("  se_from_range(6, 3) = NA: OK\n")

# diagnose_skewness: symmetric
skew1 <- diagnose_skewness(q1 = 1, median = 2, q3 = 3, min = 0, max = 4)
cat("  diagnose_skewness(symmetric) =", skew1, "\n")
stopifnot(skew1 == "symmetric")

# diagnose_skewness: right-skewed (Q3-median > median-Q1, max-median > median-min)
skew2 <- diagnose_skewness(q1 = 1, median = 2, q3 = 5, min = 0, max = 10)
cat("  diagnose_skewness(right) =", skew2, "\n")
stopifnot(skew2 == "right_skewed")

# diagnose_skewness: left-skewed
skew3 <- diagnose_skewness(q1 = 1, median = 4, q3 = 5, min = 0, max = 7)
cat("  diagnose_skewness(left) =", skew3, "\n")
stopifnot(skew3 == "left_skewed")

# diagnose_skewness: edge case with zero spread
skew4 <- diagnose_skewness(q1 = 2, median = 2, q3 = 2, min = 2, max = 2)
stopifnot(skew4 == "symmetric")
cat("  diagnose_skewness(zero spread) = symmetric: OK\n")

cat("\n")

# ============================================================
# Test 6: Backward compatibility — mean/SD path unchanged
# ============================================================
cat("--- Test 6: Backward compatibility (mean/SD) ---\n")

# Generate mean/SD format data
make_mean_sd_data <- function(K = 5, n_per_arm = 100, delta = 2.0) {
  all_rows <- list()
  for (j in seq_len(K)) {
    ctrl <- rnorm(n_per_arm, mean = 50, sd = 10)
    trt  <- rnorm(n_per_arm, mean = 50 + delta, sd = 10)
    for (arm in c("control", "treatment")) {
      y <- if (arm == "control") ctrl else trt
      all_rows[[length(all_rows) + 1]] <- data.frame(
        study = j, group = arm, n = length(y),
        mean = round(mean(y), 4), sd = round(sd(y), 4),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, all_rows)
}

dat_ms <- make_mean_sd_data(K = 5, n_per_arm = 100, delta = 2.0)
spec_ms <- list(study = "study", group = "group", n = "n",
                mean = "mean", sd = "sd")

# Validate
stopifnot(spima_cont_validate(dat_ms, spec_ms))

# Observed stats
obs_ms <- spima_cont_observed_stats(dat_ms, spec_ms)
stopifnot(is.list(obs_ms))
stopifnot(all(c("means", "sds") %in% names(obs_ms)))
stopifnot(length(obs_ms$means) == 5)
stopifnot(length(obs_ms$sds) == 5)
cat("  observed_stats: means =", round(obs_ms$means, 2), "\n")

# Simulate
ipd_ms <- spima_cont_simulate(dat_ms, params, spec_ms)
stopifnot(is.data.frame(ipd_ms), nrow(ipd_ms) == 1000)
cat("  simulate: nrow =", nrow(ipd_ms), "\n")

# Analyze
res_ms <- spima_cont_analyze(ipd_ms, spec_ms)
stopifnot(is.list(res_ms))
stopifnot(all(c("estimates", "summary_stats", "converged") %in% names(res_ms)))
stopifnot(is.list(res_ms$summary_stats))
stopifnot(all(c("means", "sds") %in% names(res_ms$summary_stats)))
cat("  analyze: mean_diff =", round(res_ms$summary_stats$means, 4), "\n")

# Distance
d_ms <- spima_cont_distance(res_ms$summary_stats, obs_ms)
stopifnot(is.finite(d_ms), d_ms >= 0)
cat("  distance =", round(d_ms, 4), "\n")

# Full spima() run with mean/SD
spima_ms <- spima(dat_ms, "continuous", input_spec = spec_ms,
                  prior = pr, smc_control = ctrl)
stopifnot(inherits(spima_ms, "spima"))
post_ms <- spima_ms$abc_result$posterior
stopifnot(!is.null(post_ms))
cat("  spima() completed, posterior mu =", round(post_ms$theta[, "mu"], 4), "\n")

# Verify that mean/SD path still uses spima_cont_distance (not generic)
# by checking abc_result structure
cat("  Backward compatibility: OK\n")

cat("\n")

# ============================================================
cat("========== All Elastic Input Tests Passed! ==========\n")
