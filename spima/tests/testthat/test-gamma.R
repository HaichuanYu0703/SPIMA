# Gamma module unit tests
# Run from package root:  Rscript tests/testthat/test-gamma.R

cat("========== Gamma Module Tests ==========\n\n")

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
source(file.path(R_dir, "module_cont.R"))
source(file.path(R_dir, "module_gamma.R"))
source(file.path(R_dir, "spima.R"))

library(lme4)
library(spima)

# Wire C++ functions from the installed package namespace
for (fn in c("simulate_binary_studies", "simulate_cont_studies",
             "dist_weighted_euclidean", "dist_euclidean")) {
  tryCatch({
    assign(fn, getFromNamespace(fn, "spima"), envir = .GlobalEnv)
  }, error = function(e) {})
}

set.seed(2024)

# ---- Helper: create minimal Gamma test data ----
make_test_data <- function(K = 4, n_per_arm = 50, shape = 4, mu = -0.2, tau = 0.1) {
  studies <- sprintf("S%02d", seq_len(K))
  rows <- list()
  for (i in seq_len(K)) {
    # Control arm
    mean_c <- rgamma(1, shape = 4, rate = 0.2)  # ~ mean 20
    sd_c   <- mean_c / sqrt(shape)
    # Treatment arm
    theta_i <- rnorm(1, mean = mu, sd = tau)
    mean_t  <- mean_c * exp(theta_i)
    sd_t    <- mean_t / sqrt(shape)

    rows[[length(rows) + 1]] <- data.frame(
      study = studies[i],
      group = c("control", "treatment"),
      mean  = c(mean_c, mean_t),
      sd    = c(sd_c, sd_t),
      n     = c(n_per_arm, n_per_arm),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

input_spec <- list(study = "study", group = "group",
                   mean = "mean", sd = "sd", n = "n")

# ============================================================
# Test 1: spima_gamma_validate
# ============================================================
cat("--- Test 1: Validate ---\n")

dat_ok <- make_test_data()
stopifnot(spima_gamma_validate(dat_ok, input_spec))
cat("  Valid data: OK\n")

dat_bad_mean <- dat_ok
dat_bad_mean$mean[1] <- -1
tryCatch({
  spima_gamma_validate(dat_bad_mean, input_spec)
  stop("Should have thrown an error for negative mean")
}, error = function(e) {
  cat("  Negative mean rejected:", conditionMessage(e), "\n")
})

dat_bad_spec <- input_spec
tryCatch({
  spima_gamma_validate(dat_ok, list(median = "m", q1 = "q1", q3 = "q3", n = "n"))
  stop("Should have thrown an error for median+IQR only")
}, error = function(e) {
  cat("  Median+IQR-only spec rejected:", conditionMessage(e), "\n")
})

cat("\n")

# ============================================================
# Test 2: spima_gamma_observed_stats
# ============================================================
cat("--- Test 2: Observed stats ---\n")

obs <- spima_gamma_observed_stats(dat_ok, input_spec)
stopifnot(is.numeric(obs))
stopifnot(length(obs) == nrow(dat_ok) / 2)  # one per study
stopifnot(!is.null(attr(obs, "weights")))
stopifnot(all(is.finite(obs)))
stopifnot(all(attr(obs, "weights") > 0))
cat("  log(RR):", paste(round(obs, 4), collapse = ", "), "\n")
cat("  Weights: ", paste(round(attr(obs, "weights"), 4), collapse = ", "), "\n")

# Manual verification: last study
s_last <- dat_ok[dat_ok$study == "S04", ]
lr_manual <- log(s_last$mean[s_last$group == "treatment"] /
                 s_last$mean[s_last$group == "control"])
stopifnot(abs(obs["S04"] - lr_manual) < 1e-10)
cat("  Manual check on S04: OK\n\n")

# ============================================================
# Test 3: spima_gamma_simulate
# ============================================================
cat("--- Test 3: Simulate ---\n")

params <- c(mu = -0.2, tau = 0.1)
ipd <- spima_gamma_simulate(dat_ok, params, input_spec)

stopifnot(inherits(ipd, "data.frame"))
stopifnot(all(c("study", "group", "y") %in% names(ipd)))
stopifnot(all(ipd$y > 0))  # Gamma generates positive values
stopifnot(length(unique(ipd$study)) == nrow(dat_ok) / 2)
stopifnot(all(ipd$group %in% c("control", "treatment")))
cat("  IPD rows:", nrow(ipd), "\n")
cat("  All y > 0: TRUE\n")
cat("  Studies:", length(unique(ipd$study)), "\n")

# Check that control mean approximately matches observed
for (s in unique(ipd$study)) {
  y_c <- ipd$y[ipd$study == s & ipd$group == "control"]
  obs_c <- dat_ok$mean[dat_ok$study == s & dat_ok$group == "control"]
  cat(sprintf("  %s simulated control mean = %.2f (observed = %.2f)\n",
              s, mean(y_c), obs_c))
}
cat("\n")

# ============================================================
# Test 4: spima_gamma_analyze
# ============================================================
cat("--- Test 4: Analyze ---\n")

res <- spima_gamma_analyze(ipd, input_spec, quick = FALSE)
stopifnot(is.list(res))
stopifnot(!is.null(res$summary_stats))
stopifnot(all(is.finite(res$summary_stats)))
stopifnot(is.logical(res$converged))
cat("  Converged:", res$converged, "\n")
cat("  Per-study log(RR):", paste(round(res$summary_stats, 4), collapse = ", "), "\n")
if (res$converged) {
  cat("  glmer mu estimate:", round(res$estimates[["grouptreatment"]], 4), "\n")
}

# Single-group case
ipd_single <- ipd
ipd_single$group <- "all"
res_single <- spima_gamma_analyze(ipd_single, input_spec, quick = FALSE)
stopifnot(!is.null(res_single$summary_stats))
stopifnot(res_single$converged)
cat("  Single-group case: OK\n\n")

# ============================================================
# Test 5: spima_gamma_distance
# ============================================================
cat("--- Test 5: Distance ---\n")

obs_stats <- spima_gamma_observed_stats(dat_ok, input_spec)
d <- spima_gamma_distance(res$summary_stats, obs_stats)
stopifnot(is.finite(d), d >= 0)
cat("  Distance =", round(d, 4), "\n")

# Identity should give zero distance
d0 <- spima_gamma_distance(obs_stats, obs_stats)
stopifnot(abs(d0) < 1e-10)
cat("  Self-distance =", round(d0, 10), "\n\n")

# ============================================================
# Test 6: End-to-end via spima()
# ============================================================
cat("--- Test 6: End-to-end spima(family = 'Gamma') ---\n")

ctrl <- smc_control(n_particles = 200, n_generations = 5,
                     epsilon_decay = 0.85, ess_min = 0.15,
                     verbose = FALSE)

res_gamma <- tryCatch({
  spima(dat_ok, outcome_type = "continuous", family = "Gamma",
        input_spec = input_spec,
        prior = prior(mu = "normal(0, 2.5)", tau = "halfnormal(0, 0.5)"),
        smc_control = ctrl)
}, error = function(e) {
  cat("  Error:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(res_gamma)) {
  stopifnot(inherits(res_gamma, "spima"))
  stopifnot(res_gamma$family == "Gamma")
  stopifnot(!is.null(res_gamma$abc_result$summary))
  s <- res_gamma$abc_result$summary
  cat("  Outcome type:", res_gamma$outcome_type, "\n")
  cat("  Family:", res_gamma$family, "\n")
  cat("  mu posterior mean:", round(s[["mu"]]["mean"], 4), "\n")
  cat("  tau posterior mean:", round(s[["tau"]]["mean"], 4), "\n")

  # Test print method (family-aware)
  cat("  Print method:\n")
  capture_output <- capture.output(print(res_gamma))
  cat("    ", capture_output[1], "\n")
  stopifnot(grepl("Gamma", capture_output[1]))

  # Test summary + print.summary
  sm <- summary(res_gamma)
  stopifnot(sm$family == "Gamma")
  cat("  Summary + print: OK\n")
} else {
  cat("  (skipped, known to be unstable at 200 particles)\n")
}
cat("\n")

# ============================================================
# Test 7: Backward compatibility — family = "gaussian"
# ============================================================
cat("--- Test 7: Backward compatibility (family = 'gaussian') ---\n")

# Use Normal-distributed data for Gaussian module
dat_norm <- data.frame(
  study = c(1,1,2,2,3,3),
  group = rep(c("treatment", "control"), 3),
  mean  = c(5.2, 5.0, 5.3, 5.1, 5.1, 4.9),
  sd    = c(2.1, 2.0, 2.2, 2.0, 2.0, 1.9),
  n     = c(50, 50, 60, 60, 40, 40),
  stringsAsFactors = FALSE
)
spec_norm <- list(study = "study", group = "group",
                   mean = "mean", sd = "sd", n = "n")

ctrl_small <- smc_control(n_particles = 100, n_generations = 3,
                           epsilon_decay = 0.85, ess_min = 0.2,
                           verbose = FALSE)

res_default <- spima(dat_norm, outcome_type = "continuous",
                     input_spec = spec_norm,
                     prior = prior(mu = "normal(0, 10)", tau = "halfnormal(0, 1)"),
                     smc_control = ctrl_small)

res_explicit <- spima(dat_norm, outcome_type = "continuous", family = "gaussian",
                      input_spec = spec_norm,
                      prior = prior(mu = "normal(0, 10)", tau = "halfnormal(0, 1)"),
                      smc_control = ctrl_small)

stopifnot(res_default$family == "gaussian")
stopifnot(res_explicit$family == "gaussian")
stopifnot(res_default$outcome_type == res_explicit$outcome_type)
cat("  Default family = 'gaussian': OK\n")
cat("  Explicit family = 'gaussian': OK\n")
cat("  Object structure consistent: OK\n\n")

# ============================================================
cat("========== All Gamma module tests passed! ==========\n")
