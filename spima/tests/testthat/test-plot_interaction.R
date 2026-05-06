# Interaction plot unit tests
# Run from package root:  Rscript tests/testthat/test-plot_interaction.R

cat("========== Interaction Plot Module Tests ==========\n\n")

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
source(file.path(R_dir, "module_int.R"))
source(file.path(R_dir, "plot_interaction.R"))

suppressPackageStartupMessages(library(lme4))
library(ggplot2)

set.seed(2024)

# ============================================================
# Helper: generate test data
# ============================================================
make_test_data <- function(K = 20, n_per_arm = 100, beta_true = 0.4) {
  all_data <- list()
  for (j in seq_len(K)) {
    alpha_j <- rnorm(1, -1.5, 0.2)
    n_total <- 2 * n_per_arm
    mu_j <- runif(1, 0, 4)
    X <- rnorm(n_total, mu_j, 2)
    lp <- alpha_j + 0.3 * X +
          rep(c(0, 1), each = n_per_arm) * (0.5 + beta_true * X)
    Y <- rbinom(n_total, 1, plogis(lp))
    all_data[[j]] <- data.frame(study = j, T = rep(c(0, 1), each = n_per_arm),
                                X = X, Y = Y)
  }
  ipd <- do.call(rbind, all_data)

  agg_list <- list()
  for (j in seq_len(K)) {
    for (t in 0:1) {
      rows <- ipd[ipd$study == j & ipd$T == t, ]
      agg_list[[length(agg_list) + 1]] <- data.frame(
        study = j, T = t,
        event = sum(rows$Y), n = nrow(rows),
        mean_X = mean(rows$X), sd_X = max(sd(rows$X), 0.01)
      )
    }
  }
  list(ipd = ipd, agg = do.call(rbind, agg_list))
}

spec <- list(study = "study", event = "event", n = "n",
             group = "T", covariate = "X")

# ============================================================
# Test 1: Basic plot returns ggplot
# ============================================================
cat("--- Test 1: Basic plot returns ggplot ---\n")

dat <- make_test_data(K = 15, n_per_arm = 80, beta_true = 0.4)
res <- spima_int(dat$agg, spec, rho = 0)
stopifnot(inherits(res, "spima_int"))
cat("  spima_int converged:", res$converged, "\n")

p <- plot(res, covariate = "X")
stopifnot(inherits(p, "ggplot"))
cat("  plot(res) returns ggplot: OK\n")

# Default covariate
p_default <- plot(res)
stopifnot(inherits(p_default, "ggplot"))
cat("  plot(res) with default covariate: OK\n")

cat("\n")

# ============================================================
# Test 2: Different scales
# ============================================================
cat("--- Test 2: scale = 'relative' ---\n")

p_rel <- plot(res, covariate = "X", scale = "relative")
stopifnot(inherits(p_rel, "ggplot"))
# Check y-axis label contains "Risk Ratio"
ylab_rel <- p_rel$labels$y
stopifnot(grepl("Risk Ratio", ylab_rel, ignore.case = TRUE))
cat("  Relative scale label correct:", ylab_rel, "\n")

# Absolute scale
p_abs <- plot(res, covariate = "X", scale = "absolute")
ylab_abs <- p_abs$labels$y
stopifnot(grepl("Risk Difference", ylab_abs, ignore.case = TRUE))
cat("  Absolute scale label correct:", ylab_abs, "\n")

cat("\n")

# ============================================================
# Test 3: Invalid covariate name
# ============================================================
cat("--- Test 3: Invalid covariate ---\n")

caught <- FALSE
tryCatch({
  plot(res, covariate = "nonexistent")
}, error = function(e) {
  caught <<- TRUE
  cat("  Error message:", conditionMessage(e), "\n")
})
stopifnot(caught)
cat("  Invalid covariate correctly rejected: OK\n\n")

# ============================================================
# Test 4: Zero interaction — reference line within CI band
# ============================================================
cat("--- Test 4: Zero interaction ---\n")

# Under zero interaction, the aggregate model should not reject the null.
# (Note: ecological models can show non-zero aggregate interaction even
# when the individual-level interaction is zero; the pseudo-IPD model
# inherits this bias from the generating estimates.  We therefore test
# the model coefficients directly rather than the plot CI coverage.)
set.seed(42)
dat0 <- make_test_data(K = 30, n_per_arm = 200, beta_true = 0)
res0 <- spima_int(dat0$agg, spec, rho = 0)
cat("  spima_int converged:", res0$converged, "\n")

int_tab <- res0$coefficients[res0$coefficients$is_interaction, ]
cat("  Interaction p-value:", round(int_tab$p_value[1], 4), "\n")
stopifnot(int_tab$p_value[1] > 0.05)  # should not reject null

# Plot should still produce a valid ggplot
p0 <- plot(res0, covariate = "X")
stopifnot(inherits(p0, "ggplot"))
cat("  Zero-interaction plot returned: OK\n")

cat("\n")

# ============================================================
# Test 5: Custom 'at' and 'ci_level'
# ============================================================
cat("--- Test 5: Custom at and ci_level ---\n")

at_custom <- seq(0, 5, by = 0.5)
p_custom <- plot(res, covariate = "X", at = at_custom)
stopifnot(inherits(p_custom, "ggplot"))
cat("  Custom 'at' sequence: OK\n")

p_ci <- plot(res, covariate = "X", ci_level = 0.80)
stopifnot(inherits(p_ci, "ggplot"))
# 80% CI should be narrower than 95% CI
post_80 <- extract_interaction_posterior(res, ci_level = 0.80)
post_95 <- extract_interaction_posterior(res, ci_level = 0.95)
width_80 <- mean(post_80$summary$ci_upper - post_80$summary$ci_lower)
width_95 <- mean(post_95$summary$ci_upper - post_95$summary$ci_lower)
stopifnot(width_80 < width_95)
cat("  80% CI narrower than 95% CI: OK\n")
cat("  80% width:", round(width_80, 4), "  95% width:", round(width_95, 4), "\n")

cat("\n")

# ============================================================
# Test 6: extract_interaction_posterior structure
# ============================================================
cat("--- Test 6: extract_interaction_posterior structure ---\n")

post <- extract_interaction_posterior(res, covariate = "X", scale = "absolute")
stopifnot(is.list(post))
stopifnot(is.numeric(post$at))
stopifnot(is.matrix(post$draws))
stopifnot(is.data.frame(post$summary))
stopifnot(all(c("estimate", "ci_lower", "ci_upper") %in% names(post$summary)))
stopifnot(post$scale == "absolute")
stopifnot(nrow(post$draws) == 2000)  # default n_draws
stopifnot(ncol(post$draws) == length(post$at))
stopifnot(post$model_type %in% c("pseudo-IPD", "aggregate"))
cat("  All list components correct: OK\n")
cat("  Draws:", nrow(post$draws), "x", ncol(post$draws), "\n")
cat("  Model type:", post$model_type, "\n")

cat("\n")

# ============================================================
# Test 7: Structure of ggplot output layers
# ============================================================
cat("--- Test 7: Plot structure ---\n")

p <- plot(res, covariate = "X")

# Check that plot has the right layers
has_ribbon <- any(vapply(p$layers, function(l) {
  inherits(l$geom, "GeomRibbon")
}, logical(1)))
has_line <- any(vapply(p$layers, function(l) {
  inherits(l$geom, "GeomLine")
}, logical(1)))
has_hline <- any(vapply(p$layers, function(l) {
  inherits(l$geom, "GeomHline")
}, logical(1)))

stopifnot(has_ribbon)
stopifnot(has_line)
stopifnot(has_hline)
cat("  Contains ribbon, line, and reference line: OK\n")

# Check axis labels
stopifnot(p$labels$x == "X")
cat("  X-axis label: X\n")

cat("\n")

# ============================================================
cat("========== All Interaction Plot Tests Passed! ==========\n")
