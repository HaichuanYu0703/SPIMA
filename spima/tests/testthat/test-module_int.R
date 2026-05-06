# Stage V integration test: Interaction Analysis Module
# Run from package root: Rscript tests/testthat/test-module_int.R

cat("========== spima Interaction Module Test ==========\n\n")

# ---- Source dependencies ----
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
source(file.path(R_dir, "module_bin.R"))
source(file.path(R_dir, "module_int.R"))

suppressPackageStartupMessages(library(lme4))
set.seed(2024)

# ---- 1. Test validation ----
cat("\n--- Test 1: Validation ---\n")

dat1 <- data.frame(
  study   = c(1,1,2,2,3,3),
  group   = c(0,1,0,1,0,1),
  event   = c(10, 8, 15, 12, 20, 18),
  n       = c(50, 50, 60, 60, 80, 80),
  mean_X1 = c(45, 47, 52, 53, 60, 58),
  sd_X1   = c(10, 11, 9, 10, 10, 10)
)

spec1 <- list(
  study = "study", event = "event", n = "n",
  group = "group", covariate = "X1"
)
spima_int_validate(dat1, spec1)
cat("  Validation passed (single covariate)\n")

dat2 <- data.frame(
  study   = c(1,1,2,2,3,3),
  group   = c(0,1,0,1,0,1),
  event   = c(10, 8, 15, 12, 20, 18),
  n       = c(50, 50, 60, 60, 80, 80),
  mean_X1 = c(45, 47, 52, 53, 60, 58),
  sd_X1   = c(10, 11, 9, 10, 10, 10),
  mean_X2 = c(130, 128, 135, 136, 140, 142),
  sd_X2   = c(15, 16, 14, 15, 15, 14)
)
spec2 <- list(
  study = "study", event = "event", n = "n",
  group = "group", covariate = c("X1", "X2")
)
spima_int_validate(dat2, spec2)
cat("  Validation passed (two covariates)\n")

caught <- FALSE
tryCatch({
  spima_int_validate(dat1, list(event = "event", n = "n", group = "group"))
}, error = function(e) { caught <<- TRUE })
stopifnot(caught)
cat("  Validation catches missing covariate spec\n")

caught <- FALSE
dat_bad <- dat1
dat_bad$group[1] <- 2
tryCatch({
  spima_int_validate(dat_bad, spec1)
}, error = function(e) { caught <<- TRUE })
stopifnot(caught)
cat("  Validation catches invalid group values\n")

# ---- 2. Test X generation and Y generation utilities ----
cat("\n--- Test 2: X and Y generation utilities ---\n")

set.seed(2024)
K <- 10
n_per_arm <- 50

all_data <- list()
for (j in 1:K) {
  alpha_j <- rnorm(1, -1.5, 0.2)
  n_total <- 2 * n_per_arm
  mu_j <- runif(2, 0, 5)
  X <- cbind(rnorm(n_total, mu_j[1], 2), rnorm(n_total, mu_j[2], 1))
  Y <- rbinom(n_total, 1, plogis(alpha_j + 0.2*X[,1] + 0.3*X[,2]))
  all_data[[j]] <- data.frame(study = j, T = rep(c(0,1), each = n_per_arm),
                              X1 = X[,1], X2 = X[,2], Y = Y)
}
ipd <- do.call(rbind, all_data)

agg_list <- list()
for (j in 1:K) {
  for (t in 0:1) {
    rows <- ipd[ipd$study == j & ipd$T == t, ]
    agg_list[[length(agg_list) + 1]] <- data.frame(
      study = j, T = t,
      event = sum(rows$Y), n = nrow(rows),
      mean_X1 = mean(rows$X1), sd_X1 = max(sd(rows$X1), 0.01),
      mean_X2 = mean(rows$X2), sd_X2 = max(sd(rows$X2), 0.01)
    )
  }
}
agg <- do.call(rbind, agg_list)
spec_agg <- list(
  study = "study", event = "event", n = "n",
  group = "T", covariate = c("X1", "X2")
)

x_test <- .int_gen_x(agg, spec_agg, rho = 0.3)
cat("  Generated", nrow(x_test), "X rows\n")
stopifnot(all(c("study", spec_agg$group, "X1", "X2") %in% names(x_test)))

pseudo_y <- .int_gen_y(x_test, agg, spec_agg, c(X1=0.2, X2=0.3), c(X1=0, X2=0))
stopifnot("Y" %in% names(pseudo_y))
stopifnot(all(pseudo_y$Y %in% c(0, 1)))
cat("  Y generation OK (", sum(pseudo_y$Y), "events in", nrow(pseudo_y), "rows)\n")

# ---- 3. Full spima_int (known interaction scenario) ----
cat("\n--- Test 3: spima_int interaction detection ---\n")

set.seed(2024)
K <- 30
n_per_arm <- 200

# Realistic scenario:
#   X1: biomarker z-score (mean varies 0-4 across studies)
#   gamma: small-moderate (0.2)
#   beta: strong interaction (0.5)
#   Within-study SD = 2 (moderate)

all_data <- list()
for (j in 1:K) {
  alpha_j <- rnorm(1, -1.5, 0.2)
  n_total <- 2 * n_per_arm
  mu_j <- runif(2, 0, 4)  # study-specific means
  X <- cbind(rnorm(n_total, mu_j[1], 2), rnorm(n_total, mu_j[2], 1.5))
  lp <- alpha_j + 0.2*X[,1] + 0.3*X[,2] +
        rep(c(0,1), each = n_per_arm) * 0.5 * X[,1]
  Y <- rbinom(n_total, 1, plogis(lp))
  all_data[[j]] <- data.frame(study = j, T = rep(c(0,1), each = n_per_arm),
                              X1 = X[,1], X2 = X[,2], Y = Y)
}
ipd <- do.call(rbind, all_data)

agg_list <- list()
for (j in 1:K) {
  for (t in 0:1) {
    rows <- ipd[ipd$study == j & ipd$T == t, ]
    agg_list[[length(agg_list) + 1]] <- data.frame(
      study = j, T = t,
      event = sum(rows$Y), n = nrow(rows),
      mean_X1 = mean(rows$X1), sd_X1 = max(sd(rows$X1), 0.01),
      mean_X2 = mean(rows$X2), sd_X2 = max(sd(rows$X2), 0.01)
    )
  }
}
agg <- do.call(rbind, agg_list)

spec_agg <- list(
  study = "study", event = "event", n = "n",
  group = "T", covariate = c("X1", "X2")
)

# IPD truth for verification
ipd_fit <- glmer(Y ~ T * X1 + T * X2 + (1 | study),
                 data = ipd, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))
ipd_beta <- fixef(ipd_fit)["T:X1"]
cat(sprintf("  IPD Truth T:X1 = %.4f (true = 0.5)\n", ipd_beta))

# SPI-MA interaction analysis
res <- spima_int(agg, spec_agg, rho = 0)

stopifnot(inherits(res, "spima_int"))
cat("  Aggregate model converged:", res$converged, "\n")
cat("  Pseudo-IPD converged:", res$pseudo_converged, "\n")

if (res$converged) {
  int_tab <- res$coefficients[res$coefficients$is_interaction, ]
  for (i in seq_len(nrow(int_tab))) {
    cat(sprintf("  %s: est = %.4f, SE = %.4f, p = %.4f\n",
                int_tab$term[i], int_tab$estimate[i],
                int_tab$se[i], int_tab$p_value[i]))
  }

  # The aggregate interaction should be positive and significant
  agg_x1 <- int_tab$estimate[grep("mean_X1", int_tab$term)]
  agg_p  <- int_tab$p_value[grep("mean_X1", int_tab$term)]
  if (length(agg_x1) > 0) {
    cat(sprintf("  Aggregate T:mean_X1 = %.4f, p = %.4f\n", agg_x1[1], agg_p[1]))
    # Should detect positive interaction (even with attenuation)
    stopifnot(agg_x1[1] > 0)
    stopifnot(agg_p[1] < 0.05)
  }
}

# ---- 4. S3 methods ----
cat("\n--- Test 4: S3 methods ---\n")

captured <- capture.output(print(res))
cat("  print() OK (", length(captured), "lines)\n")

tab <- as.data.frame(res)
stopifnot(inherits(tab, "data.frame"))
stopifnot("term" %in% names(tab))
cat("  as.data.frame() OK (", nrow(tab), "rows)\n")

cat("\n========== All Interaction Module Tests Passed ==========\n")
