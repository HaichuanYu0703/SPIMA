#' Validate Continuous Outcome Input
#'
#' @param data A data frame with means, SDs, and sample sizes.
#' @param input_spec A named list, e.g.
#'   \code{list(mean = "mean_t", sd = "sd_t", n = "n_t", group = "trt")}.
#'   For median+IQR input, use \code{median} and \code{q1}, \code{q3}.
#' @return \code{TRUE} invisibly.
#' @export
spima_cont_validate <- function(data, input_spec) {
  has_mean_sd    <- all(c("mean", "sd") %in% names(input_spec))
  has_median_iqr <- all(c("median", "q1", "q3") %in% names(input_spec))
  has_range      <- all(c("median", "min", "max") %in% names(input_spec))
  has_five       <- all(c("median", "q1", "q3", "min", "max") %in% names(input_spec))

  if (!has_mean_sd && !has_median_iqr && !has_range && !has_five) {
    stop("input_spec must contain 'mean'+'sd', 'median'+'q1'+'q3', ",
         "'median'+'min'+'max' (range), or all five (five-number summary).")
  }

  if (!"n" %in% names(input_spec)) {
    stop("input_spec must contain 'n' (sample size column).")
  }

  # When both mean/SD and quantile columns are present, warn and prefer mean/SD
  if (has_mean_sd && (has_median_iqr || has_range)) {
    message("Both mean/SD and quantile columns detected in input_spec. ",
            "Using mean/SD format.")
  }

  for (col in unlist(input_spec)) {
    if (!col %in% names(data)) {
      stop("Column '", col, "' not found in data.")
    }
  }

  n_col <- data[[input_spec[["n"]]]]
  if (any(n_col <= 0, na.rm = TRUE))
    stop("Sample sizes must be positive.")

  if (has_mean_sd) {
    sd_col <- data[[input_spec[["sd"]]]]
    if (any(sd_col < 0, na.rm = TRUE))
      stop("SD must be non-negative.")
  }

  if (has_median_iqr || has_five) {
    q1 <- data[[input_spec[["q1"]]]]
    q3 <- data[[input_spec[["q3"]]]]
    if (any(q1 > q3, na.rm = TRUE))
      stop("Q1 must be <= Q3.")
  }

  if (has_range || has_five) {
    min_val <- data[[input_spec[["min"]]]]
    max_val <- data[[input_spec[["max"]]]]
    if (any(min_val > max_val, na.rm = TRUE))
      stop("Min must be <= Max.")
    n_vals <- n_col
    if (any(n_vals > 1 & min_val == max_val, na.rm = TRUE))
      stop("Min must be < Max when n > 1.")
  }

  invisible(TRUE)
}

#' Simulate Pseudo-IPD for Continuous Outcome
#'
#' For each study, individual data are drawn from a normal (or skew-normal)
#' distribution matching the observed mean and SD. The treatment group mean
#' is shifted by a study-specific effect drawn from N(mu, tau^2).
#'
#' @inheritParams spima_bin_simulate
#' @return A data frame with columns \code{study}, \code{group}, \code{y}.
#' @export
spima_cont_simulate <- function(study_spec, params, input_spec) {
  mu  <- params["mu"]
  tau <- params["tau"]

  has_mean_sd    <- all(c("mean", "sd") %in% names(input_spec))
  has_median_iqr <- all(c("median", "q1", "q3") %in% names(input_spec))
  has_min_max    <- all(c("min", "max") %in% names(input_spec))
  has_five       <- has_median_iqr && has_min_max

  n_col   <- input_spec[["n"]]
  grp_col <- input_spec[["group"]]

  mapped <- resolve_study_col(study_spec, input_spec)
  study_spec <- mapped$data
  study_col  <- mapped$col
  studies    <- mapped$ids

  # ---- Two-group mean/SD path: C++ accelerated ----
  if (has_mean_sd && !is.null(grp_col) && length(unique(study_spec[[grp_col]])) >= 2) {
    n_studies <- length(studies)
    n_c_vec <- numeric(n_studies)
    n_t_vec <- numeric(n_studies)
    mc_vec  <- numeric(n_studies)
    sc_vec  <- numeric(n_studies)
    st_vec  <- numeric(n_studies)

    for (i in seq_len(n_studies)) {
      sid <- studies[i]
      rows <- study_spec[study_spec[[study_col]] == sid, , drop = FALSE]
      glev <- unique(rows[[grp_col]])
      ctrl_rows <- rows[rows[[grp_col]] == glev[1], , drop = FALSE]
      trt_rows  <- rows[rows[[grp_col]] == glev[2], , drop = FALSE]

      n_c_vec[i] <- sum(ctrl_rows[[n_col]], na.rm = TRUE)
      n_t_vec[i] <- sum(trt_rows[[n_col]], na.rm = TRUE)
      mc_vec[i]  <- mean(ctrl_rows[[input_spec[["mean"]]]], na.rm = TRUE)
      sc_vec[i]  <- mean(ctrl_rows[[input_spec[["sd"]]]],   na.rm = TRUE)
      st_vec[i]  <- mean(trt_rows[[input_spec[["sd"]]]],    na.rm = TRUE)
    }

    theta_i_vec <- rnorm(n_studies, mean = mu, sd = max(tau, 1e-6))
    mat <- simulate_cont_studies(n_c_vec, n_t_vec, mc_vec, sc_vec, st_vec, theta_i_vec)

    return(data.frame(
      study = studies[mat[, 1]],
      group = ifelse(mat[, 2] == 0, "control", "treatment"),
      y     = mat[, 3],
      stringsAsFactors = FALSE
    ))
  }

  # ---- Single-group and quantile-format path (R) ----
  all_data <- list()

  for (sid in studies) {
    rows <- study_spec[study_spec[[study_col]] == sid, , drop = FALSE]
    n_val <- sum(rows[[n_col]], na.rm = TRUE)

    # Extract location (mu) and scale (sigma) from the available format
    if (has_mean_sd) {
      mean_val <- mean(rows[[input_spec[["mean"]]]], na.rm = TRUE)
      sd_val   <- mean(rows[[input_spec[["sd"]]]],   na.rm = TRUE)
    } else {
      median_val <- mean(rows[[input_spec[["median"]]]], na.rm = TRUE)
      if (has_five || has_median_iqr) {
        q1_val <- mean(rows[[input_spec[["q1"]]]], na.rm = TRUE)
        q3_val <- mean(rows[[input_spec[["q3"]]]], na.rm = TRUE)
        sd_val <- (q3_val - q1_val) / 1.35
      } else {
        # Range-only: estimate sigma via Wan et al. 2014
        min_val   <- mean(rows[[input_spec[["min"]]]], na.rm = TRUE)
        max_val   <- mean(rows[[input_spec[["max"]]]], na.rm = TRUE)
        range_val <- max_val - min_val
        p <- max(0.001, min(0.999, (n_val - 0.375) / (n_val + 0.25)))
        sd_val <- range_val / (2 * stats::qnorm(p))
      }
      mean_val <- median_val
    }

    if (n_val <= 0) next

    # Study-specific true mean
    theta_i <- rnorm(1, mean = mu, sd = max(tau, 1e-6))

    if (!is.null(grp_col) && length(unique(rows[[grp_col]])) >= 2) {
      groups <- unique(rows[[grp_col]])
      ctrl_rows <- rows[rows[[grp_col]] == groups[1], , drop = FALSE]
      trt_rows  <- rows[rows[[grp_col]] == groups[2], , drop = FALSE]

      if (has_mean_sd) {
        mean_c <- mean(ctrl_rows[[input_spec[["mean"]]]], na.rm = TRUE)
        sd_c   <- mean(ctrl_rows[[input_spec[["sd"]]]],   na.rm = TRUE)
        mean_t <- mean(trt_rows[[input_spec[["mean"]]]],  na.rm = TRUE)
        sd_t   <- mean(trt_rows[[input_spec[["sd"]]]],    na.rm = TRUE)
      } else {
        med_c <- mean(ctrl_rows[[input_spec[["median"]]]], na.rm = TRUE)
        med_t <- mean(trt_rows[[input_spec[["median"]]]], na.rm = TRUE)
        if (has_five || has_median_iqr) {
          q1_c <- mean(ctrl_rows[[input_spec[["q1"]]]], na.rm = TRUE)
          q3_c <- mean(ctrl_rows[[input_spec[["q3"]]]], na.rm = TRUE)
          sd_c <- (q3_c - q1_c) / 1.35
          q1_t <- mean(trt_rows[[input_spec[["q1"]]]], na.rm = TRUE)
          q3_t <- mean(trt_rows[[input_spec[["q3"]]]], na.rm = TRUE)
          sd_t <- (q3_t - q1_t) / 1.35
        } else {
          # Range-only
          min_c <- mean(ctrl_rows[[input_spec[["min"]]]], na.rm = TRUE)
          max_c <- mean(ctrl_rows[[input_spec[["max"]]]], na.rm = TRUE)
          min_t <- mean(trt_rows[[input_spec[["min"]]]], na.rm = TRUE)
          max_t <- mean(trt_rows[[input_spec[["max"]]]], na.rm = TRUE)
          n_c_arm <- sum(ctrl_rows[[n_col]], na.rm = TRUE)
          n_t_arm <- sum(trt_rows[[n_col]], na.rm = TRUE)
          p_c <- max(0.001, min(0.999, (n_c_arm - 0.375) / (n_c_arm + 0.25)))
          p_t <- max(0.001, min(0.999, (n_t_arm - 0.375) / (n_t_arm + 0.25)))
          sd_c <- (max_c - min_c) / (2 * stats::qnorm(p_c))
          sd_t <- (max_t - min_t) / (2 * stats::qnorm(p_t))
        }
        mean_c <- med_c
        mean_t <- med_t
      }

      n_c <- sum(ctrl_rows[[n_col]], na.rm = TRUE)
      n_t <- sum(trt_rows[[n_col]], na.rm = TRUE)

      # Generate: control ~ observed mean; treatment ~ shifted by theta_i
      y_c <- rnorm(n_c, mean = mean_c, sd = sd_c)
      y_t <- rnorm(n_t, mean = mean_c + theta_i, sd = sd_t)

      all_data[[length(all_data) + 1]] <- rbind(
        data.frame(study = sid, group = "control",    y = y_c, stringsAsFactors = FALSE),
        data.frame(study = sid, group = "treatment",  y = y_t, stringsAsFactors = FALSE)
      )
    } else {
      # Single group
      y <- rnorm(n_val, mean = mean_val, sd = sd_val)
      all_data[[length(all_data) + 1]] <- data.frame(
        study = sid, group = "all", y = y, stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, all_data)
}

#' Analyze Pseudo-IPD for Continuous Outcome
#'
#' Computes per-study mean differences from the simulated pseudo-IPD.
#' Also attempts a linear mixed model (\code{lmer}) for overall estimate.
#'
#' @inheritParams spima_bin_analyze
#' @return A list with \code{estimates}, \code{summary_stats} (per-study
#'   mean differences and pooled SDs, matching observed_stats format),
#'   and \code{converged}.
#' @export
spima_cont_analyze <- function(pseudo_ipd, input_spec) {
  has_group <- length(unique(pseudo_ipd[["group"]])) > 1

  if (!has_group) {
    # Single group: per-study means and SDs
    study_ids <- unique(pseudo_ipd$study)
    means <- sds <- numeric(length(study_ids))
    for (i in seq_along(study_ids)) {
      y <- pseudo_ipd$y[pseudo_ipd$study == study_ids[i]]
      means[i] <- mean(y)
      sds[i]   <- sd(y)
    }
    names(means) <- names(sds) <- as.character(study_ids)
    return(list(
      estimates     = NULL,
      summary_stats = list(means = means, sds = sds),
      converged     = TRUE,
      fit           = NULL
    ))
  }

  # Two-group case: per-study mean difference (treatment - control) + pooled SD
  study_ids <- unique(pseudo_ipd$study)
  mean_diff <- numeric(length(study_ids))
  pooled_sd <- numeric(length(study_ids))

  for (i in seq_along(study_ids)) {
    rows <- pseudo_ipd[pseudo_ipd$study == study_ids[i], ]
    groups <- unique(rows$group)
    if (length(groups) < 2) {
      mean_diff[i] <- 0
      pooled_sd[i] <- sd(rows$y)
      next
    }
    ctrl <- rows$y[rows$group == groups[1]]
    trt  <- rows$y[rows$group == groups[2]]

    if (length(ctrl) < 2 || length(trt) < 2) {
      mean_diff[i] <- mean(trt) - mean(ctrl)
      pooled_sd[i] <- sqrt(var(c(ctrl, trt)))
    } else {
      mean_diff[i] <- mean(trt) - mean(ctrl)
      s1 <- var(ctrl)
      s2 <- var(trt)
      pooled_sd[i] <- sqrt(((length(ctrl) - 1) * s1 +
                            (length(trt) - 1) * s2) /
                           (length(ctrl) + length(trt) - 2))
    }
  }

  names(mean_diff) <- names(pooled_sd) <- as.character(study_ids)

  # Also attempt mixed model
  fit <- tryCatch({
    lme4::lmer(y ~ group + (1 | study), data = pseudo_ipd)
  }, error = function(e) NULL)

  converged <- !is.null(fit) &&
    length(fit@optinfo$conv$lme4$messages) == 0

  list(
    estimates     = if (!is.null(fit)) lme4::fixef(fit) else NULL,
    summary_stats = list(means = mean_diff, sds = pooled_sd),
    converged     = converged,
    fit           = fit
  )
}

#' Compute Observed Summary Statistics for Continuous Data
#'
#' @param data Original data frame.
#' @param input_spec Column mapping.
#' @return A list with \code{means} and \code{sds} (named vectors).
spima_cont_observed_stats <- function(data, input_spec) {
  has_mean_sd   <- all(c("mean", "sd") %in% names(input_spec))
  has_median_iqr <- all(c("median", "q1", "q3") %in% names(input_spec))
  grp_col <- input_spec[["group"]]

  mapped <- resolve_study_col(data, input_spec)
  data      <- mapped$data
  study_col <- mapped$col
  studies   <- mapped$ids

  if (study_col == ".spima_id") {
    # Each row is a separate study
    if (has_mean_sd) {
      means <- setNames(data[[input_spec[["mean"]]]], studies)
      sds   <- setNames(data[[input_spec[["sd"]]]],   studies)
    } else {
      means <- setNames(data[[input_spec[["median"]]]], studies)
      q1    <- data[[input_spec[["q1"]]]]
      q3    <- data[[input_spec[["q3"]]]]
      sds   <- setNames((q3 - q1) / 1.35, studies)
    }
    return(list(means = means, sds = sds))
  }

  means <- sds <- numeric(length(studies))
  names(means) <- names(sds) <- studies

  for (i in seq_along(studies)) {
    sid <- studies[i]
    rows <- data[data[[study_col]] == sid, , drop = FALSE]

    if (has_mean_sd) {
      means[i] <- mean(rows[[input_spec[["mean"]]]], na.rm = TRUE)
      sds[i]   <- mean(rows[[input_spec[["sd"]]]],   na.rm = TRUE)
    } else {
      means[i] <- mean(rows[[input_spec[["median"]]]], na.rm = TRUE)
      q1_val   <- mean(rows[[input_spec[["q1"]]]], na.rm = TRUE)
      q3_val   <- mean(rows[[input_spec[["q3"]]]], na.rm = TRUE)
      sds[i]   <- (q3_val - q1_val) / 1.35
    }

    # Adjust for two groups: compute raw effect size
    if (!is.null(grp_col) && length(unique(rows[[grp_col]])) >= 2) {
      groups <- unique(rows[[grp_col]])
      ctrl <- rows[rows[[grp_col]] == groups[1], ]
      trt  <- rows[rows[[grp_col]] == groups[2], ]

      if (has_mean_sd) {
        means[i] <- mean(trt[[input_spec[["mean"]]]], na.rm = TRUE) -
                    mean(ctrl[[input_spec[["mean"]]]], na.rm = TRUE)
        # Pooled SD
        n1 <- sum(ctrl[[input_spec[["n"]]]], na.rm = TRUE)
        n2 <- sum(trt[[input_spec[["n"]]]], na.rm = TRUE)
        s1 <- mean(ctrl[[input_spec[["sd"]]]], na.rm = TRUE)
        s2 <- mean(trt[[input_spec[["sd"]]]], na.rm = TRUE)
        sds[i] <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
      }
    }
  }

  list(means = means, sds = sds)
}


#' Compute Observed Effect Sizes from Quantile-Format Data
#'
#' Converts median/IQR, median/range, or five-number summary data into
#' per-study effect sizes (median differences) with inverse-variance
#' weights, matching the \code{spima_generic_observed_stats} output format.
#'
#' @param data Original data frame with arm-level quantiles.
#' @param input_spec Column mapping (must contain \code{median},
#'   optionally \code{q1}/\code{q3} and/or \code{min}/\code{max}).
#' @return A named numeric vector of per-study effect sizes with an
#'   \code{attr("weights")} attribute.
#' @keywords internal
spima_cont_quantile_observed_stats <- function(data, input_spec) {
  has_median_iqr <- all(c("median", "q1", "q3") %in% names(input_spec))
  has_min_max    <- all(c("min", "max") %in% names(input_spec))
  has_five       <- has_median_iqr && has_min_max
  grp_col <- input_spec[["group"]]

  mapped <- resolve_study_col(data, input_spec)
  data      <- mapped$data
  study_col <- mapped$col
  studies   <- mapped$ids

  effects <- numeric(length(studies))
  weights <- numeric(length(studies))
  names(effects) <- names(weights) <- studies

  for (i in seq_along(studies)) {
    sid <- studies[i]
    rows <- data[data[[study_col]] == sid, , drop = FALSE]

    if (has_five || has_median_iqr) {
      q1_col <- input_spec[["q1"]]
      q3_col <- input_spec[["q3"]]
    }
    if (has_min_max) {
      min_col <- input_spec[["min"]]
      max_col <- input_spec[["max"]]
    }

    # ---- One-group ----
    if (is.null(grp_col) || length(unique(rows[[grp_col]])) < 2) {
      effects[i] <- mean(rows[[input_spec[["median"]]]], na.rm = TRUE)
      n_arm <- sum(rows[[input_spec[["n"]]]], na.rm = TRUE)
      if (has_five || has_median_iqr) {
        q1v <- mean(rows[[q1_col]], na.rm = TRUE)
        q3v <- mean(rows[[q3_col]], na.rm = TRUE)
        sigma <- (q3v - q1v) / 1.35
      } else {
        minv <- mean(rows[[min_col]], na.rm = TRUE)
        maxv <- mean(rows[[max_col]], na.rm = TRUE)
        p <- max(0.001, min(0.999, (n_arm - 0.375) / (n_arm + 0.25)))
        sigma <- (maxv - minv) / (2 * stats::qnorm(p))
      }
      se <- if (n_arm >= 2 && is.finite(sigma) && sigma > 0)
              sqrt(1.570796 * sigma^2 / n_arm) else NA_real_
      weights[i] <- if (is.finite(se) && se > 0) 1 / se^2 else 0
      next
    }

    # ---- Two-group ----
    groups <- unique(rows[[grp_col]])
    ctrl <- rows[rows[[grp_col]] == groups[1], , drop = FALSE]
    trt  <- rows[rows[[grp_col]] == groups[2], , drop = FALSE]

    med_c <- mean(ctrl[[input_spec[["median"]]]], na.rm = TRUE)
    med_t <- mean(trt[[input_spec[["median"]]]],   na.rm = TRUE)
    n_c   <- sum(ctrl[[input_spec[["n"]]]],        na.rm = TRUE)
    n_t   <- sum(trt[[input_spec[["n"]]]],          na.rm = TRUE)

    # Sigma estimation per arm
    if (has_five || has_median_iqr) {
      q1_c <- mean(ctrl[[q1_col]], na.rm = TRUE)
      q3_c <- mean(ctrl[[q3_col]], na.rm = TRUE)
      q1_t <- mean(trt[[q1_col]],  na.rm = TRUE)
      q3_t <- mean(trt[[q3_col]],  na.rm = TRUE)
      sigma_c <- (q3_c - q1_c) / 1.35
      sigma_t <- (q3_t - q1_t) / 1.35
    } else {
      min_c <- mean(ctrl[[min_col]], na.rm = TRUE)
      max_c <- mean(ctrl[[max_col]], na.rm = TRUE)
      min_t <- mean(trt[[min_col]],  na.rm = TRUE)
      max_t <- mean(trt[[max_col]],  na.rm = TRUE)
      p_c <- max(0.001, min(0.999, (n_c - 0.375) / (n_c + 0.25)))
      p_t <- max(0.001, min(0.999, (n_t - 0.375) / (n_t + 0.25)))
      sigma_c <- (max_c - min_c) / (2 * stats::qnorm(p_c))
      sigma_t <- (max_t - min_t) / (2 * stats::qnorm(p_t))
    }

    # Effect size = median difference
    effects[i] <- med_t - med_c

    # SE(median_diff) = sqrt(pi/2 * sigma_c^2/n_c + pi/2 * sigma_t^2/n_t)
    se_c <- if (n_c >= 2 && is.finite(sigma_c) && sigma_c > 0)
              sqrt(1.570796 * sigma_c^2 / n_c) else NA_real_
    se_t <- if (n_t >= 2 && is.finite(sigma_t) && sigma_t > 0)
              sqrt(1.570796 * sigma_t^2 / n_t) else NA_real_
    se_diff <- if (is.finite(se_c) && is.finite(se_t)) sqrt(se_c^2 + se_t^2) else NA_real_

    weights[i] <- if (is.finite(se_diff) && se_diff > 0) 1 / se_diff^2 else 0
  }

  # Filter valid
  ok <- is.finite(effects) & is.finite(weights) & weights > 0
  effects <- effects[ok]
  w <- weights[ok]

  if (length(effects) == 0) {
    return(structure(numeric(0), weights = numeric(0)))
  }

  structure(effects, weights = w / sum(w, na.rm = TRUE))
}


#' Analyze Pseudo-IPD for Quantile-Format Data
#'
#' Computes per-study median differences (two-group) or per-study medians
#' (one-group) from simulated pseudo-IPD.  Returns a flat named vector of
#' effect sizes matching the output of
#' \code{spima_cont_quantile_observed_stats()} for use with
#' \code{spima_generic_distance()}.
#'
#' @param pseudo_ipd Data frame with columns \code{study}, \code{group},
#'   \code{y} from \code{spima_cont_simulate()}.
#' @param input_spec The original column mapping (only \code{group} is used).
#' @return A list with components \code{estimates} (NULL),
#'   \code{summary_stats} (named numeric vector), \code{converged} (TRUE),
#'   and \code{fit} (NULL).
#' @keywords internal
spima_cont_quantile_analyze <- function(pseudo_ipd, input_spec) {
  grp_col <- input_spec[["group"]]
  has_group <- !is.null(grp_col) &&
               length(unique(pseudo_ipd[[grp_col]])) > 1

  study_ids <- unique(pseudo_ipd$study)
  effects <- numeric(length(study_ids))
  names(effects) <- as.character(study_ids)

  if (!has_group) {
    # One-group: per-study median
    for (i in seq_along(study_ids)) {
      y <- pseudo_ipd$y[pseudo_ipd$study == study_ids[i]]
      effects[i] <- median(y)
    }
  } else {
    # Two-group: per-study median difference
    for (i in seq_along(study_ids)) {
      rows <- pseudo_ipd[pseudo_ipd$study == study_ids[i], ]
      groups <- unique(rows[[grp_col]])
      if (length(groups) < 2) {
        effects[i] <- 0
        next
      }
      y_c <- rows$y[rows[[grp_col]] == groups[1]]
      y_t <- rows$y[rows[[grp_col]] == groups[2]]
      effects[i] <- median(y_t) - median(y_c)
    }
  }

  list(
    estimates     = NULL,
    summary_stats = effects,
    converged     = TRUE,
    fit           = NULL
  )
}
