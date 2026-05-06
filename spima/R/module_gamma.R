#' Validate Gamma Outcome Input
#'
#' Delegates to \code{spima_cont_validate} (same data format: mean, sd, n
#' per arm) and additionally checks that all means are positive (Gamma
#' distribution is supported on the positive real line).
#'
#' @inheritParams spima_cont_validate
#' @return \code{TRUE} invisibly.
#' @export
spima_gamma_validate <- function(data, input_spec) {
  spima_cont_validate(data, input_spec)

  # Gamma support is (0, Inf)
  for (key in intersect(names(input_spec), c("mean", "median"))) {
    col <- input_spec[[key]]
    if (any(data[[col]] <= 0, na.rm = TRUE)) {
      stop("All ", key, " values must be positive for Gamma family ",
           "(Gamma is supported on the positive real line).")
    }
  }

  invisible(TRUE)
}


#' Simulate Pseudo-IPD for Gamma Outcome
#'
#' For each study, individual data are drawn from a Gamma distribution
#' matching the observed mean and SD via method-of-moments. The treatment
#' group mean is shifted by a multiplicative factor \code{exp(theta_i)}
#' where \code{theta_i ~ N(mu, tau^2)} — this encodes the log-Rate Ratio
#' treatment effect on the original scale.  The shape parameter is held
#' constant within each study, preserving the variance structure implied
#' by the Gamma GLM with log link.
#'
#' @inheritParams spima_cont_simulate
#' @return A data frame with columns \code{study}, \code{group}, \code{y}.
#' @export
spima_gamma_simulate <- function(study_spec, params, input_spec) {
  mu  <- params["mu"]
  tau <- params["tau"]

  if (!all(c("mean", "sd") %in% names(input_spec))) {
    stop("Gamma simulation requires 'mean' and 'sd' in input_spec.")
  }

  n_col   <- input_spec[["n"]]
  grp_col <- input_spec[["group"]]

  mapped <- resolve_study_col(study_spec, input_spec)
  study_spec <- mapped$data
  study_col  <- mapped$col
  studies    <- mapped$ids

  all_data <- list()

  for (sid in studies) {
    rows <- study_spec[study_spec[[study_col]] == sid, , drop = FALSE]
    n_val <- sum(rows[[n_col]], na.rm = TRUE)
    if (n_val <= 0) next

    if (!is.null(grp_col) && length(unique(rows[[grp_col]])) >= 2) {
      groups <- unique(rows[[grp_col]])
      ctrl_rows <- rows[rows[[grp_col]] == groups[1], , drop = FALSE]
      trt_rows  <- rows[rows[[grp_col]] == groups[2], , drop = FALSE]

      mean_c <- mean(ctrl_rows[[input_spec[["mean"]]]], na.rm = TRUE)
      sd_c   <- mean(ctrl_rows[[input_spec[["sd"]]]],   na.rm = TRUE)
      mean_t <- mean(trt_rows[[input_spec[["mean"]]]],  na.rm = TRUE)
      sd_t   <- mean(trt_rows[[input_spec[["sd"]]]],    na.rm = TRUE)

      n_c <- sum(ctrl_rows[[n_col]], na.rm = TRUE)
      n_t <- sum(trt_rows[[n_col]], na.rm = TRUE)

      # Method-of-moments: shape and rate from control arm
      if (mean_c <= 0 || sd_c <= 0 || mean_t <= 0 || sd_t <= 0) next
      shape <- mean_c^2 / sd_c^2
      rate  <- mean_c   / sd_c^2

      # Study-specific log-RR
      theta_i <- rnorm(1, mean = mu, sd = max(tau, 1e-6))

      # Control: Gamma(shape, rate)  →  mean = shape / rate = mean_c
      y_c <- rgamma(n_c, shape = shape, rate = rate)

      # Treatment: adjust rate so that mean = mean_c * exp(theta_i)
      #   mean_t_new = shape / rate_trt = mean_c * exp(theta_i)
      #   =>  rate_trt = rate / exp(theta_i)
      rate_trt <- rate / exp(theta_i)
      y_t <- rgamma(n_t, shape = shape, rate = rate_trt)

      all_data[[length(all_data) + 1]] <- rbind(
        data.frame(study = sid, group = "control",   y = y_c,
                   stringsAsFactors = FALSE),
        data.frame(study = sid, group = "treatment", y = y_t,
                   stringsAsFactors = FALSE)
      )
    } else {
      # Single group: generate from observed mean/SD, no treatment shift
      mean_val <- mean(rows[[input_spec[["mean"]]]], na.rm = TRUE)
      sd_val   <- mean(rows[[input_spec[["sd"]]]],   na.rm = TRUE)
      if (mean_val <= 0 || sd_val <= 0) next
      shape <- mean_val^2 / sd_val^2
      rate  <- mean_val   / sd_val^2
      y <- rgamma(n_val, shape = shape, rate = rate)
      all_data[[length(all_data) + 1]] <- data.frame(
        study = sid, group = "all", y = y, stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, all_data)
}


#' Analyze Pseudo-IPD for Gamma Outcome
#'
#' Fits a one-stage Gamma GLMM via \code{glmer(y ~ group + (1 | study),
#' family = Gamma(link = "log"))} on the simulated pseudo-IPD.  Also
#' returns per-study log-Rate Ratio values for use as summary statistics
#' in the ABC distance computation.
#'
#' Use \code{quick = TRUE} (default) during ABC-SMC sampling where only
#' the per-study summary statistics are needed for distance computation.
#' Set \code{quick = FALSE} to additionally fit the full GLMM (useful for
#' external diagnostics).
#'
#' @inheritParams spima_cont_analyze
#' @param quick If \code{TRUE} (default), skip the GLMM fit and only
#'   compute per-study log-RR summary statistics.
#' @return A list with components:
#'   \describe{
#'     \item{\code{estimates}}{Named vector of fixed effects from the
#'       Gamma GLMM (or \code{NULL} if \code{quick = TRUE} or the model
#'       does not converge).}
#'     \item{\code{summary_stats}}{Named vector of per-study log-RR values.}
#'     \item{\code{converged}}{Logical indicating GLMM convergence
#'       (\code{TRUE} when \code{quick = TRUE}).}
#'     \item{\code{fit}}{The \code{glmer} fit object (or \code{NULL}).}
#'   }
#' @export
spima_gamma_analyze <- function(pseudo_ipd, input_spec, quick = TRUE) {
  has_group <- length(unique(pseudo_ipd[["group"]])) > 1

  if (!has_group) {
    # Single group: per-study log-means
    study_ids <- unique(pseudo_ipd$study)
    log_means <- numeric(length(study_ids))
    for (i in seq_along(study_ids)) {
      y <- pseudo_ipd$y[pseudo_ipd$study == study_ids[i]]
      log_means[i] <- log(mean(y))
    }
    names(log_means) <- as.character(study_ids)
    return(list(
      estimates     = NULL,
      summary_stats = log_means,
      converged     = TRUE,
      fit           = NULL
    ))
  }

  # Two-group case: per-study log-RR (always needed)
  study_ids <- unique(pseudo_ipd$study)
  logRR <- numeric(length(study_ids))
  names(logRR) <- as.character(study_ids)

  for (i in seq_along(study_ids)) {
    rows <- pseudo_ipd[pseudo_ipd$study == study_ids[i], ]
    groups <- unique(rows$group)
    if (length(groups) < 2) {
      logRR[i] <- 0
      next
    }
    y_c <- rows$y[rows$group == groups[1]]
    y_t <- rows$y[rows$group == groups[2]]
    m_c <- mean(y_c)
    m_t <- mean(y_t)
    logRR[i] <- if (m_c > 0 && m_t > 0) log(m_t / m_c) else NA
  }
  logRR <- logRR[is.finite(logRR)]

  if (quick) {
    return(list(
      estimates     = NULL,
      summary_stats = logRR,
      converged     = TRUE,
      fit           = NULL
    ))
  }

  # Full GLMM fit (slow, for external diagnostics only)
  fit <- tryCatch({
    lme4::glmer(y ~ group + (1 | study), data = pseudo_ipd,
                family = Gamma(link = "log"),
                control = lme4::glmerControl(
                  check.conv.grad    = "ignore",
                  check.conv.singular = "ignore"
                ))
  }, error = function(e) NULL)

  converged <- !is.null(fit) &&
    length(fit@optinfo$conv$lme4$messages) == 0

  list(
    estimates     = if (converged) lme4::fixef(fit) else NULL,
    summary_stats = logRR,
    converged     = converged,
    fit           = fit
  )
}


#' Compute Observed Summary Statistics for Gamma Data
#'
#' For each study, computes the observed log-Rate Ratio and its
#' delta-method variance for inverse-variance weighting.
#'
#' @inheritParams spima_cont_observed_stats
#' @param data Original data frame with arm-level means, SDs, and sample
#'   sizes.
#' @param input_spec Column mapping.
#' @return A named vector of per-study log-RR values with an attribute
#'   \code{"weights"} containing inverse-variance weights.
#' @export
spima_gamma_observed_stats <- function(data, input_spec) {
  if (!all(c("mean", "sd") %in% names(input_spec))) {
    stop("Gamma observed stats require 'mean' and 'sd' in input_spec.")
  }

  grp_col <- input_spec[["group"]]

  mapped <- resolve_study_col(data, input_spec)
  data      <- mapped$data
  study_col <- mapped$col
  studies   <- mapped$ids

  if (is.null(grp_col) || length(unique(data[[grp_col]])) < 2) {
    # Single group: return log-mean per study
    log_means <- numeric(length(studies))
    names(log_means) <- studies
    for (i in seq_along(studies)) {
      sid <- studies[i]
      rows <- data[data[[study_col]] == sid, , drop = FALSE]
      m <- mean(rows[[input_spec[["mean"]]]], na.rm = TRUE)
      log_means[i] <- if (m > 0) log(m) else NA
    }
    return(structure(log_means, weights = rep(1, length(log_means))))
  }

  # Two-group case: per-study log-RR with delta-method weights
  logRR <- numeric(length(studies))
  w     <- numeric(length(studies))
  names(logRR) <- names(w) <- studies

  for (i in seq_along(studies)) {
    sid <- studies[i]
    rows <- data[data[[study_col]] == sid, , drop = FALSE]
    groups <- unique(rows[[grp_col]])
    if (length(groups) < 2) {
      logRR[i] <- NA
      w[i]     <- 0
      next
    }
    ctrl <- rows[rows[[grp_col]] == groups[1], ]
    trt  <- rows[rows[[grp_col]] == groups[2], ]

    m_c <- mean(ctrl[[input_spec[["mean"]]]], na.rm = TRUE)
    sd_c <- mean(ctrl[[input_spec[["sd"]]]],   na.rm = TRUE)
    n_c <- sum(ctrl[[input_spec[["n"]]]],      na.rm = TRUE)
    m_t <- mean(trt[[input_spec[["mean"]]]],   na.rm = TRUE)
    sd_t <- mean(trt[[input_spec[["sd"]]]],    na.rm = TRUE)
    n_t <- sum(trt[[input_spec[["n"]]]],       na.rm = TRUE)

    if (m_c <= 0 || m_t <= 0 || n_c < 1 || n_t < 1) {
      logRR[i] <- NA
      w[i] <- 0
      next
    }

    logRR[i] <- log(m_t / m_c)

    # Delta-method variance: var(log(RR)) ≈ sd_c²/(n_c·m_c²) + sd_t²/(n_t·m_t²)
    v <- sd_c^2 / (n_c * m_c^2) + sd_t^2 / (n_t * m_t^2)
    w[i] <- if (is.finite(v) && v > 0) 1 / v else 1
  }

  # Keep only valid entries
  ok <- is.finite(logRR) & is.finite(w)
  logRR <- logRR[ok]
  w <- w[ok]

  structure(logRR, weights = w / sum(w, na.rm = TRUE))
}


#' Distance Function for Gamma Outcome
#'
#' Weighted Euclidean distance on the per-study log-Rate Ratio vector.
#' Delegates to \code{spima_generic_distance}.
#'
#' @inheritParams spima_generic_distance
#' @export
spima_gamma_distance <- function(sim_stats, obs_stats) {
  spima_generic_distance(sim_stats, obs_stats)
}
