#' Validate Binary Outcome Input
#'
#' @param data A data frame with columns for events and sample sizes.
#' @param input_spec A named list specifying column mappings, e.g.
#'   \code{list(event = "event_t", n = "n_t", group = "trt")} for
#'   two-group data, or \code{list(event = "event", n = "n")} for
#'   single-group data.
#' @return \code{TRUE} invisibly; stops with a message on failure.
#' @export
spima_bin_validate <- function(data, input_spec) {
  required <- c("event", "n")
  if (!all(required %in% names(input_spec))) {
    stop("input_spec must contain 'event' and 'n'.")
  }

  for (col in unlist(input_spec)) {
    if (!col %in% names(data)) {
      stop("Column '", col, "' not found in data.")
    }
  }

  # Check counts are plausible
  ev_col <- data[[input_spec[["event"]]]]
  n_col  <- data[[input_spec[["n"]]]]

  if (any(ev_col < 0, na.rm = TRUE))
    stop("Event counts must be non-negative.")
  if (any(n_col <= 0, na.rm = TRUE))
    stop("Sample sizes must be positive.")
  if (any(ev_col > n_col, na.rm = TRUE))
    stop("Event counts cannot exceed sample sizes.")

  invisible(TRUE)
}

#' Simulate Pseudo-IPD for Binary Outcome
#'
#' For each study, the control-group proportion is taken from observed data,
#' and the treatment-group log-odds are shifted by a study-specific effect
#' drawn from N(mu, tau^2). Individual Bernoulli outcomes are then generated.
#'
#' @param study_spec A data frame (subset for one study) containing the
#'   observed counts.
#' @param params Named vector \code{c(mu = ..., tau = ...)}.
#' @param input_spec Column mapping (passed through from \code{spima}).
#' @return A data frame with columns \code{study}, \code{group}, \code{y}.
#' @export
spima_bin_simulate <- function(study_spec, params, input_spec) {
  mu <- params["mu"]
  tau <- params["tau"]

  ev_col  <- input_spec[["event"]]
  n_col   <- input_spec[["n"]]
  grp_col <- input_spec[["group"]]

  # Resolve study column
  mapped <- resolve_study_col(study_spec, input_spec)
  study_spec <- mapped$data
  study_col  <- mapped$col
  studies    <- mapped$ids

  # Check if two-group
  if (is.null(grp_col) || length(unique(study_spec[[grp_col]])) < 2) {
    # ---- Single-group path (R) ----
    all_data <- list()
    for (sid in studies) {
      rows <- study_spec[study_spec[[study_col]] == sid, , drop = FALSE]
      ev <- sum(rows[[ev_col]], na.rm = TRUE)
      n  <- sum(rows[[n_col]], na.rm = TRUE)
      p <- ev / max(n, 1)

      theta_i <- rnorm(1, mean = mu, sd = max(tau, 1e-6))
      logit_p <- qlogis(p + 1e-8) + theta_i
      p_sim <- plogis(logit_p)

      y <- rbinom(n, size = 1, prob = p_sim)
      all_data[[length(all_data) + 1]] <- data.frame(
        study = sid, group = "all", y = y, stringsAsFactors = FALSE
      )
    }
    return(do.call(rbind, all_data))
  }

  # ---- Two-group path (C++) ----
  n_studies <- length(studies)
  n_ctrl_vec <- numeric(n_studies)
  n_trt_vec  <- numeric(n_studies)
  p_ctrl_vec <- numeric(n_studies)

  for (i in seq_len(n_studies)) {
    sid <- studies[i]
    rows <- study_spec[study_spec[[study_col]] == sid, , drop = FALSE]
    glev <- unique(rows[[grp_col]])
    ctrl <- rows[rows[[grp_col]] == glev[1], , drop = FALSE]
    trt  <- rows[rows[[grp_col]] == glev[2], , drop = FALSE]

    ev_c <- sum(ctrl[[ev_col]], na.rm = TRUE)
    n_c  <- sum(ctrl[[n_col]], na.rm = TRUE)
    n_t  <- sum(trt[[n_col]], na.rm = TRUE)

    n_ctrl_vec[i] <- n_c
    n_trt_vec[i]  <- n_t
    p_ctrl_vec[i] <- ev_c / max(n_c, 1)
  }

  theta_i_vec <- rnorm(n_studies, mean = mu, sd = max(tau, 1e-6))

  # C++ generates flat matrix: [study_idx, group(0/1), y]
  mat <- simulate_binary_studies(n_ctrl_vec, n_trt_vec, p_ctrl_vec, theta_i_vec)

  data.frame(
    study = studies[mat[, 1]],
    group = ifelse(mat[, 2] == 0, "control", "treatment"),
    y     = mat[, 3],
    stringsAsFactors = FALSE
  )
}

#' Analyze Pseudo-IPD for Binary Outcome
#'
#' Computes per-study log odds ratios from the simulated pseudo-IPD by
#' constructing 2x2 tables. Also fits a one-stage logistic mixed model
#' (\code{glmer}) for the overall treatment effect estimate.
#'
#' @param pseudo_ipd A data frame from \code{spima_bin_simulate}.
#' @param input_spec Column mapping.
#' @return A list with \code{estimates} (mixed-model fixed effects),
#'   \code{summary_stats} (named vector of per-study log ORs, matching
#'   the format from \code{spima_bin_observed_stats}), and
#'   \code{converged} (logical).
#' @export
spima_bin_analyze <- function(pseudo_ipd, input_spec) {
  grp_col <- input_spec[["group"]]
  has_group <- !is.null(grp_col) &&
               length(unique(pseudo_ipd[["group"]])) > 1

  if (!has_group) {
    # Single-group case: compute log-odds per study
    study_ids <- unique(pseudo_ipd$study)
    logodds <- numeric(length(study_ids))
    for (i in seq_along(study_ids)) {
      rows <- pseudo_ipd[pseudo_ipd$study == study_ids[i], ]
      ev <- sum(rows$y)
      n  <- nrow(rows)
      p  <- (ev + 0.5) / (n + 1)
      logodds[i] <- log(p / (1 - p))
    }
    names(logodds) <- as.character(study_ids)
    return(list(estimates = NULL, summary_stats = logodds,
                converged = TRUE, fit = NULL))
  }

  # Two-group case: compute per-study log ORs like observed_stats does
  study_ids <- unique(pseudo_ipd$study)
  lor <- numeric(length(study_ids))
  names(lor) <- study_ids

  for (i in seq_along(study_ids)) {
    rows <- pseudo_ipd[pseudo_ipd$study == study_ids[i], ]
    groups <- unique(rows$group)
    # Identify control and treatment (use first group level as control)
    # The group column is character; "control" < "treatment" alphabetically
    if (length(groups) < 2) {
      lor[i] <- 0
      next
    }
    # First group = control, second = treatment
    ctrl <- rows[rows$group == groups[1], ]
    trt  <- rows[rows$group == groups[2], ]

    a <- sum(ctrl$y)                        # events in control
    b <- nrow(ctrl) - a                     # non-events in control
    c <- sum(trt$y)                         # events in treatment
    d <- nrow(trt) - c                      # non-events in treatment

    # Log OR with Haldane correction
    lor[i] <- log((c + 0.5) * (b + 0.5) / ((a + 0.5) * (d + 0.5)))
  }

  # Also attempt mixed model for overall estimate (optional; summary_stats
  # uses direct per-study log ORs for distance computation)
  fit <- tryCatch({
    lme4::glmer(y ~ group + (1 | study),
                data = pseudo_ipd,
                family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa",
                                       calc.derivs = FALSE))
  }, error = function(e) NULL)

  glmer_ok <- !is.null(fit) &&
    length(fit@optinfo$conv$lme4$messages) == 0

  list(
    estimates     = if (!is.null(fit)) lme4::fixef(fit) else NULL,
    summary_stats = lor,
    converged     = glmer_ok,
    fit           = fit
  )
}

#' Compute Observed Summary Statistics for Binary Data
#'
#' @param data Original data frame per blueprint.
#' @param input_spec Column mapping.
#' @return Named vector of observed log-ORs (or log-odds) with optional
#'   inverse-variance weights as attribute.
spima_bin_observed_stats <- function(data, input_spec) {
  ev_col  <- input_spec[["event"]]
  n_col   <- input_spec[["n"]]
  grp_col <- input_spec[["group"]]

  mapped <- resolve_study_col(data, input_spec)
  data       <- mapped$data
  study_col  <- mapped$col
  studies    <- mapped$ids

  stats <- numeric(length(studies))
  wts   <- numeric(length(studies))
  names(stats) <- studies

  for (i in seq_along(studies)) {
    sid <- studies[i]
    rows <- data[data[[study_col]] == sid, , drop = FALSE]

    if (!is.null(grp_col) && length(unique(rows[[grp_col]])) >= 2) {
      groups <- unique(rows[[grp_col]])
      ctrl <- rows[rows[[grp_col]] == groups[1], ]
      trt  <- rows[rows[[grp_col]] == groups[2], ]

      a <- sum(ctrl[[ev_col]], na.rm = TRUE)
      b <- sum(ctrl[[n_col]], na.rm = TRUE) - a
      c <- sum(trt[[ev_col]], na.rm = TRUE)
      d <- sum(trt[[n_col]], na.rm = TRUE) - c

      # Log OR with Haldane correction
      lor <- log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))
      se2 <- 1/(a + 0.5) + 1/(b + 0.5) + 1/(c + 0.5) + 1/(d + 0.5)
      stats[i] <- lor
      wts[i]   <- 1 / se2
    } else {
      # Single group: log-odds
      ev <- sum(rows[[ev_col]], na.rm = TRUE)
      n  <- sum(rows[[n_col]], na.rm = TRUE)
      p  <- (ev + 0.5) / (n + 1)
      stats[i] <- log(p / (1 - p))
      wts[i]   <- n  # weight proportional to sample size
    }
  }

  names(wts) <- studies
  attr(stats, "weights") <- wts
  stats
}
