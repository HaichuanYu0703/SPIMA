#' Validate Generic Effect-Size Input
#'
#' @param data A data frame with effect-size and SE columns.
#' @param input_spec A named list, e.g.
#'   \code{list(yi = "yi", sei = "sei")}.
#' @return \code{TRUE} invisibly.
#' @export
spima_generic_validate <- function(data, input_spec) {
  required <- c("yi", "sei")
  missing <- setdiff(required, names(input_spec))
  if (length(missing) > 0) {
    stop("input_spec must contain: ", paste(missing, collapse = ", "))
  }

  for (col in unlist(input_spec)) {
    if (!col %in% names(data)) {
      stop("Column '", col, "' not found in data.")
    }
  }

  sei <- data[[input_spec[["sei"]]]]
  if (any(sei <= 0, na.rm = TRUE))
    stop("Standard errors must be positive.")

  invisible(TRUE)
}

#' Simulate Pseudo-Data for Generic Effect-Size
#'
#' No individual-level data is generated. Instead, study-level effect sizes
#' are drawn from the random-effects model:
#'   theta_i ~ N(mu, tau^2)
#'   y_i* ~ N(theta_i, sei_i^2)
#'
#' The returned vector can be treated as "pseudo-IPD" since it directly
#' represents the summary statistics needed for distance computation.
#'
#' @inheritParams spima_bin_simulate
#' @return A named numeric vector of simulated effect sizes (one per study).
#' @export
spima_generic_simulate <- function(study_spec, params, input_spec) {
  mu  <- params["mu"]
  tau <- params["tau"]

  yi_col  <- input_spec[["yi"]]
  sei_col <- input_spec[["sei"]]

  n_studies <- nrow(study_spec)
  sei <- study_spec[[sei_col]]

  # Per-study true effects
  theta_i <- rnorm(n_studies, mean = mu, sd = max(tau, 1e-6))

  # Generate "observed" effect sizes
  yi_sim <- theta_i + rnorm(n_studies, mean = 0, sd = sei)

  # Name with study IDs
  mapped <- resolve_study_col(study_spec, input_spec)
  names(yi_sim) <- mapped$ids

  yi_sim
}

#' Analyze Pseudo-Data for Generic Effect-Size
#'
#' For the generic module, the "pseudo-IPD" is already the summary
#' statistics (a vector of effect sizes). This function simply passes
#' them through with converged = TRUE.
#'
#' @inheritParams spima_bin_analyze
#' @return A list with \code{estimates = NULL},
#'   \code{summary_stats} (the effect-size vector), and
#'   \code{converged = TRUE}.
#' @export
spima_generic_analyze <- function(pseudo_ipd, input_spec) {
  list(
    estimates     = NULL,
    summary_stats = pseudo_ipd,
    converged     = TRUE,
    fit           = NULL
  )
}

#' Compute Observed Summary Statistics for Generic Effect-Size
#'
#' @inheritParams spima_cont_observed_stats
#' @return Named vector of effect sizes with \code{"weights"} attribute
#'   (inverse-variance: 1/sei^2).
#' @export
spima_generic_observed_stats <- function(data, input_spec) {
  yi  <- data[[input_spec[["yi"]]]]
  sei <- data[[input_spec[["sei"]]]]

  mapped <- resolve_study_col(data, input_spec)
  names(yi) <- mapped$ids

  wts <- 1 / (sei^2 + 1e-10)
  names(wts) <- mapped$ids

  attr(yi, "weights") <- wts
  yi
}

#' Distance for Generic Effect-Size
#'
#' Weighted Euclidean distance using inverse-variance weights.
#' @inheritParams spima_bin_distance
#' @export
spima_generic_distance <- function(sim_stats, obs_stats) {
  nm <- intersect(names(sim_stats), names(obs_stats))
  if (length(nm) == 0) return(Inf)

  d <- sim_stats[nm] - obs_stats[nm]
  w <- if (!is.null(attr(obs_stats, "weights"))) {
    attr(obs_stats, "weights")[nm]
  } else {
    1 / (obs_stats[nm]^2 + 1e-8)
  }
  w <- w / sum(w, na.rm = TRUE)
  sqrt(sum(w * d^2, na.rm = TRUE))
}
