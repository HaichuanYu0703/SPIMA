#' Distance Functions for ABC-SMC
#'
#' Each module exports a distance function that compares simulated
#' summary statistics to the observed summary statistics.
#'
#' @param sim_stats Simulated summary statistics (vector or list).
#' @param obs_stats Observed summary statistics (same structure).
#' @return A non-negative scalar distance.
#' @name distance_functions
NULL

#' @describeIn distance_functions Binary outcome: Euclidean distance on
#'   (possibly weighted) log-odds scale.
spima_bin_distance <- function(sim_stats, obs_stats) {
  nm <- intersect(names(sim_stats), names(obs_stats))
  if (length(nm) == 0) return(Inf)
  w <- if (!is.null(attr(obs_stats, "weights"))) {
    attr(obs_stats, "weights")[nm]
  } else {
    rep(1, length(nm))
  }
  dist_weighted_euclidean(sim_stats[nm], obs_stats[nm], w)
}

#' @describeIn distance_functions Continuous outcome: inverse-variance weighted
#'   Euclidean distance on study-level mean differences.
spima_cont_distance <- function(sim_stats, obs_stats) {
  nm_m <- intersect(names(sim_stats$means), names(obs_stats$means))
  if (length(nm_m) == 0) return(Inf)

  # Use pre-computed precision weights if available (from cont_observed_stats),
  # otherwise fall back to pooled SD-based weights
  w <- if (!is.null(attr(obs_stats$sds, "weights"))) {
    attr(obs_stats$sds, "weights")[nm_m]
  } else {
    1 / (obs_stats$sds[nm_m]^2 + 1e-8)
  }
  dist_weighted_euclidean(sim_stats$means[nm_m], obs_stats$means[nm_m], w)
}

#' Generic (effect-size) distance: weighted Euclidean distance on effects
spima_generic_distance <- function(sim_stats, obs_stats) {
  nm <- intersect(names(sim_stats), names(obs_stats))
  if (length(nm) == 0) return(Inf)
  w <- if (!is.null(attr(obs_stats, "weights"))) {
    attr(obs_stats, "weights")[nm]
  } else {
    1 / (obs_stats[nm]^2 + 1e-8)
  }
  dist_weighted_euclidean(sim_stats[nm], obs_stats[nm], w)
}
