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

#' @describeIn distance_functions Continuous outcome: sum of squared
#'   standardised differences in means and SDs.
spima_cont_distance <- function(sim_stats, obs_stats) {
  nm_m <- intersect(names(sim_stats$means), names(obs_stats$means))
  nm_s <- intersect(names(sim_stats$sds),   names(obs_stats$sds))
  if (length(nm_m) == 0) return(Inf)

  d_mean <- (sim_stats$means[nm_m] - obs_stats$means[nm_m]) /
            abs(obs_stats$means[nm_m] + 1e-8)
  d_sd <- if (length(nm_s) > 0) {
    (sim_stats$sds[nm_s] - obs_stats$sds[nm_s]) /
      abs(obs_stats$sds[nm_s] + 1e-8)
  } else numeric(0)

  all_d <- c(d_mean, d_sd)
  dist_euclidean(all_d, numeric(length(all_d)))
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
