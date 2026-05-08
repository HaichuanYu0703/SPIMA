#' spima: Simulated Pseudo-Individual Data Meta-Analysis
#'
#' The main entry point. Dispatches to the appropriate module based on
#' \code{outcome_type} and runs ABC-SMC for meta-analytic inference.
#'
#' @param data A data frame of study-level summary statistics; one row per
#'   study (or per study-arm when \code{group} is specified in
#'   \code{input_spec}).
#' @param outcome_type Outcome type: \code{"binary"}, \code{"continuous"},
#'   or \code{"generic"}.
#' @param input_spec A named list mapping column names to roles. The
#'   required entries depend on \code{outcome_type}:
#'   \describe{
#'     \item{binary}{\code{event}, \code{n}, optionally \code{group} and
#'       \code{study}.}
#'     \item{continuous}{\code{mean}, \code{sd}, \code{n}, optionally
#'       \code{group} and \code{study}.  Alternative formats:
#'       \code{median}, \code{q1}, \code{q3}, \code{n} (median+IQR);
#'       \code{median}, \code{min}, \code{max}, \code{n} (range);
#'       or \code{median}, \code{q1}, \code{q3}, \code{min}, \code{max},
#'       \code{n} (five-number summary).}
#'     \item{generic}{\code{yi} (effect size), \code{sei} (standard error),
#'       optionally \code{study}.}
#'   }
#' @param prior A \code{spima_prior} object created by \code{prior()}.
#' @param smc_control An \code{smc_control} list created by
#'   \code{smc_control()}.
#' @param parallel Logical; if \code{TRUE}, use \code{parallel::mclapply}
#'   for particle simulations (Unix only).
#' @param subgroup Optional column name for subgroup analysis.  When
#'   specified, the analysis is run separately for each level of this
#'   variable.
#' @param family Distributional family for the pseudo-IPD likelihood.
#'   Only used when \code{outcome_type = "continuous"}.
#'   \code{"gaussian"} (default) assumes normally-distributed outcomes
#'   and estimates Mean Differences. \code{"Gamma"} assumes
#'   Gamma-distributed outcomes and estimates Rate Ratios
#'   on the log scale (log-RR).
#' @param ... Additional arguments passed to module functions.
#' @return A \code{spima} object with components:
#'   \item{call}{The matched call.}
#'   \item{outcome_type}{The outcome type.}
#'   \item{abc_result}{Full ABC-SMC output (generations, posterior, etc.).}
#'   \item{data}{The input data.}
#'   \item{input_spec}{The column mapping.}
#' @export
#'
#' @examples
#' \dontrun{
#' # Binary outcome meta-analysis (two-arm per study)
#' data_bin <- data.frame(
#'   study = 1:4,
#'   group = c(0, 1, 0, 1, 0, 1, 0, 1),
#'   event = c(30, 45, 28, 32, 40, 58, 18, 22),
#'   n     = c(100, 100, 80, 80, 120, 120, 60, 60)
#' )
#' res <- spima(data_bin, "binary",
#'              input_spec = list(study = "study", event = "event",
#'                                n = "n", group = "group"),
#'              prior = prior(mu = "normal(0, 10)", tau = "halfnormal(0, 1)"),
#'              smc_control = smc_control(n_particles = 500, n_generations = 5))
#' }
spima <- function(data, outcome_type = c("binary", "continuous", "generic"),
                  input_spec, prior, smc_control, parallel = FALSE,
                  subgroup = NULL, family = c("gaussian", "Gamma"), ...) {

  outcome_type <- match.arg(outcome_type)
  family <- match.arg(family)

  # Resolve study column default
  if (is.null(input_spec[["study"]])) {
    if ("study" %in% names(data)) {
      input_spec[["study"]] <- "study"
    }
  }

  # ---- Subgroup analysis ----
  if (!is.null(subgroup)) {
    if (!subgroup %in% names(data)) {
      stop("Subgroup column '", subgroup, "' not found in data.")
    }
    levels <- unique(data[[subgroup]])
    if (length(levels) < 2) {
      stop("Subgroup variable must have at least 2 levels.")
    }

    results <- list()
    for (lv in levels) {
      if (smc_control$verbose) {
        cat("--- Subgroup: ", subgroup, " = ", lv, " ---\n", sep = "")
      }
      sub_data <- data[data[[subgroup]] == lv, , drop = FALSE]

      sub_result <- switch(outcome_type,
        binary       = run_module_bin(sub_data, input_spec, prior, smc_control, ...),
        continuous   = run_module_cont(sub_data, input_spec, prior, smc_control,
                                        family = family, ...),
        generic      = run_module_generic(sub_data, input_spec, prior, smc_control, ...)
      )

      if (!is.null(sub_result) && length(sub_result$generations) > 0) {
        results[[as.character(lv)]] <- sub_result
      } else if (smc_control$verbose) {
        cat("  (no results)\n")
      }
    }

    if (length(results) < 1) {
      stop("No subgroups produced valid results.")
    }

    return(structure(
      list(
        call         = match.call(),
        outcome_type = outcome_type,
        family       = family,
        results      = results,
        data         = data,
        input_spec   = input_spec,
        subgroup_col = subgroup
      ),
      class = "spima_subgroup"
    ))
  }

  # ---- Standard (non-subgroup) dispatch ----
  result <- switch(outcome_type,
    binary      = run_module_bin(data, input_spec, prior, smc_control, ...),
    continuous  = run_module_cont(data, input_spec, prior, smc_control,
                                   family = family, ...),
    generic     = run_module_generic(data, input_spec, prior, smc_control, ...)
  )

  structure(
    list(
      call         = match.call(),
      outcome_type = outcome_type,
      family       = family,
      abc_result   = result,
      data         = data,
      input_spec   = input_spec
    ),
    class = "spima"
  )
}

# ---------- Module dispatchers ----------

run_module_bin <- function(data, input_spec, prior_obj, ctrl, ...) {
  # Validate
  spima_bin_validate(data, input_spec)

  # Compute observed summary stats once
  obs_stats <- spima_bin_observed_stats(data, input_spec)

  # Build simulation wrapper
  sim_fn <- function(theta, ...) {
    ipd <- spima_bin_simulate(data, theta, input_spec)
    if (is.null(ipd) || nrow(ipd) == 0) return(NULL)
    res <- spima_bin_analyze(ipd, input_spec)
    if (is.null(res$summary_stats)) return(NULL)
    res$summary_stats
  }

  run_abc_smc(prior_obj, sim_fn, spima_bin_distance, obs_stats, ctrl)
}

run_module_cont <- function(data, input_spec, prior_obj, ctrl, ...,
                             family = "gaussian") {

  if (family == "Gamma") {
    # ---- Gamma module ----
    spima_gamma_validate(data, input_spec)
    obs_stats <- spima_gamma_observed_stats(data, input_spec)

    sim_fn <- function(theta, ...) {
      ipd <- spima_gamma_simulate(data, theta, input_spec)
      if (is.null(ipd) || nrow(ipd) == 0) return(NULL)
      res <- spima_gamma_analyze(ipd, input_spec, quick = TRUE)
      if (is.null(res$summary_stats)) return(NULL)
      res$summary_stats
    }

    run_abc_smc(prior_obj, sim_fn, spima_gamma_distance, obs_stats, ctrl)

  } else {
    # ---- Gaussian module ----
    has_mean_sd <- all(c("mean", "sd") %in% names(input_spec))
    has_mean_range <- all(c("mean", "min", "max") %in% names(input_spec)) && !"sd" %in% names(input_spec)

    if (has_mean_sd || has_mean_range) {
      # ---- Standard mean/SD path ----
      spima_cont_validate(data, input_spec)
      obs_stats <- spima_cont_observed_stats(data, input_spec)

      sim_fn <- function(theta, ...) {
        ipd <- spima_cont_simulate(data, theta, input_spec)
        if (is.null(ipd) || nrow(ipd) == 0) return(NULL)
        res <- spima_cont_analyze(ipd, input_spec)
        if (is.null(res$summary_stats)) return(NULL)
        res$summary_stats
      }

      run_abc_smc(prior_obj, sim_fn, spima_cont_distance, obs_stats, ctrl)
    } else {
      # ---- Quantile-format path (median/IQR, range, five-number) ----
      spima_cont_validate(data, input_spec)
      obs_stats <- spima_cont_quantile_observed_stats(data, input_spec)

      sim_fn <- function(theta, ...) {
        ipd <- spima_cont_simulate(data, theta, input_spec)
        if (is.null(ipd) || nrow(ipd) == 0) return(NULL)
        res <- spima_cont_quantile_analyze(ipd, input_spec)
        if (is.null(res$summary_stats)) return(NULL)
        res$summary_stats
      }

      run_abc_smc(prior_obj, sim_fn, spima_generic_distance, obs_stats, ctrl)
    }
  }
}

run_module_generic <- function(data, input_spec, prior_obj, ctrl, ...) {
  spima_generic_validate(data, input_spec)
  obs_stats <- spima_generic_observed_stats(data, input_spec)

  sim_fn <- function(theta, ...) {
    yi_sim <- spima_generic_simulate(data, theta, input_spec)
    res <- spima_generic_analyze(yi_sim, input_spec)
    res$summary_stats
  }

  run_abc_smc(prior_obj, sim_fn, spima_generic_distance, obs_stats, ctrl)
}

