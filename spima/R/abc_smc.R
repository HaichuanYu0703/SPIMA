#' Control Parameters for ABC-SMC
#'
#' @param n_particles Number of particles (simulations) per generation.
#' @param n_generations Maximum number of SMC generations.
#' @param epsilon_init Initial acceptance threshold. If \code{NULL}, it is
#'   set to the median of distances from an initial pilot run.
#' @param epsilon_decay Multiplicative factor applied to epsilon each
#'   generation (0 < decay < 1).
#' @param ess_min Minimum effective-sample-size ratio (relative to
#'   \code{n_particles}); algorithm halts when ESS drops below this.
#' @param kernel Perturbation kernel type: \code{"gaussian"} (default).
#' @param accept_rate_target Target acceptance rate used for adaptive
#'   epsilon tuning.
#' @param verbose Print progress information?
#' @param parallel Logical; if \code{TRUE}, run particle simulations in
#'   parallel using \code{parallel::mclapply} (Unix) or a PSOCK cluster
#'   (Windows).
#' @param n_cores Number of CPU cores for parallel execution.  If
#'   \code{NULL} (default), uses \code{getOption("mc.cores")} or
#'   \code{parallel::detectCores() - 1}.
#' @return A list of class \code{smc_control}.
#' @export
#'
#' @examples
#' smc_control(n_particles = 500, n_generations = 8)
smc_control <- function(n_particles      = 2000,
                        n_particles_max  = 10000,
                        n_generations    = 10,
                        epsilon_init     = NULL,
                        epsilon_decay    = 0.85,
                        ess_min          = 0.3,
                        kernel           = "gaussian",
                        accept_rate_target = 0.2,
                        verbose          = TRUE,
                        parallel         = FALSE,
                        n_cores          = NULL) {
  stopifnot(n_particles > 0, n_particles_max >= n_particles)
  stopifnot(n_generations > 0)
  stopifnot(epsilon_decay > 0, epsilon_decay < 1)
  stopifnot(ess_min > 0, ess_min <= 1)
  stopifnot(accept_rate_target > 0, accept_rate_target < 1)

  # Determine default cores
  if (parallel && is.null(n_cores)) {
    n_cores <- getOption("mc.cores", parallel::detectCores())
    n_cores <- max(1, n_cores - 1)  # leave one free
  }

  structure(
    list(
      n_particles         = n_particles,
      n_particles_max     = n_particles_max,
      n_generations       = n_generations,
      epsilon_init        = epsilon_init,
      epsilon_decay       = epsilon_decay,
      ess_min             = ess_min,
      kernel              = kernel,
      accept_rate_target  = accept_rate_target,
      verbose             = verbose,
      parallel            = parallel,
      n_cores             = n_cores
    ),
    class = "smc_control"
  )
}

#' Run ABC-SMC Inference
#'
#' @param prior_obj A \code{spima_prior} object.
#' @param sim_fn Simulation function: \code{function(theta, ...)} returning
#'   simulated summary statistics.
#' @param distance_fn Distance function: \code{function(sim, obs)} returning
#'   a scalar.
#' @param obs_stats Observed (target) summary statistics.
#' @param ctrl An \code{smc_control} list.
#' @param ... Additional arguments passed to \code{sim_fn}.
#' @return A list of class \code{spima_abc} containing posterior samples,
#'   weights, diagnostics, and generation records.
run_abc_smc <- function(prior_obj, sim_fn, distance_fn, obs_stats, ctrl, ...) {

  P <- ctrl$n_particles
  n_par <- length(prior_obj)
  par_names <- names(prior_obj)

  # --- Setup parallel cluster (Windows-compatible) ---
  if (ctrl$parallel && .Platform$OS.type == "windows") {
    cl <- tryCatch({
      cl_tmp <- parallel::makeCluster(ctrl$n_cores)
      invisible(parallel::clusterEvalQ(cl_tmp, {
        require(spima, quietly = TRUE)
        require(lme4, quietly = TRUE)
        NULL
      }))
      cl_tmp
    }, error = function(e) {
      if (ctrl$verbose) cat("  Parallel cluster failed (", e$message,
                            "); reverting to sequential.\n")
      NULL
    })
    if (!is.null(cl)) {
      on.exit(parallel::stopCluster(cl))
      ctrl$._cl <- cl
    } else {
      ctrl$parallel <- FALSE
    }
  }

  # --- Generation 1: prior sampling ---
  if (ctrl$verbose) cat("ABC-SMC generation 1 /", ctrl$n_generations, "...\n")

  theta <- sample_prior(P, prior_obj)
  dists <- numeric(P)
  accept <- rep(TRUE, P)

  # If epsilon_init is NULL, pilot run to determine initial epsilon
  if (is.null(ctrl$epsilon_init)) {
    if (ctrl$verbose) cat("  Pilot run to determine initial epsilon...\n")
    d_pilot <- .run_simulations(theta, sim_fn, distance_fn, obs_stats, ctrl, ...)
    # Remove infinite / NA distances
    d_ok <- d_pilot[is.finite(d_pilot)]
    if (length(d_ok) < 10) {
      stop("Pilot run: fewer than 10 particles with finite distance. ",
           "Check priors or simulation function.")
    }
    # Set epsilon as ~30th percentile of pilot distances (roughly 30% accept)
    eps <- quantile(d_ok, probs = 0.3, na.rm = TRUE)
    if (ctrl$verbose) cat("  Initial epsilon set to", signif(eps, 3), "\n")
  } else {
    eps <- ctrl$epsilon_init
  }

  # Generation 1 with initial epsilon
  dists <- .run_simulations(theta, sim_fn, distance_fn, obs_stats, ctrl, ...)

  weights <- rep(1 / P, P)
  idx_ok <- is.finite(dists) & (dists <= eps)
  n_ok <- sum(idx_ok)

  if (n_ok < 5) {
    # Relax epsilon if too few accepted
    eps <- quantile(dists[is.finite(dists)], probs = 0.5, na.rm = TRUE)
    idx_ok <- is.finite(dists) & (dists <= eps)
    n_ok <- sum(idx_ok)
    if (ctrl$verbose) cat("  Relaxed epsilon to", signif(eps, 3),
                          "(", n_ok, "accepted)\n")
  }

  if (n_ok == 0) stop("No particles accepted in generation 1.")

  # Reweight
  weights[!idx_ok] <- 0
  weights <- weights / sum(weights)

  generations <- list()
  generations[[1]] <- list(
    theta   = theta,
    dists   = dists,
    weights = weights,
    epsilon = eps,
    n_accept = n_ok,
    ess     = effective_sample_size(weights)
  )

  # --- Subsequent generations ---
  for (gen in seq_len(ctrl$n_generations - 1) + 1) {

    if (ctrl$verbose) {
      cat("ABC-SMC generation", gen, "/", ctrl$n_generations, "...  ")
    }

    prev <- generations[[gen - 1]]

    # ESS check
    ess <- effective_sample_size(prev$weights)
    if (ess < ctrl$ess_min * P) {
      if (ctrl$verbose) cat("ESS =", round(ess, 1),
                            "<", ctrl$ess_min * P, " -> stopping.\n")
      break
    }

    # Update epsilon
    prev_eps <- prev$epsilon
    new_eps <- prev_eps * ctrl$epsilon_decay

    # Covariance for perturbation kernel (from weighted previous gen)
    idx_active <- prev$weights > 0
    if (sum(idx_active) < 2 * n_par) break
    theta_prev <- prev$theta[idx_active, , drop = FALSE]
    w_prev     <- prev$weights[idx_active]

    # Weighted covariance, scaled by 2 (standard ABC-SMC practice)
    Sigma <- weighted_cov(theta_prev, w_prev) * 2
    # Regularise if singular
    if (any(is.na(Sigma)) || rcond(Sigma) < 1e-12) {
      Sigma <- diag(diag(Sigma) + 1e-6, n_par)
      diag(Sigma) <- pmax(diag(Sigma), 1e-6)
    }

    # Sample particles from previous generation
    new_theta <- matrix(NA, P, n_par)
    colnames(new_theta) <- par_names
    new_dists   <- rep(Inf, P)
    new_weights <- rep(0, P)
    n_accept <- 0

    # Step 1: generate all candidate particles (resample + perturb)
    candidates <- matrix(NA, P, n_par)
    colnames(candidates) <- par_names
    for (i in seq_len(P)) {
      idx <- sample(seq_len(P), size = 1, prob = prev$weights)
      candidate <- prev$theta[idx, ]
      candidates[i, ] <- perturb_particle(candidate, Sigma, prior_obj, ctrl)
    }

    # Step 2: simulate all candidates (parallel if enabled)
    all_dists <- .run_simulations(candidates, sim_fn, distance_fn,
                                  obs_stats, ctrl, ...)

    # Step 3: process results (accept/reject + importance weights)
    for (i in seq_len(P)) {
      perturbed <- candidates[i, ]
      d <- all_dists[i]

      if (any(is.na(perturbed)) || !is.finite(d) || d > new_eps) next

      # Accept
      new_theta[i, ]  <- perturbed
      new_dists[i]    <- d
      n_accept        <- n_accept + 1

      # Compute importance weight
      prior_dens <- exp(log_prior_density(perturbed, prior_obj))

      kern_sum <- 0
      for (j in seq_len(P)) {
        if (prev$weights[j] > 0) {
          kern_sum <- kern_sum + prev$weights[j] * kernel_density(
            perturbed, prev$theta[j, ], Sigma, ctrl$kernel
          )
        }
      }
      new_weights[i] <- if (kern_sum > 0) prior_dens / kern_sum else 0
    }

    if (n_accept < 5) {
      # Hard threshold: if too few accepted, keep previous epsilon and retry
      if (ctrl$verbose) cat("Only", n_accept, "accepted -> stopping.\n")
      break
    }

    # Adaptive particle count: increase N if acceptance rate is too low
    accept_rate <- n_accept / P
    if (accept_rate < 0.5 * ctrl$accept_rate_target && P < ctrl$n_particles_max) {
      P_new <- min(P * 2, ctrl$n_particles_max)
      if (ctrl$verbose) cat("  Low acceptance (", round(accept_rate, 3),
                            "); increasing particles", P, "->", P_new, "\n")
      P <- P_new
    }

    # Normalise weights
    new_weights <- new_weights / sum(new_weights)
    new_ess <- effective_sample_size(new_weights)

    generations[[gen]] <- list(
      theta    = new_theta,
      dists    = new_dists,
      weights  = new_weights,
      epsilon  = new_eps,
      n_accept = n_accept,
      ess      = new_ess
    )

    if (ctrl$verbose) {
      cat("eps =", signif(new_eps, 3), "  accepted =", n_accept,
          "  ESS =", round(new_ess, 1), "\n")
    }

    # Early stop if ESS collapses
    if (new_ess < ctrl$ess_min * P) {
      if (ctrl$verbose) cat("ESS below threshold -> stopping.\n")
      break
    }
  }

  # Collate final posterior (last generation, or best of all)
  last <- generations[[length(generations)]]
  idx_final <- last$weights > 0

  posterior <- list(
    theta      = last$theta[idx_final, , drop = FALSE],
    weights    = last$weights[idx_final],
    distances  = last$dists[idx_final],
    epsilon    = last$epsilon,
    ess        = last$ess
  )

  # Posterior summary -- handle arbitrary parameter names
  pnames <- colnames(posterior$theta)
  summary_list <- list()
  for (pn in pnames) {
    summary_list[[pn]] <- weighted_summary(posterior$theta[, pn, drop = TRUE],
                                           posterior$weights)
  }

  result <- list(
    posterior   = posterior,
    generations = generations,
    summary     = summary_list,
    control     = ctrl,
    call        = match.call()
  )
  class(result) <- c("spima_abc", "list")
  result
}

# ---------- Parallel dispatch ----------

#' Worker function for PSOCK cluster (package-level, avoids .GlobalEnv capture)
#'
#' @param i Particle index (row of \code{theta}).
#' @param theta Matrix of parameter vectors, one row per particle.
#' @param sim_fn Simulation function \code{function(theta, ...)}.
#' @param distance_fn Distance function \code{function(sim, obs)}.
#' @param obs_stats Observed summary statistics.
#' @param extra_args List of additional arguments passed to \code{sim_fn}.
#' @return A scalar distance, or \code{Inf} if simulation fails.
#' @keywords internal
.run_psock_worker <- function(i, theta, sim_fn, distance_fn, obs_stats, extra_args) {
  if (any(is.na(theta[i, ]))) return(Inf)
  sim <- do.call(sim_fn, c(list(theta[i, , drop = TRUE]), extra_args))
  if (is.null(sim)) return(Inf)
  distance_fn(sim, obs_stats)
}

#' Run multiple simulations, optionally in parallel
#'
#' @param theta Matrix of parameter vectors (one row per particle).
#' @param sim_fn Simulation function \code{function(theta, ...)}.
#' @param distance_fn Distance function \code{function(sim, obs)}.
#' @param obs_stats Observed summary statistics.
#' @param ctrl An \code{smc_control} list (controls parallel dispatch).
#' @param ... Additional arguments passed to \code{sim_fn}.
#' @return Numeric vector of distances, one per particle.
#' @keywords internal
.run_simulations <- function(theta, sim_fn, distance_fn, obs_stats, ctrl, ...) {
  P <- nrow(theta)
  if (ctrl$parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (!is.null(ctrl$._cl)) {
      # Windows PSOCK cluster parallel
      extra_args <- list(...)
      parallel::clusterSetRNGStream(ctrl$._cl)
      results <- parallel::parLapply(ctrl$._cl, seq_len(P),
                                      .run_psock_worker,
                                      theta, sim_fn, distance_fn,
                                      obs_stats, extra_args)
      unlist(results)
    } else {
      # Unix forking parallel
      results <- parallel::mclapply(seq_len(P), function(i) {
        if (any(is.na(theta[i, ]))) return(Inf)
        sim <- sim_fn(theta[i, ], ...)
        if (is.null(sim)) return(Inf)
        distance_fn(sim, obs_stats)
      }, mc.cores = ctrl$n_cores, mc.set.seed = TRUE)
      unlist(results)
    }
  } else {
    dists <- numeric(P)
    for (i in seq_len(P)) {
      if (any(is.na(theta[i, ]))) {
        dists[i] <- Inf
        next
      }
      sim <- sim_fn(theta[i, ], ...)
      dists[i] <- if (is.null(sim)) Inf else distance_fn(sim, obs_stats)
    }
    dists
  }
}

# ---------- helpers ----------

effective_sample_size <- function(w) {
  w <- w[w > 0]
  if (length(w) < 2) return(1)
  1 / sum(w^2)
}

weighted_cov <- function(x, w) {
  if (nrow(x) < 2) return(diag(ncol(x)))
  w <- w / sum(w)
  center <- colSums(w * x)
  xc <- sweep(x, 2, center)
  t(xc) %*% diag(w) %*% xc
}

weighted_summary <- function(x, w) {
  w <- w / sum(w)
  m <- sum(w * x)
  v <- sum(w * (x - m)^2)
  c(mean = m, sd = sqrt(v),
    q2.5 = weighted_quantile(x, w, 0.025),
    q25  = weighted_quantile(x, w, 0.25),
    q50  = weighted_quantile(x, w, 0.50),
    q75  = weighted_quantile(x, w, 0.75),
    q97.5 = weighted_quantile(x, w, 0.975))
}

weighted_quantile <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w)
  cs <- cumsum(w)
  idx <- which.max(cs >= p)
  if (idx == 1) x[1] else x[idx]
}

perturb_particle <- function(theta, Sigma, prior_obj, ctrl) {
  if (ctrl$kernel == "gaussian") {
    theta_new <- as.vector(rmvn(1, mu = theta, Sigma = Sigma))
  } else {
    theta_new <- theta + rnorm(length(theta), 0, sqrt(diag(Sigma)))
  }
  names(theta_new) <- names(theta)

  # Reject if outside prior support (for bounded priors)
  for (nm in names(theta_new)) {
    if (nm %in% names(prior_obj)) {
      d <- prior_obj[[nm]]$dfun(theta_new[nm])
      if (is.na(d) || d <= 0) return(rep(NA, length(theta_new)))
    }
  }

  theta_new
}

kernel_density <- function(x, mean, Sigma, kernel_type) {
  if (kernel_type == "gaussian") {
    k <- length(x)
    det_sig <- det(Sigma)
    if (det_sig <= 0) return(0)
    diff <- x - mean
    exp_val <- try(-0.5 * t(diff) %*% solve(Sigma, diff), silent = TRUE)
    if (inherits(exp_val, "try-error")) return(0)
    (2 * pi)^(-k / 2) * det_sig^(-0.5) * exp(as.numeric(exp_val))
  } else {
    # fallback: product of independent normals
    prod(dnorm(x, mean = mean, sd = sqrt(diag(Sigma)) + 1e-8))
  }
}
