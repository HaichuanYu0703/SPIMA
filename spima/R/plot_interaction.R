# ---- plot_interaction.R ----
# Treatment Effect Modification Visualization for spima_int
#
# Provides S3 plot() method for spima_int objects, generating a curve
# of predicted treatment effect (risk difference or risk ratio) across
# the range of a continuous effect-modifier.

#' Extract Posterior Treatment Effects for Interaction Visualization
#'
#' Internal helper that extracts model coefficients and their uncertainty
#' from a \code{spima_int} object, then computes the predicted treatment
#' effect (absolute risk difference or risk ratio) across a range of
#' covariate values.  Uncertainty is propagated by sampling from the
#' multivariate normal approximation of the fixed effects.
#'
#' @param object A \code{spima_int} object (from \code{\link{spima_int}}).
#' @param covariate Character; covariate name to visualise.  If \code{NULL},
#'   uses the first covariate.
#' @param at Numeric vector of covariate values at which to evaluate the
#'   treatment effect.  If \code{NULL}, 50 equally-spaced points spanning
#'   the observed data range are generated.
#' @param scale \code{"absolute"} for risk difference, \code{"relative"}
#'   for risk ratio.
#' @param ci_level Confidence level (default 0.95).
#' @param n_draws Number of multivariate-normal draws (default 2000).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{at}}{Covariate values (length \code{n_at}).}
#'     \item{\code{draws}}{Matrix of posterior treatment effects
#'       (\code{n_draws} rows x \code{n_at} columns).}
#'     \item{\code{summary}}{Data frame with columns \code{estimate},
#'       \code{ci_lower}, \code{ci_upper}.}
#'     \item{\code{scale}}{\code{"absolute"} or \code{"relative"}.}
#'     \item{\code{model_type}}{\code{"pseudo-IPD"} or \code{"aggregate"}.}
#'   }
#' @keywords internal
extract_interaction_posterior <- function(object, covariate = NULL,
                                           at = NULL,
                                           scale = c("absolute", "relative"),
                                           ci_level = 0.95,
                                           n_draws = 2000) {
  scale <- match.arg(scale)

  # ---- Determine which fitted model to use ----
  if (!is.null(object$fit_pseudo) && isTRUE(object$pseudo_converged)) {
    fit <- object$fit_pseudo
    prefix <- ""           # pseudo-IPD: covariates are stored as-is
  } else if (!is.null(object$fit_aggregate) && isTRUE(object$converged)) {
    fit <- object$fit_aggregate
    prefix <- "mean_"      # aggregate: covariates are study-level means
  } else {
    stop("No converged model available for plotting. ",
         "The spima_int analysis did not produce a usable fit.")
  }

  # ---- Validate covariate ----
  cov_names <- object$input_spec$covariate
  if (is.null(covariate)) covariate <- cov_names[1]
  if (!covariate %in% cov_names)
    stop("Covariate '", covariate, "' not found in the interaction model. ",
         "Available: ", paste(cov_names, collapse = ", "))
  other_covs <- setdiff(cov_names, covariate)
  grp_col <- object$input_spec$group

  # ---- Extract fixed effects and vcov ----
  if (inherits(fit, "merMod")) {
    fe <- lme4::fixef(fit)
    vc <- as.matrix(stats::vcov(fit))
  } else {
    fe <- stats::coef(fit)
    vc <- stats::vcov(fit)
  }

  # Robust coefficient-name lookup (handles T:X vs X:T ordering)
  .find_name <- function(nm) {
    if (nm %in% names(fe)) return(nm)
    if (grepl(":", nm, fixed = TRUE)) {
      parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
      alt <- paste(rev(parts), collapse = ":")
      if (alt %in% names(fe)) return(alt)
    }
    NULL
  }

  cov_nm   <- .find_name(paste0(prefix, covariate))
  int_nm   <- .find_name(paste0(grp_col, ":", prefix, covariate))
  other_nms <- vapply(other_covs, function(cn)
    .find_name(paste0(prefix, cn)) %||% NA_character_, character(1))

  # Build the set of coefficients to include in MVN sampling
  coef_keep <- c("(Intercept)", grp_col)
  if (!is.null(cov_nm))  coef_keep <- c(coef_keep, cov_nm)
  if (!is.null(int_nm))  coef_keep <- c(coef_keep, int_nm)
  other_nms <- other_nms[!is.na(other_nms)]
  coef_keep <- c(coef_keep, other_nms)
  coef_keep <- intersect(coef_keep, names(fe))

  # Point estimates for all coefficients (including non-sampled ones)
  alpha_val <- fe["(Intercept)"]
  tau_val   <- fe[grp_col]
  gamma_val <- if (!is.null(cov_nm)) fe[cov_nm] else 0
  beta_val  <- if (!is.null(int_nm)) fe[int_nm] else 0

  other_vals <- setNames(rep(0, length(other_covs)), other_covs)
  for (cn in other_covs) {
    onm <- .find_name(paste0(prefix, cn))
    if (!is.null(onm)) other_vals[cn] <- fe[onm]
  }

  # ---- Covariate range ----
  data <- object$data
  if (prefix == "mean_") {
    raw_x <- data[[paste0("mean_", covariate)]]
    x_rng <- range(raw_x, na.rm = TRUE)
  } else {
    mcol <- paste0("mean_", covariate)
    scol <- paste0("sd_", covariate)
    if (mcol %in% names(data) && scol %in% names(data)) {
      m <- mean(data[[mcol]], na.rm = TRUE)
      s <- mean(data[[scol]], na.rm = TRUE)
      raw_x <- data[[mcol]]
      x_rng <- range(raw_x, na.rm = TRUE)
      x_rng <- c(min(x_rng[1], m - 2 * s), max(x_rng[2], m + 2 * s))
    } else {
      raw_x <- data[[covariate]]
      x_rng <- range(raw_x, na.rm = TRUE)
    }
  }

  if (is.null(at)) at <- seq(x_rng[1], x_rng[2], length.out = 50)
  n_at <- length(at)

  # Means for other covariates (held constant)
  other_means <- vapply(other_covs, function(cn) {
    mcol <- paste0("mean_", cn)
    if (mcol %in% names(data)) mean(data[[mcol]], na.rm = TRUE) else 0
  }, numeric(1))

  # ---- MVN sampling ----
  mu_samp <- fe[coef_keep]
  Sigma_samp <- vc[coef_keep, coef_keep, drop = FALSE]

  # Regularise if singular
  if (rcond(Sigma_samp) < 1e-12)
    Sigma_samp <- Sigma_samp + diag(1e-8, nrow(Sigma_samp))

  draws <- rmvn(n_draws, mu_samp, Sigma_samp)
  colnames(draws) <- coef_keep

  # Indices
  i_alpha <- match("(Intercept)", coef_keep)
  i_tau   <- match(grp_col, coef_keep)
  i_gamma <- if (!is.null(cov_nm)) match(cov_nm, coef_keep) else NA
  i_beta  <- if (!is.null(int_nm)) match(int_nm, coef_keep) else NA

  # Offset from other covariates
  other_offset <- rep(0, n_draws)
  for (cn in other_covs) {
    onm <- .find_name(paste0(prefix, cn))
    if (!is.null(onm) && onm %in% coef_keep) {
      j <- match(onm, coef_keep)
      other_offset <- other_offset + draws[, j] * other_means[cn]
    }
  }

  # ---- Compute treatment effect at each covariate value ----
  effect_mat <- matrix(NA, n_draws, n_at)
  for (j in seq_len(n_at)) {
    xj <- at[j]
    gx <- if (!is.na(i_gamma)) draws[, i_gamma] * xj else gamma_val * xj
    bx <- if (!is.na(i_beta))  draws[, i_beta] * xj  else beta_val * xj

    lp_c <- draws[, i_alpha] + other_offset + gx
    lp_t <- draws[, i_alpha] + draws[, i_tau] + other_offset + gx + bx

    p_c <- stats::plogis(lp_c)
    p_t <- stats::plogis(lp_t)

    if (scale == "absolute") {
      effect_mat[, j] <- p_t - p_c
    } else {
      effect_mat[, j] <- p_t / p_c
    }
  }
  effect_mat[!is.finite(effect_mat)] <- NA

  # ---- Summarise ----
  lo_prob <- (1 - ci_level) / 2
  hi_prob <- 1 - lo_prob

  summary_df <- data.frame(
    estimate  = colMeans(effect_mat, na.rm = TRUE),
    ci_lower  = apply(effect_mat, 2, stats::quantile, probs = lo_prob, na.rm = TRUE),
    ci_upper  = apply(effect_mat, 2, stats::quantile, probs = hi_prob, na.rm = TRUE)
  )

  list(
    at         = at,
    draws      = effect_mat,
    summary    = summary_df,
    scale      = scale,
    model_type = if (prefix == "") "pseudo-IPD" else "aggregate"
  )
}


#' Plot Treatment Effect Modification
#'
#' Generates a plot showing how the predicted treatment effect (absolute risk
#' difference or risk ratio) varies across the range of a continuous
#' covariate, based on interaction estimates from \code{\link{spima_int}}.
#'
#' The underlying model is either the pseudo-IPD individual-level GLMM
#' (preferred) or the aggregate ecological GLMM (fallback).  Uncertainty
#' is propagated by sampling from the multivariate normal approximation
#' of the fixed effects.
#'
#' @param x A \code{spima_int} object from \code{\link{spima_int}}.
#' @param covariate Character; name of the covariate to plot.  If
#'   \code{NULL}, the first covariate is used.
#' @param ci_level Confidence level for the uncertainty band (default 0.95).
#' @param at Numeric vector of covariate values at which to evaluate the
#'   treatment effect.  If \code{NULL}, 50 points are generated from the
#'   observed data range.
#' @param scale \code{"absolute"} (default) for risk difference,
#'   \code{"relative"} for risk ratio.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- spima_int(data, input_spec)
#' plot(res, covariate = "X1", scale = "absolute")
#' }
plot.spima_int <- function(x, covariate = NULL, ci_level = 0.95,
                            at = NULL, scale = c("absolute", "relative"),
                            ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.  Install it with ",
         "install.packages(\"ggplot2\").")

  scale <- match.arg(scale)

  post <- extract_interaction_posterior(
    object    = x,
    covariate = covariate,
    at        = at,
    scale     = scale,
    ci_level  = ci_level
  )

  # ---- Build the plot ----
  ref_line <- if (scale == "absolute") 0 else 1
  y_lab   <- if (scale == "absolute") "Risk Difference" else "Risk Ratio"
  title   <- paste("Treatment Effect Modification by", covariate %||% x$input_spec$covariate[1])

  plot_data <- cbind(data.frame(x = post$at), post$summary)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = "#0072B2", alpha = 0.25
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$estimate),
      colour = "#0072B2", linewidth = 1.1
    ) +
    ggplot2::geom_hline(
      yintercept = ref_line, linetype = "dashed",
      colour = "grey40", linewidth = 0.6
    ) +
    ggplot2::labs(x = covariate %||% x$input_spec$covariate[1],
                  y = y_lab, title = title) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = 13),
      axis.text    = ggplot2::element_text(size = 11),
      plot.title   = ggplot2::element_text(size = 14, face = "plain"),
      plot.margin  = ggplot2::margin(12, 14, 10, 10)
    )

  p
}
