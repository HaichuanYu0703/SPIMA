# ---- module_int.R ----
# SPI-MA Interaction Analysis Module
#
# Tests whether continuous covariate(s) modify the treatment effect,
# using only aggregate data (per-arm summary statistics).
#
# Primary method: mixed-effects logistic regression on aggregate data
#   glmer(cbind(event, n-event) ~ T * mean_X + (1 | study), family = binomial)
#
# Sensitivity: pseudo-IPD generated with aggregate estimates + individual glmer
#
# Input: data frame, one row per study-arm, with columns:
#   study, group (0/1), event, n, mean_X<k>, sd_X<k>

# =============================================================
# Main function
# =============================================================

#' SPI-MA Interaction Analysis
#'
#' Tests whether continuous covariate(s) modify the treatment effect using
#' aggregate data only.  The primary method fits a mixed-effects logistic
#' regression on the aggregate data.  A sensitivity analysis generates
#' pseudo-IPD and fits an individual-level model.
#'
#' @param data Data frame, one row per study-arm.
#' @param input_spec Named list:
#'   \describe{
#'     \item{study}{Study identifier column name.}
#'     \item{event}{Event count column name.}
#'     \item{n}{Sample size column name.}
#'     \item{group}{Treatment group column (0 = control, 1 = treatment).}
#'     \item{covariate}{Character vector of covariate name(s), e.g.
#'       \code{c("X1","X2")}.  The function looks for columns
#'       \code{mean_<name>} and \code{sd_<name>} in \code{data}.}
#'   }
#' @param rho Assumed between-covariate correlation for pseudo-IPD
#'   generation.  Default 0.
#' @param ... Additional arguments to \code{\link[lme4]{glmer}}.
#' @return \code{spima_int} object.
#' @export
spima_int <- function(data, input_spec, rho = 0, ...) {

  spima_int_validate(data, input_spec)
  cov_names <- input_spec[["covariate"]]

  mapped <- resolve_study_col(data, input_spec)
  data      <- mapped$data
  study_col <- mapped$col

  # ---- Step 1: Primary — aggregate GLMM ----
  agg_result <- .int_fit_aggregate(data, input_spec, study_col, ...)
  beta_agg  <- agg_result$beta
  gamma_agg <- agg_result$gamma

  # ---- Step 2: Sensitivity — pseudo-IPD glmer ----
  pseudo_result <- NULL
  if (agg_result$converged) {
    pseudo_result <- .int_pseudo_sensitivity(data, input_spec, study_col,
                                              agg_result$gamma_for_gen,
                                              agg_result$beta_for_gen,
                                              rho, ...)
  }

  # ---- Compose result ----
  result <- list(
    call        = match.call(),
    data        = data,
    input_spec  = input_spec,
    rho         = rho,
    K           = length(unique(data[[study_col]])),
    n_cov       = length(cov_names),

    # Aggregate model (primary)
    coefficients     = agg_result$coefficients,
    converged        = agg_result$converged,
    fit_aggregate    = agg_result$fit,
    tau_aggregate    = agg_result$tau,

    # Pseudo-IPD model (sensitivity)
    pseudo_coefficients = pseudo_result$coefficients,
    pseudo_converged    = pseudo_result$converged,
    fit_pseudo          = pseudo_result$fit,
    tau_pseudo          = pseudo_result$tau,

    # Gamma estimates
    gamma         = agg_result$gamma,
    beta_aggregate = beta_agg
  )

  class(result) <- "spima_int"
  result
}

# =============================================================
# Validation
# =============================================================

#' @rdname spima_int
#' @export
spima_int_validate <- function(data, input_spec) {
  required <- c("event", "n", "group", "covariate")
  miss <- setdiff(required, names(input_spec))
  if (length(miss)) stop("input_spec missing: ", paste(miss, collapse = ", "))

  cols <- c(input_spec$event, input_spec$n, input_spec$group)
  if (!is.null(input_spec$study)) cols <- c(cols, input_spec$study)
  for (cn in input_spec$covariate)
    cols <- c(cols, paste0("mean_", cn), paste0("sd_", cn))

  miss <- setdiff(unique(cols), names(data))
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "))

  ev <- data[[input_spec$event]]
  n  <- data[[input_spec$n]]
  if (any(ev < 0, na.rm = TRUE))  stop("Event counts must be non-negative.")
  if (any(n <= 0, na.rm = TRUE))  stop("Sample sizes must be positive.")
  if (any(ev > n, na.rm = TRUE))  stop("Events cannot exceed sample size.")

  ug <- unique(data[[input_spec$group]])
  if (!all(ug %in% c(0, 1)))
    stop("Group must be 0 (control) or 1 (treatment).")
  invisible(TRUE)
}

# =============================================================
# Step 1: Aggregate GLMM
# =============================================================

#' Fit mixed-effects logistic regression on aggregate data
#'
#' Fits glmer(cbind(event, n-event) ~ T * mean_X + (1|study))
#' on the aggregate data.  This is an ecological regression —
#' coefficients are attenuated toward zero, but the interaction
#' test is valid (conservative).
#'
#' @return list(beta, gamma, coefficients, converged, fit, tau,
#'              gamma_for_gen, beta_for_gen)
#' @keywords internal
.int_fit_aggregate <- function(data, input_spec, study_col, ...) {
  cov_names <- input_spec[["covariate"]]
  grp_col   <- input_spec[["group"]]
  ev_col    <- input_spec[["event"]]
  n_col     <- input_spec[["n"]]
  mean_cols <- paste0("mean_", cov_names)

  # Build formula: cbind(event, n-event) ~ T * mean_X1 * mean_X2 + (1|study)
  rhs_agg <- paste(c(sprintf("%s * (%s)", grp_col,
                              paste(mean_cols, collapse = " + ")),
                      sprintf("(1 | %s)", study_col)), collapse = " + ")
  form_agg <- as.formula(paste("cbind(", ev_col, ", ", n_col, " - ",
                                ev_col, ") ~", rhs_agg))

  fit <- tryCatch(
    lme4::glmer(form_agg, data = data, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa",
                                              calc.derivs = FALSE), ...),
    error = function(e) NULL)

  if (is.null(fit)) {
    # Fallback: fixed-effects glm
    rhs_fe <- paste(grp_col, "*", paste(mean_cols, collapse = " + "))
    form_fe <- as.formula(paste("cbind(", ev_col, ", ", n_col, " - ",
                                 ev_col, ") ~", rhs_fe))
    fit <- tryCatch(glm(form_fe, data = data, family = binomial),
                    error = function(e) NULL)
    if (is.null(fit))
      return(list(coefficients = NULL, converged = FALSE, fit = NULL,
                  beta = rep(NA, length(cov_names)),
                  gamma = rep(NA, length(cov_names)),
                  gamma_for_gen = rep(0, length(cov_names)),
                  beta_for_gen = rep(0, length(cov_names))))

    se <- sqrt(diag(vcov(fit)))
    coef_tab <- .int_coef_table(fit, se)
    all_fe <- fixef(fit)

    # Extract gamma (mean_X main effects) and beta (T:mean_X interactions)
    gamma <- sapply(mean_cols, function(m)
      if (m %in% names(all_fe)) all_fe[m] else 0)
    names(gamma) <- cov_names
    int_terms <- paste0(grp_col, ":", mean_cols)
    beta <- sapply(int_terms, function(t)
      if (t %in% names(all_fe)) all_fe[t] else 0)
    names(beta) <- cov_names

    return(list(coefficients = coef_tab, converged = TRUE,
                fit = fit, tau = NA,
                beta = beta, gamma = gamma,
                gamma_for_gen = gamma, beta_for_gen = beta))
  }

  # Refit if convergence issues
  if (length(fit@optinfo$conv$lme4$messages)) {
    fit2 <- tryCatch(
      lme4::glmer(form_agg, data = data, family = binomial,
                  start = lme4::getME(fit, "theta"),
                  control = lme4::glmerControl(
                    optimizer = "bobyqa", calc.derivs = FALSE,
                    optCtrl = list(maxfun = 2e5))),
      error = function(e) fit)
    if (length(fit2@optinfo$conv$lme4$messages) <
        length(fit@optinfo$conv$lme4$messages))
      fit <- fit2
  }

  se <- tryCatch(sqrt(diag(vcov(fit))), error = function(e) NULL)
  if (is.null(se))
    return(list(coefficients = NULL, converged = FALSE, fit = fit,
                beta = rep(NA, length(cov_names)),
                gamma = rep(NA, length(cov_names)),
                gamma_for_gen = rep(0, length(cov_names)),
                beta_for_gen = rep(0, length(cov_names))))

  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  tau <- if (!is.null(vc) && nrow(vc) > 0) vc$sdcor[1] else NA

  coef_tab <- .int_coef_table(fit, se)
  all_fe <- fixef(fit)

  # Extract gamma and beta
  gamma <- sapply(mean_cols, function(m)
    if (m %in% names(all_fe)) all_fe[m] else 0)
  names(gamma) <- cov_names

  int_sep <- paste0(grp_col, ":", mean_cols)
  int_comb <- paste0(mean_cols, ":", grp_col)
  beta <- numeric(length(cov_names)); names(beta) <- cov_names
  for (i in seq_along(cov_names)) {
    if (int_sep[i] %in% names(all_fe)) beta[i] <- all_fe[int_sep[i]]
    else if (int_comb[i] %in% names(all_fe)) beta[i] <- all_fe[int_comb[i]]
    else beta[i] <- 0
  }

  list(coefficients = coef_tab, converged = TRUE, fit = fit, tau = tau,
       beta = beta, gamma = gamma,
       gamma_for_gen = gamma, beta_for_gen = beta)
}

# =============================================================
# Step 2: Pseudo-IPD sensitivity analysis
# =============================================================

#' Generate pseudo-IPD and refit glmer for sensitivity analysis
#'
#' Uses gamma and beta from the aggregate model to generate pseudo-IPD,
#' then fits an individual-level glmer.
#'
#' @keywords internal
.int_pseudo_sensitivity <- function(data, input_spec, study_col,
                                     gamma, beta, rho, ...) {
  cov_names <- input_spec[["covariate"]]
  n_cov <- length(cov_names)

  # Generate X
  x_data <- .int_gen_x(data, input_spec, rho)
  if (is.null(x_data) || nrow(x_data) == 0)
    return(list(coefficients = NULL, converged = FALSE, fit = NULL))

  # Generate Y with aggregate-model estimates
  pseudo <- .int_gen_y(x_data, data, input_spec, gamma, beta)
  if (is.null(pseudo) || nrow(pseudo) == 0)
    return(list(coefficients = NULL, converged = FALSE, fit = NULL))

  # Fit individual-level glmer
  grp_col <- input_spec[["group"]]
  rhs <- paste(c(sprintf("%s * (%s)", grp_col,
                          paste(cov_names, collapse = " + ")),
                  sprintf("(1 | %s)", study_col)), collapse = " + ")
  form <- as.formula(paste("Y ~", rhs))

  fit <- tryCatch(
    lme4::glmer(form, data = pseudo, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa",
                                              calc.derivs = FALSE), ...),
    error = function(e) NULL)

  if (is.null(fit)) {
    rhs_fe <- paste(grp_col, "*", paste(cov_names, collapse = " + "))
    form_fe <- as.formula(paste("Y ~", rhs_fe))
    fit <- tryCatch(glm(form_fe, data = pseudo, family = binomial),
                    error = function(e) NULL)
    if (is.null(fit))
      return(list(coefficients = NULL, converged = FALSE, fit = NULL))
    se <- sqrt(diag(vcov(fit)))
    return(list(coefficients = .int_coef_table(fit, se), converged = TRUE,
                fit = fit, tau = NA))
  }

  if (length(fit@optinfo$conv$lme4$messages)) {
    fit2 <- tryCatch(
      lme4::glmer(form, data = pseudo, family = binomial,
                  start = lme4::getME(fit, "theta"),
                  control = lme4::glmerControl(
                    optimizer = "bobyqa", calc.derivs = FALSE,
                    optCtrl = list(maxfun = 2e5))),
      error = function(e) fit)
    fit <- fit2
  }

  se <- tryCatch(sqrt(diag(vcov(fit))), error = function(e) NULL)
  if (is.null(se))
    return(list(coefficients = NULL, converged = FALSE, fit = fit))

  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  tau <- if (!is.null(vc) && nrow(vc) > 0) vc$sdcor[1] else NA

  list(coefficients = .int_coef_table(fit, se), converged = TRUE,
       fit = fit, tau = tau)
}

# =============================================================
# X generation (shared utility)
# =============================================================

#' @keywords internal
.int_gen_x <- function(data, input_spec, rho) {
  cov_names <- input_spec[["covariate"]]
  n_cov  <- length(cov_names)
  mcols  <- paste0("mean_", cov_names)
  scol   <- paste0("sd_", cov_names)
  grp_col <- input_spec[["group"]]
  n_col   <- input_spec[["n"]]

  mapped <- resolve_study_col(data, input_spec)
  data      <- mapped$data
  study_col <- mapped$col

  build_Sigma <- function(sds) {
    sds <- pmax(sds, 0.01)
    S <- diag(sds^2, n_cov)
    if (n_cov > 1 && abs(rho) > 0) {
      for (i in 1:(n_cov - 1)) for (j in (i + 1):n_cov) {
        S[i, j] <- rho * sds[i] * sds[j]
        S[j, i] <- S[i, j]
      }
    }
    if (rcond(S) < 1e-12) S <- S + diag(1e-6, n_cov)
    S
  }

  pieces <- vector("list", nrow(data))
  for (i in seq_len(nrow(data))) {
    arm <- data[i, ]
    ni <- arm[[n_col]]
    if (ni <= 0) next

    mu <- sapply(mcols, function(cn) arm[[cn]])
    sd <- sapply(scol,  function(cn) arm[[cn]])

    Xi <- if (ni > 1) rmvn(ni, mu = mu, Sigma = build_Sigma(sd))
          else matrix(mu, 1L)
    df <- data.frame(.row_id = rep.int(i, ni),
                     study = rep.int(arm[[study_col]], ni),
                     stringsAsFactors = FALSE)
    df[[grp_col]] <- rep.int(arm[[grp_col]], ni)
    for (j in seq_len(n_cov)) df[[cov_names[j]]] <- Xi[, j]
    pieces[[i]] <- df
  }
  do.call(rbind, pieces)
}

# =============================================================
# Y generation for pseudo-IPD
# =============================================================

#' @keywords internal
.int_gen_y <- function(x_data, agg_data, input_spec, gamma, beta) {
  cov_names <- names(gamma)
  grp_col  <- input_spec[["group"]]
  ev_col   <- input_spec[["event"]]
  n_col    <- input_spec[["n"]]

  mapped <- resolve_study_col(agg_data, input_spec)
  agg      <- mapped$data
  study_col <- mapped$col
  studies   <- unique(x_data$study)

  pieces <- vector("list", length(studies))
  ii <- 1L

  for (sid in studies) {
    ca <- agg[agg[[study_col]] == sid & agg[[grp_col]] == 0, , drop = FALSE]
    if (nrow(ca) > 0) {
      p_c <- ca[[ev_col]][1] / ca[[n_col]][1]
      mxc <- sapply(paste0("mean_", cov_names), function(cn) ca[[cn]][1])
      alpha_j <- qlogis(p_c + 1e-8) - sum(gamma * mxc)
    } else {
      alpha_j <- 0
    }

    idx <- which(x_data$study == sid)
    sdat <- x_data[idx, , drop = FALSE]

    # Control arm
    ci <- which(sdat[[grp_col]] == 0)
    if (length(ci)) {
      Xc <- as.matrix(sdat[ci, cov_names, drop = FALSE])
      sdat$Y <- NA_real_
      sdat$Y[ci] <- rbinom(length(ci), 1, plogis(c(alpha_j + Xc %*% gamma)))
    }

    # Treatment arm
    ti <- which(sdat[[grp_col]] == 1)
    if (length(ti)) {
      Xt <- as.matrix(sdat[ti, cov_names, drop = FALSE])
      sdat$Y[ti] <- rbinom(length(ti), 1, plogis(c(alpha_j + Xt %*% (gamma + beta))))
    }

    sdat$Y[is.na(sdat$Y)] <- 0L
    pieces[[ii]] <- sdat; ii <- ii + 1L
  }

  do.call(rbind, pieces)
}

# =============================================================
# Coefficient table helper
# =============================================================

#' @keywords internal
.int_coef_table <- function(fit, se) {
  fe <- fixef(fit)
  tr <- names(fe)
  int <- grep(":", tr, value = TRUE)
  tab <- data.frame(term = tr, estimate = unname(fe), se = unname(se[tr]),
                    stringsAsFactors = FALSE)
  tab$ci_l    <- tab$estimate - 1.96 * tab$se
  tab$ci_u    <- tab$estimate + 1.96 * tab$se
  tab$z_value <- tab$estimate / tab$se
  tab$p_value <- 2 * pnorm(-abs(tab$z_value))
  tab$is_interaction <- tab$term %in% int
  tab
}

# =============================================================
# S3 methods
# =============================================================

#' @export
print.spima_int <- function(x, ...) {
  cat("SPI-MA Interaction Analysis\n")
  cat("  Studies:", x$K, "\n")
  cat("  Covariates:", paste(x$input_spec$covariate, collapse = ", "), "\n")
  cat("  Assumed correlation ?:", x$rho, "\n")
  if (!is.null(x$tau_aggregate)) cat("  ? (heterogeneity SD):", round(x$tau_aggregate, 4), "\n")

  cat("\n--- Aggregate model (primary) ---\n")
  if (!is.null(x$coefficients) && nrow(x$coefficients)) {
    int <- x$coefficients[x$coefficients$is_interaction, ]
    for (i in seq_len(nrow(int))) {
      r <- int[i, ]
      cat(sprintf("  %s:\n", r$term))
      cat(sprintf("    Estimate: %.4f   SE: %.4f\n", r$estimate, r$se))
      cat(sprintf("    95%% CI: [%.4f, %.4f]   p = %.4f\n", r$ci_l, r$ci_u, r$p_value))
    }
  } else {
    cat("  (model did not converge)\n")
  }

  if (!is.null(x$pseudo_coefficients) && nrow(x$pseudo_coefficients)) {
    cat("\n--- Pseudo-IPD model (sensitivity) ---\n")
    int <- x$pseudo_coefficients[x$pseudo_coefficients$is_interaction, ]
    for (i in seq_len(nrow(int))) {
      r <- int[i, ]
      cat(sprintf("  %s:\n", r$term))
      cat(sprintf("    Estimate: %.4f   SE: %.4f\n", r$estimate, r$se))
      cat(sprintf("    95%% CI: [%.4f, %.4f]   p = %.4f\n", r$ci_l, r$ci_u, r$p_value))
    }
  }

  invisible(x)
}

#' @export
summary.spima_int <- function(object, ...) {
  structure(object, class = c("summary.spima_int", "spima_int"))
}
#' @export
print.summary.spima_int <- function(x, ...) NextMethod()

#' @export
as.data.frame.spima_int <- function(x, ...) {
  if (is.null(x$coefficients)) return(data.frame())
  x$coefficients[, !names(x$coefficients) %in% "is_interaction", drop = FALSE]
}
