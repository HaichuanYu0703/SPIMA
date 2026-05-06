# -------- Helper operators --------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Parameter label helper: maps "mu" to a family-specific display name
# Only meaningful for continuous outcome type.
.param_label <- function(nm, family, outcome_type = "continuous") {
  if (outcome_type != "continuous" || is.null(family)) return(nm)
  if (nm == "mu") {
    switch(family,
      Gamma    = "Rate Ratio (log scale)",
      gaussian = "Mean Difference",
      nm)
  } else {
    nm
  }
}

# -------- Study column resolution --------

#' Resolve the study column: if input_spec has "study", use it;
#' else if data has a "study" column, use it; else assign row IDs.
#'
#' @param data A data frame with study-level summary statistics.
#' @param input_spec A named list mapping column names to roles; may
#'   contain \code{"study"}.
#' @return A list with components \code{data}, \code{col} (resolved
#'   column name), and \code{ids} (unique study identifiers).
#' @keywords internal
resolve_study_col <- function(data, input_spec) {
  if (!is.null(input_spec[["study"]]) &&
      input_spec[["study"]] %in% names(data)) {
    col <- input_spec[["study"]]
  } else if ("study" %in% names(data)) {
    col <- "study"
  } else {
    # No study column: assign integer IDs per unique row pattern
    data$.spima_id <- seq_len(nrow(data))
    col <- ".spima_id"
  }
  ids <- unique(data[[col]])
  list(data = data, col = col, ids = ids)
}

# -------- S3 methods --------

#' @export
print.spima <- function(x, ...) {
  cat("spima meta-analysis (", x$outcome_type, ")", sep = "")
  if (!is.null(x$family) && x$outcome_type == "continuous")
    cat(" [family: ", x$family, "]", sep = "")
  cat("\n")
  cat("  ABC-SMC generations:", length(x$abc_result$generations), "\n")

  # Show DTA derived quantities first
  if (x$outcome_type == "dta" && !is.null(x$abc_result$dta_derived)) {
    dd <- x$abc_result$dta_derived
    cat("\n  Derived quantities:\n")
    cat("    Sensitivity:  ",
        "mean =", signif(dd$sens["mean"], 3),
        "  95% CI [", signif(dd$sens["q2.5"], 3), ",",
                     signif(dd$sens["q97.5"], 3), "]\n")
    cat("    Specificity:  ",
        "mean =", signif(dd$spec["mean"], 3),
        "  95% CI [", signif(dd$spec["q2.5"], 3), ",",
                     signif(dd$spec["q97.5"], 3), "]\n")
    cat("    DOR:          ",
        "mean =", signif(dd$dor["mean"], 3),
        "  95% CI [", signif(dd$dor["q2.5"], 3), ",",
                     signif(dd$dor["q97.5"], 3), "]\n")
    cat("    Positive LR:  ",
        "mean =", signif(dd$plr["mean"], 3),
        "  95% CI [", signif(dd$plr["q2.5"], 3), ",",
                     signif(dd$plr["q97.5"], 3), "]\n")
    cat("    Negative LR:  ",
        "mean =", signif(dd$nlr["mean"], 3),
        "  95% CI [", signif(dd$nlr["q2.5"], 3), ",",
                     signif(dd$nlr["q97.5"], 3), "]\n")
    cat("\n  Model parameters (logit scale):\n")
  }

  post <- x$abc_result$summary
  for (nm in names(post)) {
    label <- .param_label(nm, x$family, x$outcome_type)
    cat("  ", label, ":\n", sep = "")
    cat("    mean =", signif(post[[nm]]["mean"], 4),
        "  sd =", signif(post[[nm]]["sd"], 4),
        "  95% CI [", signif(post[[nm]]["q2.5"], 4), ",",
                     signif(post[[nm]]["q97.5"], 4), "]\n")
  }
  invisible(x)
}

#' @export
summary.spima <- function(object, ...) {
  post <- object$abc_result$summary
  gen <- object$abc_result$generations
  last <- gen[[length(gen)]]

  out <- list(
    outcome_type  = object$outcome_type,
    family        = object$family %||% "gaussian",
    n_generations = length(gen),
    n_particles   = nrow(last$theta),
    ess           = last$ess,
    final_epsilon = last$epsilon
  )
  # Include all parameter summaries
  for (nm in names(post)) {
    out[[nm]] <- post[[nm]]
  }
  class(out) <- "summary.spima"
  out
}

#' @export
print.summary.spima <- function(x, ...) {
  cat("Outcome type:", x$outcome_type, "\n")
  if (!is.null(x$family) && x$outcome_type == "continuous") cat("Family: ", x$family, "\n", sep = "")
  cat("Generations: ", x$n_generations, "\n")
  cat("Final particles: ", x$n_particles, "\n")
  cat("Final ESS: ", round(x$ess, 1), "\n")
  cat("Final epsilon: ", signif(x$final_epsilon, 4), "\n\n")

  param_names <- setdiff(names(x), c("outcome_type", "family", "n_generations",
                                      "n_particles", "ess", "final_epsilon"))
  for (nm in param_names) {
    label <- .param_label(nm, x$family, x$outcome_type)
    cat(label, ":\n", sep = "")
    cat("  Mean: ", signif(x[[nm]]["mean"], 4), "\n")
    cat("  SD:   ", signif(x[[nm]]["sd"], 4), "\n")
    cat("  95% CI: [", signif(x[[nm]]["q2.5"], 4), ", ",
                      signif(x[[nm]]["q97.5"], 4), "]\n\n")
  }
  invisible(x)
}

# internal helper: posterior-interval annotation for plot
.make_post_ann <- function(x, w) {
  w <- w / sum(w)
  m  <- sum(w * x)
  ci <- as.numeric(quantile(x, c(0.025, 0.975)))
  sprintf("mean = %.3f\n95%% CI [%.3f, %.3f]", m, ci[1], ci[2])
}

#' @export
plot.spima <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }

  post <- x$abc_result$posterior
  df_theta <- as.data.frame(post$theta)
  df_theta$weight <- post$weights

  pnames <- setdiff(names(df_theta), "weight")

  PAL  <- "#0072B2"   # blue  (Okabe-Ito)
  PAL2 <- "#D55E00"   # orange

  for (nm in pnames) {
    ann <- .make_post_ann(df_theta[[nm]], df_theta$weight)

    p <- ggplot2::ggplot(df_theta,
                         ggplot2::aes(x = .data[[nm]],
                                      weight = .data[["weight"]])) +
      ggplot2::geom_histogram(bins = 40, fill = PAL, colour = NA,
                              alpha = 0.88) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = ann, hjust = 1.05, vjust = 1.4,
                        size = 3.2, colour = "grey30",
                        fontface = "plain") +
      ggplot2::labs(x = nm, y = "Weighted count") +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 13),
        axis.text  = ggplot2::element_text(size = 11),
        plot.margin = ggplot2::margin(12, 14, 10, 10)
      )
    suppressWarnings(print(p))
  }

  if (length(pnames) == 2) {
    p2 <- ggplot2::ggplot(df_theta,
                          ggplot2::aes(x = .data[[pnames[1]]],
                                       y = .data[[pnames[2]]],
                                       size = .data[["weight"]])) +
      ggplot2::geom_point(alpha = 0.55, colour = PAL2, stroke = 0) +
      ggplot2::scale_size_continuous(guide = "none") +
      ggplot2::labs(x = pnames[1], y = pnames[2]) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 13),
        axis.text  = ggplot2::element_text(size = 11),
        plot.margin = ggplot2::margin(12, 14, 10, 10)
      )
    suppressWarnings(print(p2))
  }

  invisible(x)
}

#' @export
print.spima_abc <- function(x, ...) {
  cat("spima ABC-SMC result\n")
  cat("Generations:", length(x$generations), "\n")
  cat("Final epsilon:", signif(x$posterior$epsilon, 4), "\n")
  cat("Effective sample size:", round(x$posterior$ess, 1), "\n")
  invisible(x)
}

# -------- Subgroup S3 methods --------

#' @export
print.spima_subgroup <- function(x, ...) {
  cat("spima subgroup analysis (", x$outcome_type, ")\n", sep = "")
  cat("  Subgroup variable:", x$subgroup_col, "\n")
  cat("  Number of subgroups:", length(x$results), "\n\n")

  for (nm in names(x$results)) {
    cat("--- ", nm, " ---\n", sep = "")
    r <- x$results[[nm]]
    gen <- r$generations
    last <- gen[[length(gen)]]
    cat("  Generations:", length(gen), "  Particles:", nrow(last$theta),
        "  ESS:", round(last$ess, 1), "\n")
    for (pn in names(r$summary)) {
      cat("  ", pn, ": mean = ", signif(r$summary[[pn]]["mean"], 4),
          "  sd = ", signif(r$summary[[pn]]["sd"], 4),
          "  95% CI [", signif(r$summary[[pn]]["q2.5"], 4), ", ",
                       signif(r$summary[[pn]]["q97.5"], 4), "]\n", sep = "")
    }
    cat("\n")
  }
  invisible(x)
}

#' @export
plot.spima_subgroup <- function(x, parameter = NULL, ...) {
  first <- x$results[[1]]$summary
  avail <- names(first)

  if (is.null(parameter)) {
    parameter <- if ("mu" %in% avail) "mu" else avail[1]
  }
  if (!parameter %in% avail) {
    stop("Parameter '", parameter, "' not found. Available: ",
         paste(avail, collapse = ", "))
  }

  subgroup_names <- names(x$results)
  est <- sapply(subgroup_names, function(nm) x$results[[nm]]$summary[[parameter]]["mean"])
  lo  <- sapply(subgroup_names, function(nm) x$results[[nm]]$summary[[parameter]]["q2.5"])
  hi  <- sapply(subgroup_names, function(nm) x$results[[nm]]$summary[[parameter]]["q97.5"])

  k <- length(est)
  y <- k:1

  xlim <- range(c(lo, hi, 0))
  xpad <- diff(xlim) * 0.35
  xlim <- xlim + c(-xpad * 0.6, xpad * 0.4)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(
    mar = c(4.2, 7, 3, 2),
    family = "sans",
    cex.axis = 1.05,
    cex.lab  = 1.1
  )

  plot(est, y, xlim = xlim, ylim = c(0.3, k + 0.7),
       yaxt = "n", xlab = "", ylab = "",
       main = "",
       pch = 22, cex = 1.5, bg = "#0072B2", col = "#0072B2",
       frame.plot = FALSE, ...)

  graphics::title(
    main = paste("Forest plot of", parameter, "by", x$subgroup_col),
    xlab = parameter, line = 2.5
  )

  graphics::axis(2, at = y, labels = subgroup_names, las = 1,
                 tick = FALSE, line = -0.5)

  # reference at 0
  graphics::abline(v = 0, lty = "13", col = "grey40", lwd = 0.8)

  # CI line segments
  graphics::segments(lo, y, hi, y, lwd = 2.2, col = "#0072B2")
  graphics::segments(lo, y - 0.12, lo, y + 0.12, lwd = 2.0, col = "#0072B2")
  graphics::segments(hi, y - 0.12, hi, y + 0.12, lwd = 2.0, col = "#0072B2")

  # annotation on the right
  ann <- sprintf("%.3f  [ %.3f ,  %.3f ]", est, lo, hi)
  ann_x <- graphics::grconvertX(0.97, "nfc", "user")
  graphics::text(ann_x, y, labels = ann, adj = 0, cex = 0.82,
                 col = "grey20", family = "sans")

  invisible(x)
}

# -------- as.data.frame methods for result table output --------

#' @export
as.data.frame.spima <- function(x, ...) {
  s <- x$abc_result$summary
  pnames <- names(s)
  tab <- data.frame(
    parameter = pnames,
    mean = sapply(pnames, function(p) s[[p]]["mean"]),
    sd   = sapply(pnames, function(p) s[[p]]["sd"]),
    q2.5 = sapply(pnames, function(p) s[[p]]["q2.5"]),
    q50  = sapply(pnames, function(p) s[[p]]["q50"]),
    q97.5 = sapply(pnames, function(p) s[[p]]["q97.5"]),
    row.names = NULL, stringsAsFactors = FALSE
  )
  tab
}

#' @export
as.data.frame.spima_subgroup <- function(x, ...) {
  rows <- list()
  for (nm in names(x$results)) {
    s <- x$results[[nm]]$summary
    for (pn in names(s)) {
      rows[[length(rows) + 1]] <- data.frame(
        subgroup = nm, parameter = pn,
        mean = s[[pn]]["mean"], sd = s[[pn]]["sd"],
        q2.5 = s[[pn]]["q2.5"], q50 = s[[pn]]["q50"],
        q97.5 = s[[pn]]["q97.5"],
        row.names = NULL, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# -------- Multivariate normal sampler (avoid MASS dependency) --------

rmvn <- function(n, mu, Sigma) {
  p <- length(mu)
  eig <- eigen(Sigma, symmetric = TRUE)
  if (any(eig$values < 0)) {
    eig$values <- pmax(eig$values, 1e-10)
  }
  A <- eig$vectors %*% diag(sqrt(eig$values), p)
  matrix(mu, n, p, byrow = TRUE) +
    matrix(rnorm(n * p), n, p) %*% t(A)
}

# Replace MASS::mvrnorm in abc_smc.R with rmvn
# (Note: abc_smc.R was written using MASS::mvrnorm; see the perturb_particle function)
