#' Define Prior Distributions for ABC-SMC Parameters
#'
#' @param mu Prior specification for overall effect \eqn{\mu},
#'   e.g. \code{"normal(0, 10)"}.
#' @param tau Prior specification for heterogeneity \eqn{\tau},
#'   e.g. \code{"halfnormal(0, 1)"}.
#' @param ... Additional named priors (e.g. \code{gamma = "uniform(0, 5)"}).
#' @return A list of class \code{spima_prior} with elements \code{name},
#'   \code{pars}, \code{rfun} (random generation), \code{dfun} (density),
#'   and \code{default}.
#' @export
#'
#' @examples
#' prior(mu = "normal(0, 10)", tau = "halfnormal(0, 1)")
prior <- function(mu = "normal(0, 10)", tau = "halfnormal(0, 1)", ...) {
  extras <- list(...)

  # Detect outcome-specific parameters (mu1, tau1, etc.)
  # When present, drop defaults to avoid name collisions.
  extra_names <- names(extras)
  has_numeric <- any(grepl("^(mu|tau)[0-9]", extra_names))

  if (has_numeric) {
    specs <- extras
  } else {
    specs <- c(list(mu = mu, tau = tau), extras)
  }

  out <- list()
  for (nm in names(specs)) {
    out[[nm]] <- parse_prior_spec(nm, specs[[nm]])
  }

  structure(out, class = "spima_prior")
}

parse_prior_spec <- function(name, spec) {
  spec <- gsub("\\s+", "", spec)

  # halfnormal(0, 1)
  if (grepl("^halfnormal\\(", spec, ignore.case = TRUE)) {
    inner <- gsub("^halfnormal\\(|\\)$", "", spec, ignore.case = TRUE)
    parts <- as.numeric(strsplit(inner, ",")[[1]])

    if (length(parts) != 2) stop("halfnormal requires 2 parameters: location, scale")

    list(
      name     = "halfnormal",
      pars     = c(location = parts[1], scale = parts[2]),
      rfun     = function(n, ...) {
        abs(rnorm(n, mean = parts[1], sd = parts[2]))
      },
      dfun     = function(x, ...) {
        2 * dnorm(x, mean = parts[1], sd = parts[2]) * (x >= parts[1])
      },
      default  = parts
    )

  # normal(0, 10)
  } else if (grepl("^normal\\(", spec, ignore.case = TRUE)) {
    inner <- gsub("^normal\\(|\\)$", "", spec, ignore.case = TRUE)
    parts <- as.numeric(strsplit(inner, ",")[[1]])

    if (length(parts) != 2) stop("normal requires 2 parameters: mean, sd")

    list(
      name     = "normal",
      pars     = c(mean = parts[1], sd = parts[2]),
      rfun     = function(n, ...) rnorm(n, mean = parts[1], sd = parts[2]),
      dfun     = function(x, ...) dnorm(x, mean = parts[1], sd = parts[2]),
      default  = parts
    )

  # uniform(0, 5)
  } else if (grepl("^uniform\\(", spec, ignore.case = TRUE)) {
    inner <- gsub("^uniform\\(|\\)$", "", spec, ignore.case = TRUE)
    parts <- as.numeric(strsplit(inner, ",")[[1]])

    if (length(parts) != 2) stop("uniform requires 2 parameters: min, max")

    list(
      name     = "uniform",
      pars     = c(min = parts[1], max = parts[2]),
      rfun     = function(n, ...) runif(n, min = parts[1], max = parts[2]),
      dfun     = function(x, ...) dunif(x, min = parts[1], max = parts[2]),
      default  = parts
    )

  # lnorm(0, 1)
  } else if (grepl("^lnorm\\(", spec, ignore.case = TRUE) ||
             grepl("^lognormal\\(", spec, ignore.case = TRUE)) {
    inner <- gsub("^(lnorm|lognormal)\\(|\\)$", "", spec, ignore.case = TRUE)
    parts <- as.numeric(strsplit(inner, ",")[[1]])

    if (length(parts) != 2) stop("lnorm requires 2 parameters: meanlog, sdlog")

    list(
      name     = "lnorm",
      pars     = c(meanlog = parts[1], sdlog = parts[2]),
      rfun     = function(n, ...) rlnorm(n, meanlog = parts[1], sdlog = parts[2]),
      dfun     = function(x, ...) dlnorm(x, meanlog = parts[1], sdlog = parts[2]),
      default  = parts
    )

  } else {
    stop("Unrecognised prior specification: ", spec,
         ". Use: normal, halfnormal, uniform, lnorm")
  }
}

#' Evaluate log-prior density for a parameter vector
#'
#' @param theta Named numeric vector of parameters.
#' @param prior_obj A \code{spima_prior} object.
#' @return Log-density value (summed across independent priors).
log_prior_density <- function(theta, prior_obj) {
  lp <- 0
  for (nm in names(theta)) {
    if (nm %in% names(prior_obj)) {
      lp <- lp + log(prior_obj[[nm]]$dfun(theta[nm]))
    }
  }
  lp
}

#' Sample from the joint prior
#'
#' @param n Number of samples.
#' @param prior_obj A \code{spima_prior} object.
#' @return A matrix with \code{n} rows and one column per prior.
sample_prior <- function(n, prior_obj) {
  m <- matrix(NA, nrow = n, ncol = length(prior_obj))
  colnames(m) <- names(prior_obj)
  for (nm in names(prior_obj)) {
    m[, nm] <- prior_obj[[nm]]$rfun(n)
  }
  m
}
