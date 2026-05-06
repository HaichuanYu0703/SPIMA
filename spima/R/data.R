#' Blood Pressure Continuous Outcome Data
#'
#' A dataset of study-level summary statistics for continuous outcomes
#' (blood pressure) from multiple clinical trials.  Contains mean, standard
#' deviation, and sample size per arm, suitable for the continuous module.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{group}{Treatment group indicator (0 = control, 1 = treatment).}
#'   \item{n}{Sample size per arm.}
#'   \item{mean}{Mean blood pressure.}
#'   \item{sd}{Standard deviation of blood pressure.}
#' }
#' @keywords data
"bp_cont"

#' Kidney Disease Binary Outcome Data
#'
#' A dataset of study-level summary statistics for binary outcomes
#' (kidney disease) from multiple clinical trials.  Contains event counts
#' and sample sizes per arm, suitable for the binary module.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{group}{Treatment group indicator (0 = control, 1 = treatment).}
#'   \item{n}{Sample size per arm.}
#'   \item{event}{Number of events per arm.}
#' }
#' @keywords data
"kidney_bin"

#' Generic Effect Size Data
#'
#' A dataset of study-level summary statistics for generic (continuous)
#' effect sizes.  Contains effect size estimates and their standard errors,
#' suitable for the generic module.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{yi}{Effect size estimate.}
#'   \item{sei}{Standard error of the effect size.}
#' }
#' @keywords data
"gen_effect"
