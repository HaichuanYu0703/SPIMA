# Distance function tests for spima
# Run: Rscript tests/testthat/test-distances.R

library(spima)

check <- function(expr, msg) {
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) FALSE)
  if (ok) cat("  PASS:", msg, "\n") else cat("  FAIL:", msg, "\n")
  invisible(ok)
}

# ---- spima_bin_distance ----
cat("--- spima_bin_distance ---\n")
s <- c(A=-0.5, B=0.3, C=0.8)
o <- structure(c(A=-0.3, B=0.2, C=0.9), weights = c(A=4, B=2, C=1))
d <- spima_bin_distance(s, o)
check(is.finite(d) && d > 0, "returns finite positive distance")
check(spima_bin_distance(c(A=1), c(A=1)) == 0, "zero distance when identical")

# No overlap
check(is.infinite(spima_bin_distance(c(D=1), c(A=0))), "infinite when no overlap")

# ---- spima_generic_distance ----
cat("--- spima_generic_distance ---\n")
s2 <- c(A=0.5, B=0.3)
o2 <- structure(c(A=0.4, B=0.2), weights = c(A=10, B=5))
d2 <- spima_generic_distance(s2, o2)
check(is.finite(d2) && d2 > 0, "returns finite positive distance")
check(spima_generic_distance(c(A=1), c(A=1)) == 0, "zero distance when identical")
o3 <- c(A=0.4, B=0.2)  # no weights attribute
d3 <- spima_generic_distance(s2, o3)
check(is.finite(d3) && d3 > 0, "fallback to default weights when no weight attr")

# ---- spima_cont_distance ----
cat("--- spima_cont_distance ---\n")
sim <- list(means = c(A=1.0, B=2.0), sds = c(A=1.5, B=2.0))
obs <- list(means = c(A=0.5, B=1.5), sds = c(A=1.0, B=1.5))
d4 <- spima_cont_distance(sim, obs)
check(is.finite(d4) && d4 > 0, "returns finite positive distance")

# With precision weights attribute
obs_w <- list(means = c(A=0.5, B=1.5), sds = structure(c(A=1.0, B=1.5), weights = c(A=10, B=5)))
d5 <- spima_cont_distance(sim, obs_w)
check(is.finite(d5) && d5 > 0, "uses precision weights when available")

# No overlap
check(is.infinite(spima_cont_distance(list(means=c(D=1), sds=c(D=1)), list(means=c(A=0), sds=c(A=1)))),
       "infinite when no overlap")

# ---- Edge cases ----
cat("--- Edge cases ---\n")
check(is.infinite(spima_bin_distance(numeric(0), numeric(0))), "empty input -> Inf")
check(is.infinite(spima_generic_distance(numeric(0), numeric(0))), "empty generic -> Inf")
check(is.infinite(spima_cont_distance(list(means=numeric(0), sds=numeric(0)),
                                        list(means=numeric(0), sds=numeric(0)))),
       "empty cont -> Inf")

# Negative values
check(is.finite(spima_bin_distance(c(X=-2), structure(c(X=-1), weights=c(X=1)))),
       "negative values handled")

cat("\nDone.\n")
