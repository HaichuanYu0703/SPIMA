# spima

**Simulated Pseudo-Individual Data Meta-Analysis with ABC-SMC**

spima performs meta-analysis via Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) by simulating pseudo-individual data from published group-level summary statistics. It handles binary, continuous, and generic effect-size outcomes within a one-stage mixed-model framework.

## Installation

```r
# Install from local source
install.packages("path/to/spima", repos = NULL, type = "source")
```

### Requirements

* R (>= 3.5.0)
* C++ compiler (GNU make, Rtools on Windows)
* Packages: lme4, Rcpp, RcppArmadillo

## Quick Start

```r
library(spima)

# Binary outcome meta-analysis
data_bin <- data.frame(
  study = 1:4,
  event_t = c(45, 32, 58, 22),
  n_t     = c(100, 80, 120, 60),
  event_c = c(30, 28, 40, 18),
  n_c     = c(100, 80, 120, 60)
)

res <- spima(data_bin, "binary",
  input_spec = list(study = "study", event = "event_t",
                    n = "n_t", group = "group"),
  prior = prior(mu = "normal(0, 10)", tau = "halfnormal(0, 1)"),
  smc_control = smc_control(n_particles = 500, n_generations = 5))

print(res)
as.data.frame(res)
```

## Modules

| Module | Description |
|--------|-------------|
| binary | Binary outcome (event/n), log-odds ratio |
| continuous | Continuous outcome (mean/SD, median/IQR, range, 5-number) |
| generic | Pre-computed effect sizes (yi/sei) |
| interaction | Covariate-treatment interaction analysis |

## Documentation

See the package vignettes and function documentation for details:

```r
?spima
?prior
?smc_control
?spima_int
```
