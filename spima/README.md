# spima

**Simulated Pseudo-Individual Data Meta-Analysis with ABC-SMC**

spima performs meta-analysis via Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) by simulating pseudo-individual data from published group-level summary statistics. It handles binary, continuous, and generic effect-size outcomes within a one-stage mixed-model framework.

## Installation

```r
# From GitHub (latest)
install.packages("remotes")
remotes::install_github("HaichuanYu0703/SPIMA")

# Or from local source
install.packages("path/to/spima", repos = NULL, type = "source")
```

### Requirements

* R (>= 3.5.0)
* C++ compiler (GNU make, Rtools on Windows)
* Packages: lme4, Rcpp, RcppArmadillo

## Quick Start

```r
library(spima)

# ===== Binary outcome =====
data_bin <- data.frame(
  study = 1:4, group = rep(c("c", "t"), each = 4),
  event = c(30, 100, 45, 100, 28, 80, 32, 80),
  n     = c(100, 100, 80, 80, 80, 80, 60, 60)
)
res <- spima(data_bin, "binary",
  input_spec = list(study = "study", group = "group",
                    event = "event", n = "n"),
  prior = prior(mu = "normal(0, 2.5)", tau = "halfnormal(0, 0.5)"))
print(res)
forest(res)        # forest plot
```

## Mixed-Format Data (v0.2.0+)

spima supports mixing different reporting formats in a single analysis.
Just include all possible columns; un-reported values go as `NA`:

```r
dat <- data.frame(
  study = c("A","A","B","B","C","C","D","D","E","E"),
  group = c("c","t","c","t","c","t","c","t","c","t"),
  n     = c(100,100,100,100,100,100,100,100,100,100),
  # A: mean + SD
  mean  = c(50,52, NA,NA, NA,NA, NA,NA, 48.5,51.0),
  sd    = c(10,10.5, NA,NA, NA,NA, NA,NA, NA,NA),
  # B: median + IQR
  median = c(NA,NA, 48,50.5, 49,51, 47.5,50, NA,NA),
  q1     = c(NA,NA, 42,44.5, 43,45, NA,NA, NA,NA),
  q3     = c(NA,NA, 55,57, 56,58, NA,NA, NA,NA),
  # C & D: full range
  min    = c(NA,NA, NA,NA, 30,32, 35,38, 35,37),
  max    = c(NA,NA, NA,NA, 70,72, 65,68, 62,67)
)
res <- spima(dat, "continuous",
  input_spec = list(study = "study", group = "group",
    mean = "mean", sd = "sd", median = "median",
    q1 = "q1", q3 = "q3", min = "min", max = "max", n = "n"))
print(res)
```

## Modules

| Module | Description |
|--------|-------------|
| binary | Binary outcome (event/n), log-odds ratio |
| continuous | Continuous (mean/SD, median/IQR, range, 5-number, mean+range) |
| generic | Pre-computed effect sizes (yi/sei) |
| interaction | Covariate-treatment interaction analysis |

## Documentation

```r
?spima
?prior
?smc_control
?forest.spima
```
