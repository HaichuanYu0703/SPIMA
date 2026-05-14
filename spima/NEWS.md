# spima 0.2.0

* **Distance function fix** (v0.2.0): `spima_cont_distance` changed from relative-difference + SD mixing to inverse-variance weighted Euclidean distance, eliminating systematic bias in the continuous module.
* **Mixed-format support**: all 5 input formats (mean+SD, median+IQR, median+range, mean+range, five-number summary) can be mixed in a single dataset.
* **Forest plot**: new `forest.spima()` S3 method with `log_scale` and custom label support. Alias `spima_forest()` works when `metafor` masks the generic.
* **Adaptive particles**: `smc_control()` gains `n_particles_max`; ABC-SMC doubles particle count when acceptance rate is low.
* **Parallel fallback**: graceful reversion to sequential when PSOCK cluster fails.
* **Precision weights**: `spima_cont_observed_stats` now attaches proper SE-based weights for more accurate distance computation.
* **Improved `as.data.frame()`**: new `probs` parameter for custom posterior quantiles.

# spima 0.1.0

* Initial release.
* Meta-analysis via ABC-SMC by simulating pseudo-individual data from published group-level summary statistics.
* Modules: binary, continuous (Gaussian and Gamma, with flexible input formats), generic effect-size, and interaction analysis.
* C++ accelerated simulation and distance computation.
