# EZbakR 1.0.0

While technically backwards compatible, this release radically reworks many of the
internals of EZbakR to improve computational and statistical performance, while also
adding several pieces of brand new functionality.

### New features
* **Goodness-of-fit diagnostics (`AssessGoF()`, `EZGoFPlot()`, `EZFitPlot()`)**:
  * Evaluate how well your data fit EZbakR's mutational mixture models (`EstimateFractions()`, `EstimateIsoformFractions()`, `DeconvolveFractions()`) and kinetic/ODE models (`AverageAndRegularize()`, `EZDynamics()`).
  * Flags common real-world issues like uncalled SNPs, cross-mapping paralog reads, mixed pre-mRNA/mRNA signals, outlier samples, and misspecified kinetic models.
  * Added `EZGoFPlot()` for dataset- and sample-wide QC plots, `EZFitPlot()` for inspecting individual feature fits, and a new [Goodness of fit vignette](https://isaacvock.github.io/EZbakR/articles/AssessGoF.html).
* **Easier Apache Arrow workflows (`WriteCBArrow()`, `CollectArrow()`, `ReconnectArrow()`)**:
  * `EZbakRArrowData()` now accepts folder paths, Parquet file lists, and filtered Arrow queries directly, and handles saving/reloading `.rds` objects cleanly via `ReconnectArrow()`.
  * Added `WriteCBArrow()` to convert `cB` tables or CSVs into partitioned Parquet datasets, and `CollectArrow()` to pull specific samples or features back into memory.
* **Multi-feature read filtering**: Added `remove_multi_features = TRUE` in `EstimateFractions()` to optionally drop reads ambiguously assigned to multiple features.

### Improvements & bug fixes
* **Faster runtimes across the board**: Substantially sped up `EstimateFractions()` (standard, hierarchical, and multi-label models), `EstimateIsoformFractions()`, `DeconvolveFractions()`, `AverageAndRegularize()`, `NormalizeForDropout()`, and `EZDynamics()`.
* **Progress updates for ODE fitting**: `EZDynamics()` now prints progress and estimated remaining runtime (`verbose = TRUE`).
* **More accurate fractions & uncertainties**:
  * `EstimateFractions()` now adapts its prior automatically to your experiment's overall labeling level (`EBayes = TRUE`, reducing bias at very short or long label times) and reports 95% confidence intervals (`ci_lower_*`, `ci_upper_*`).
  * Improved uncertainty propagation and extreme-value (`fn = 0` or `1`) handling across `EstimateKinetics()`, `AverageAndRegularize()`, and `CompareParameters()`.
* **Bug fixes & docs**: Improved numerical stability in `CorrectDropout()`, `NormalizeForDropout()`, and `EZDynamics()`, and cleaned up typos and examples across all vignettes.

# EZbakR 0.2.2

- Patch to address CRAN comments as well as hierarchical fraction new modeling
issue ([#33](https://github.com/isaacvock/EZbakR/issues/33))

# EZbakR 0.2.1

- Patch to address reverse compatibility notes from new arrow release (25.0.0)

# EZbakR 0.2.0

- Implemented multi-label time simulations and documented these

# EZbakR 0.1.0

* Initial CRAN submission.

- EZbakR is on CRAN now.
- First fully stable release, though much of the functionality has been around
for several months now.
- Blog post on release [here](https://isaacvock.github.io/posts/ezbakR_update/)
