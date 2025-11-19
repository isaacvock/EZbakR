# Get contrasts of estimated parameters

`CompareParameters()` calculates differences in parameters estimated by
[`AverageAndRegularize()`](https://isaacvock.github.io/EZbakR/reference/AverageAndRegularize.md)
or
[`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md)
and performs null hypothesis statistical testing, comparing their values
to a null hypothesis of 0.

## Usage

``` r
CompareParameters(
  obj,
  design_factor,
  reference,
  experimental,
  param_name,
  parameter = "log_kdeg",
  type = "averages",
  param_function,
  condition = NULL,
  features = NULL,
  exactMatch = TRUE,
  repeatID = NULL,
  formula_mean = NULL,
  sd_grouping_factors = NULL,
  fit_params = NULL,
  normalize_by_median = FALSE,
  reference_levels = NULL,
  experimental_levels = NULL,
  overwrite = TRUE
)
```

## Arguments

- obj:

  An `EZbakRFit` object, which is an `EZbakRFractions` object on which
  [`AverageAndRegularize()`](https://isaacvock.github.io/EZbakR/reference/AverageAndRegularize.md)
  has been run.

- design_factor:

  Name of factor from `metadf` whose parameter estimates at different
  factor values you would like to compare. If you specify this, you need
  to also specify `reference` and `experimental`. If `type` ==
  "dynamics", this can have multiple values, being the names of all of
  the factors you would like to stratify a group by.

- reference:

  Name of reference `design_factor` factor level value. Difference will
  be calculated as `experimental` - `reference`. If type == "dynamics",
  then this should specify the levels of all of the `design_factor`(s)
  reference group. For example, if you have multiple `design_factor`'s,
  then `reference` must be a named character vector with one element per
  `design_factor`, with elements named the corresponding
  `design_factor`. For example, if `design_factor` is c("genotype",
  "treatment"), and you would like to compare genotype = "WT" and
  treatment = "untreated" (reference) to genotype = "KO" and treatment =
  "treated", then `reference` would need to be a vector with one element
  named "genotype", equal to "WT" and one element named "treatment"
  equal to "untreated" (this example could be created with,
  `c(genotype = "WT", treatment = "untreated")`).

- experimental:

  Name of `design_factor` factor level value to compare to reference.
  Difference will be calculated as `experimental` - `reference`. If type
  == "dynamics", then this should specify the levels of all of the
  `design_factor`(s) reference group. For example, if you have multiple
  `design_factor`'s, then `experimental` must be a named character
  vector with one element per `design_factor`, with elements named the
  corresponding `design_factor`. For example, if `design_factor` is
  c("genotype", "treatment"), and you would like to compare genotype =
  "WT" and treatment = "untreated" (reference) to genotype = "KO" and
  treatment = "treated", then `experimental` would need to be a vector
  with one element named "genotype", equal to "KO" and one element named
  "treatment" equal to "treated" (this example could be created with,
  `c(genotype = "KO", treatment = "treated")`).

- param_name:

  If you want to assess the significance of a single parameter, rather
  than the comparison of two parameters, specify that one parameter's
  name here.

- parameter:

  Parameter to average across replicates of a given condition.

- type:

  Type of table to use. Can either be "averages" or "dynamics".

- param_function:

  NOT YET IMPLEMENTED. Will allow you to specify more complicated
  functions of parameters when hypotheses you need to test are
  combinations of parameters rather than individual parameters or simple
  differences in two parameters.

- condition:

  Same as `design_factor`, will be deprecated in favor of the former in
  later release.

- features:

  Character vector of the set of features you want to stratify reads by
  and estimate proportions of each RNA population. The default of "all"
  will use all feature columns in the `obj`'s cB.

- exactMatch:

  If TRUE, then `features` and `populations` have to exactly match those
  for a given fractions table for that table to be used. Means that you
  can't specify a subset of features or populations by default, since
  this is TRUE by default.

- repeatID:

  If multiple `averages` tables exist with the same metadata, then this
  is the numerical index by which they are distinguished.

- formula_mean:

  An R formula object specifying how the `parameter` of interest depends
  on the sample characteristics for the averages object you want to use.

- sd_grouping_factors:

  Metadf columns should data was grouped by when estimating standard
  deviations across replicates for the averages object you want to use.

- fit_params:

  Character vector of parameter names in the averages object you would
  like to use.

- normalize_by_median:

  If TRUE, then median difference will be set equal to 0. This can
  account for global biases in parameter estimates due to things like
  differences in effective label times. Does risk eliminating real
  signal though, so user discretion is advised.

- reference_levels:

  Same as `reference`, but exclusively parsed in case of `type` ==
  "dynamics, included for backwards compatibility.

- experimental_levels:

  Same as `experimental`, but exclusively parsed in case of `type` ==
  "dynamics, included for backwards compatibility.

- overwrite:

  If TRUE, then identical output will be overwritten if it exists.

## Value

`EZbakRData` object with an additional "comparisons" table, detailing
the result of a comparison of a parameter estimate's valules across two
different conditions.

## Details

The EZbakR website has an extensive vignette walking through various use
cases and parameters you can compare with `CompareParameters()`:
[vignette
link](https://isaacvock.github.io/EZbakR/articles/Linear-modeling.html).

There are essentially 3 scenarios that `CompareParameters()` can handle:

- Pairwise comparisons: compare `reference` to `experimental` parameter
  estimates of a specified `design_factor` from
  [`AverageAndRegularize()`](https://isaacvock.github.io/EZbakR/reference/AverageAndRegularize.md).
  log(`experimental` / `reference`) is the computed "difference" in this
  case.

- Assess the value of a single parameter, which itself should represent
  a difference between other parameters. The name of this parameter can
  be specified via the `param_name` argument. This is useful for various
  interaction models where some of the parameters of these models may
  represent things like "effect of A on condition X".

- Pairwise comparison of dynamical systems model parameter estimate:
  similar to the first case listed above, but now when
  `type == "dynamics"`. `design_factor` can now be a vector of all the
  `metadf` columns you stratified parameter estimates by.

Eventually, a 4th option via the currently non-functional
`param_function` argument will be implemented, which will allow you to
specify functions of parameters to be assessed, which can be useful for
certain interaction models.

`CompareParameters()` calculates p-values using what is essentially an
asymptotic Wald test, meaning that a standard normal distribution is
integrated. P-values are then multiple-test adjusted using the method of
Benjamini and Hochberg, implemented in R's
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) function.

## Examples

``` r
# Simulate data to analyze
simdata <- EZSimulate(30)

# Create EZbakR input
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)

# Estimate Fractions
ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output

# Estimate Kinetics
ezbdo <- EstimateKinetics(ezbdo)

# Average estimates across replicate
ezbdo <- AverageAndRegularize(ezbdo)
#> Fitting linear model
#> Estimating coverage vs. variance trend
#> Regularizing variance estimates

# Compare parameters across conditions
ezbdo <- CompareParameters(
ezbdo,
design_factor = "treatment",
reference = "treatment1",
experimental = "treatment2"
)
```
