# Average parameter estimates across replicates, and regularize variance estimates

`AverageAndRegularize` fits a generalized linear model to your data to
effectively average parameter estimates across replicates and get
overall uncertainty estimates for those parameters. The linear model to
which your data is fit is specified via an R formula object supplied to
the `formula_mean` parameter. Uncertainty estimates are regularized via
a hierarchical modeling strategy originally introduced with bakR, though
slightly improved upon since then.

## Usage

``` r
AverageAndRegularize(
  obj,
  features = NULL,
  parameter = "log_kdeg",
  type = "kinetics",
  kstrat = NULL,
  populations = NULL,
  fraction_design = NULL,
  exactMatch = TRUE,
  repeatID = NULL,
  formula_mean = NULL,
  sd_grouping_factors = NULL,
  include_all_parameters = TRUE,
  sd_reg_factor = 10,
  error_if_singular = TRUE,
  min_reads = 10,
  convert_tl_to_factor = TRUE,
  regress_se_with_abs = TRUE,
  force_lm = FALSE,
  force_optim = force_lm,
  conservative = FALSE,
  character_limit = 20,
  feature_lengths = NULL,
  feature_sample_counts = NULL,
  scale_factor_df = NULL,
  overwrite = TRUE
)
```

## Arguments

- obj:

  An `EZbakRFractions` or `EZbakRKinetics` object, which is an
  `EZbakRData` object on which
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md)
  or
  [`EstimateKinetics()`](https://isaacvock.github.io/EZbakR/reference/EstimateKinetics.md)
  has been run.

- features:

  Character vector of the set of features you want to stratify reads by
  and estimate proportions of each RNA population. The default of "all"
  will use all feature columns in the `obj`'s cB.

- parameter:

  Parameter to average across replicates of a given condition.

- type:

  What type of table is the parameter found in? Default is "kinetics",
  but can also set to "fractions".

- kstrat:

  If `type == "kinetics"`, then `kstrat` specifies the kinetic parameter
  inference strategy.

- populations:

  Character vector of the set of mutational populations that you want to
  infer the fractions of. Only relevant if type == "fractions".

- fraction_design:

  "Design matrix" specifying which RNA populations exist in your
  samples. Only relevant if type == "fractions".

- exactMatch:

  If TRUE, then `features` and `populations` have to exactly match those
  for a given fractions table for that table to be used. Means that you
  can't specify a subset of features or populations by default, since
  this is TRUE by default.

- repeatID:

  If multiple `kinetics` or `fractions` tables exist with the same
  metadata, then this is the numerical index by which they are
  distinguished.

- formula_mean:

  An R formula object specifying how the `parameter` of interest depends
  on the sample characteristics specified in `obj`'s metadf. The most
  common formula will be `~ treatment` or `~ treatment:duration`, where
  `treatment` and `duration` would be replaced with whatever you called
  the relevant sample characteristics in your metadf. `~ treatment`
  means that an average value of `parameter` should be estimated for
  each set of samples with the same value for `treatment` in the metadf.
  `~ treatment:duration` specifies that an average value of `parameter`
  should be estimated for each set of samples with the same combination
  of `treatment` and `duration` values in the metadf. An example of the
  latter case is a situation where you have two or more treatments
  (e.g., drug treated and untreated control) which were applied for
  different durations of time (e.g., 4 and 8 hours).

  NOTE: EZbakR automatically removes any intercept terms from the model.
  That way, there is no ambiguity about what parameter is defined as the
  reference.

- sd_grouping_factors:

  What metadf columns should data be grouped by when estimating standard
  deviations across replicates? If this is NULL, then EZbakR will check
  to see if the `formula_mean` specifies a formula that cleanly
  stratifies samples into disjoint groups. For example, the formula
  `~ treatment` will assign each sample to a single factor (its value
  for the metadf's `treatment` column). In this case, standard
  deviations can be calculated for sets of replicates in each
  `treatment` group. If such a stratification does not exist, a single
  standard deviation will be estimated for each feature (i.e.,
  homoskedasticity will be assumed as in standard linear modeling).

- include_all_parameters:

  If TRUE, an additional table will be saved with the prefix `fullfit_`,
  which includes all of the parameters estimated throughout the course
  of linear modeling and regularization. This can be nice for
  visualizing the regularized mean-variance trend.

- sd_reg_factor:

  Determines how strongly variance estimates are shrunk towards trend.
  Higher numbers lead to more regularization. Eventually, this will be
  replaced with estimation of how much variance there seems to be in the
  population of variances.

- error_if_singular:

  If TRUE, linear model will throw an error if parameters cannot be
  uniquely identified. This is most often caused by parameters that
  cannot be estimated from the data, e.g., due to limited replicate
  numbers or correlated sample characteristics (i.e., all treatment As
  also correspond to batch As, and all treatment Bs correspond to batch
  Bs).

- min_reads:

  Minimum number of reads in all samples for a feature to be kept.

- convert_tl_to_factor:

  If a label time variable is included in the `formula_mean`, convert
  its values to factors so as to avoid performing continuous regression
  on label times. Defaults to TRUE as including label time in the
  regression is often meant to stratify samples by their label time if,
  for example, you are averaging logit(fractions).

- regress_se_with_abs:

  If TRUE, and if `type == "fractions"`, then standard error will be
  regressed against logit fraction rather than magnitude of logit
  fraction. Makes sense to set this to FALSE if analyzing certain
  site-specific mutational probing methods when high mutation content
  things are likely low variance SNPs.

- force_lm:

  Certain formula lend them selves to efficient approximations of the
  full call to [`lm()`](https://rdrr.io/r/stats/lm.html). Namely,
  formulas that stratify samples into disjoint groups where a single
  parameter of the model is effectively estimated from each group can be
  tackled via simple averaging of data from each from group. If you
  would like to force EZbakR to fit the fully rigorous linear model
  though, set `force_lm` to `TRUE`.

- force_optim:

  Old parameter that is now passed the value `force_lm` and will be
  deprecated in later releases

- conservative:

  If TRUE, conservative variance regularation will be performed. In this
  case, variances below the trend will be regularized up to the trend,
  and variances above the trend will be left unregularized. This avoids
  underestimation of variances.

- character_limit:

  Limit on the number of characters of the name given to the output
  table. Will attempt to concatenate the parameter name with the names
  of all of the features. If this is too long, only the parameter name
  will be used.

- feature_lengths:

  Table of effective lengths for each feature combination in your data.
  For example, if your analysis includes features named GF and XF, this
  should be a data frame with columns GF, XF, and length.

- feature_sample_counts:

  Data frame with columns `<feature names>` and `nsamps`, where
  `<feature names>` are all of the feature columns in the input to
  `AverageAndRegularize()`, and `nsamps` is the number of samples that
  samples from that feature combination needs to have over the read
  count threshold.

- scale_factor_df:

  Data frame with columns "sample" and a second column of whatever name
  you please. The second column should denote scale factors by which
  read counts in that sample should be multiplied by in order to
  normalize these read counts.

- overwrite:

  If TRUE, identical, existing output will be overwritten.

## Value

`EZbakRData` object with an additional "averages" table, as well as a
fullfit table under the same list heading, which includes extra
information about the priors used for regularization purposes.

## Details

The EZbakR website has an extensive vignette walking through various use
cases and model types that you can fit with `AverageAndRegularize()`:
[vignette
link](https://isaacvock.github.io/EZbakR/articles/Linear-modeling.html).
EZbakR improves upon bakR by balancing extra conservativeness in several
steps with a more highly powered statistical testing scheme in its
[`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md)
function. In particular, the following changes to the variance
regularization scheme were made:

- Sample-specific parameter uncertainties are used to generate
  conservative estimates of feature-specific replicate variabilties. In
  addition, a small floor is set to ensure that replicate variance
  estimates are never below a certain level, for the same reason.

- Condition-wide replicate variabilities are regressed against both read
  coverage and either a) \|logit(estimate)\| when modeling average
  fraction labeled. This captures the fact thta estimates are best
  around a logit(fraction labeled) of 0 and get worse for more extreme
  fraction labeled's.; b) log(kdeg) when modeling log degradation rate
  constants. At first, I considered a strategy similar to the fraction
  labeled modeling, but found that agreement between a fully rigorous
  MCMC sampling approach and EZbakR was significantly improved by just
  regressing hee value of the log kientic parameter, likely due to the
  non-linear transformation of fraction labeled to log(kdeg); and c)
  only coverage in all other cases.

- Features with replicate variabilities below the inferred trend have
  their replicate variabilites set equal to that predicted by the trend.
  This helps limit underestimation of parameter variance. Features with
  above-trend replicate variabilties have their replicate variabilities
  regularized with a Normal prior Normal likelihood Bayesian model, as
  in bakR (so the log(variance) is the inferred mean of this
  distribution, and the known variance is inferred from the amount of
  variance about the linear dataset-wide trend).

All of this allows
[`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md)
to use a less conservative statistical test when calculating p-values,
while still controlling false discovery rates.

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
```
