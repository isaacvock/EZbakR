# Simulate NR-seq data for multiple replicates of multiple biological conditions

`SimulateMultiCondition` is a highly flexibly simulator that combines
linear modeling of log(kdeg)'s and log(ksyn)'s with `SimulateOneRep` to
simulate an NR-seq dataset. The linear model allows you to simulate
multiple distinct treatments, batch effects, interaction effects, etc.
The current downside for its flexibility is its relative complexity to
implement. Easier to use simulators are on the way to EZbakR.

## Usage

``` r
SimulateMultiCondition(
  nfeatures,
  metadf,
  mean_formula,
  param_details = NULL,
  seqdepth = nfeatures * 2500,
  label_time = 2,
  pnew = 0.05,
  pold = 0.001,
  readlength = 200,
  Ucont_alpha = 25,
  Ucont_beta = 75,
  feature_prefix = "Gene",
  dispslope = 5,
  dispint = 0.01,
  logkdegsdtrend_slope = -0.3,
  logkdegsdtrend_intercept = -2.25,
  logksynsdtrend_slope = -0.3,
  logksynsdtrend_intercept = -2.25,
  logkdeg_mean = -1.9,
  logkdeg_sd = 0.7,
  logksyn_mean = 2.3,
  logksyn_sd = 0.7,
  logkdeg_diff_avg = 0,
  logksyn_diff_avg = 0,
  logkdeg_diff_sd = 0.5,
  logksyn_diff_sd = 0.5,
  pdiff_kd = 0.1,
  pdiff_ks = 0,
  pdiff_both = 0,
  pdo = 0
)
```

## Arguments

- nfeatures:

  Number of "features" (e.g., genes) to simulate data for

- metadf:

  A data frame with the following columns:

  - sample: Names given to samples to simulate.

  - : Any number of columns with any names (not taken by other metadf
    columns) storing factors by which the samples can be stratified.
    These can be referenced in `mean_formula`, described below.

  These parameters (described more below) can also be included in metadf
  to specify sample-specific simulation parameter:

  - seqdepth

  - label_time

  - pnew

  - pold

  - readlength

  - Ucont

- mean_formula:

  A formula object that specifies the linear model used to relate the
  factors in the

  columns of `metadf` to average log(kdegs) and log(ksyns) in each
  sample.

- param_details:

  A data frame with one row for each column of the design matrix
  obtained from `model.matrix(mean_formula, metadf)` that describes how
  to simulate the linear model parameters. The columns of this data
  frame are:

  - param: Name of linear model parameter as it appears in the column
    names of the design matrix from
    `model.matrix(mean_formula, metadf)`.

  - reference: Boolean; TRUE if you want to treat that parameter as a
    "reference". This means that all other parameter values that aren't
    global parameters are set equal to this unless otherwise determined
    (see `pdiff_*` parameters for how it is determined if a parameter
    will differ from the reference).

  - global: Boolean; TRUE if you want to treat that parameter as a
    global parameter. This means that a single value is used for all
    features.

  - logkdeg_mean: If parameter is the reference, then its value for the
    log(kdeg) linear model will be drawn from a normal distribution with
    this mean. If it is a global parameter, then this value will be
    used. If it is neither of these, then its value in the log(kdeg)
    linear model will either be the reference (if there is no difference
    between this condition's value and the reference) or the reference's
    value + a normally distributed random variable centered on this
    value.

  - logkdeg_sd: sd used for draws from normal distribution as described
    for `logkdeg_mean`.

  - logksyn_mean: Same as `logkdeg_mean` but for log(ksyn) linear model.

  - logksyn_sd: Same as `logkdeg_sd` but for log(kdeg) linear model.

  - pdiff_ks: Proportion of features whose value of this parameter in
    the log(ksyn) linear model will differ from the reference's. Should
    be a number between 0 and 1, inclusive. For example, if `pdiff_ks`
    is 0.1, then for 10% of features, this parameter will equal the
    reference parameter + a normally distributed random variable with
    mean `logksyn_mean` and sd `logksyn_sd`. For the other 90% of
    features, this parameter will equal the reference.

  - pdiff_kd: Same as `pdiff_ks` but for log(kdeg) linear model.

  - pdiff_both: Proportion of features whose value for this parameter in
    BOTH the log(kdeg) and log(ksyn) linear models will differ from the
    reference. Value must be between 0 and min(c(pdiff_kd, pdiff_ks)) in
    that row.

  If param_details is not specified by the user, the first column of the
  design matrix is assumed to represent the reference parameter, all
  parameters are assumed to be non-global, logkdeg_mean and logksyn_mean
  are set to the equivalently named parameter values described below for
  the reference and `logkdeg_diff_avg` and `logksyn_diff_avg` for all
  other parameters, logkdeg_sd and logksyn_sd are set to the
  equivalently named parameter values described below for the reference
  and `logkdeg_diff_sd` and `logksyn_diff_sd` for all other parameters,
  and pdiff_kd, pdiff_ks, and pdiff_both are all set to the equivalently
  named parameter values.

- seqdepth:

  Only relevant if `read_vect` is not provided; in that case, this is
  the total number of reads to simulate.

- label_time:

  Length of s^4^U feed to simulate.

- pnew:

  Probability that a T is mutated to a C if a read is new.

- pold:

  Probability that a T is mutated to a C if a read is old.

- readlength:

  Length of simulated reads. In this simple simulation, all reads are
  simulated as being exactly this length.

- Ucont_alpha:

  Probability that a nucleotide in a simulated read from a given feature
  is a U is drawn from a beta distribution with shape1 = `Ucont_alpha`.

- Ucont_beta:

  Probability that a nucleotide in a simulated read from a given feature
  is a U is drawn from a beta distribution with shape2 = `Ucont_beta`.

- feature_prefix:

  Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows
  up in the `feature` column of the output simulated data table.

- dispslope:

  Negative binomial dispersion parameter "slope" with respect to read
  counts. See DESeq2 paper for dispersion model used.

- dispint:

  Negative binomial dispersion parameter "intercept" with respect to
  read counts. See DESeq2 paper for dispersion model used.

- logkdegsdtrend_slope:

  Slope for log10(read count) vs. log(kdeg) replicate variability trend

- logkdegsdtrend_intercept:

  Intercept for log10(read count) vs. log(kdeg) replicate variability
  trend

- logksynsdtrend_slope:

  Slope for log10(read count) vs. log(ksyn) replicate variability trend

- logksynsdtrend_intercept:

  Intercept for log10(read count) vs. log(ksyn) replicate variability
  trend

- logkdeg_mean:

  Mean of normal distribution from which reference log(kdeg) linear
  model parameter is drawn from for each feature if `param_details` is
  not provided.

- logkdeg_sd:

  Standard deviation of normal distribution from which reference
  log(kdeg) linear model parameter is drawn from for each feature if
  `param_details` is not provided.

- logksyn_mean:

  Mean of normal distribution from which reference log(ksyn) linear
  model parameter is drawn from for each feature if `param_details` is
  not provided.

- logksyn_sd:

  Standard deviation of normal distribution from which reference
  log(ksyn) linear model parameter is drawn from for each feature if
  `param_details` is not provided.

- logkdeg_diff_avg:

  Mean of normal distribution from which non-reference log(kdeg) linear
  model parameters are drawn from for each feature if `param_details` is
  not provided.

- logksyn_diff_avg:

  Mean of normal distribution from which reference log(ksyn) linear
  model parameter are drawn from for each feature if `param_details` is
  not provided.

- logkdeg_diff_sd:

  Standard deviation of normal distribution from which reference
  log(kdeg) linear model parameter are drawn from for each feature if
  `param_details` is not provided.

- logksyn_diff_sd:

  Standard deviation of normal distribution from which reference
  log(ksyn) linear model parameter are drawn from for each feature if
  `param_details` is not provided.

- pdiff_kd:

  Proportion of features for which non-reference log(kdeg) linear model
  parameters differ from the reference.

- pdiff_ks:

  Proportion of features for which non-reference log(ksyn) linear model
  parameters differ from the reference.

- pdiff_both:

  Proportion of features for which BOTH non-reference log(kdeg) and
  log(ksyn) linear model parameters differ from the reference. ksyns are
  simulated

- pdo:

  Dropout rate; think of this as the probability that a s4U containing
  molecule is lost during library preparation and sequencing. If `pdo`
  is 0 (default) then there is not dropout.

## Examples

``` r
simdata <- SimulateMultiCondition(30,
                                  data.frame(sample = c('sampleA', 'sampleB'),
                                  treatment = c('treatment1', 'treatment2')),
                                  mean_formula = ~treatment-1)
```
