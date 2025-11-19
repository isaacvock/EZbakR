# Simulate NR-seq data for multiple replicates of multiple biological conditions

`EZSimulate()` is a user friendly wrapper to
[`SimulateMultiCondition()`](https://isaacvock.github.io/EZbakR/reference/SimulateMultiCondition.md).
It sets convenient defaults so as to quickly generate easy to interpret
output. `EZSimulate()` has all of the same parameters as
[`SimulateMultiCondition()`](https://isaacvock.github.io/EZbakR/reference/SimulateMultiCondition.md),
but it also has a number of additional parameters that guide its default
behavior and allow you to simulate multi-condition data without
specifying the multiple, sometimes complex, arguments that you would
need to specify in
[`SimulateMultiCondition()`](https://isaacvock.github.io/EZbakR/reference/SimulateMultiCondition.md)
to get the same behavior. In particular, users only have to set a single
parameter, `nfeatures` (number of features to simulate data for), by
default. The `EZSimulate()`-unique parameters `ntreatments` and `nreps`
have default values that guide the simulation in the case where only
`nfeatures` is specified. In particular, `nreps` of `ntreatments`
different conditions will be simulated, with the assumed model
`log(kdeg) ~ treatment` and `log(ksyn) ~ 1`. In other words, Different
kdeg values will be simulated for each treatment level, and ksyn values
will not differ across conditions.

## Usage

``` r
EZSimulate(
  nfeatures,
  mode = c("standard", "dynamics"),
  ntreatments = ifelse(mode == "standard", 2, 1),
  nreps = 3,
  nctlreps = 1,
  metadf = NULL,
  mean_formula = NULL,
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
  pdo = 0,
  dynamics_preset = c("preRNA", "nuc2cyto", "preRNAwithPdeg", "nuc2cytowithNdeg",
    "subtlseq", "nuc2cytowithpreRNA"),
  unassigned_name = "__no_feature",
  dispersion = 1000,
  lfn_sd = 0.2,
  treatment_effects = NULL,
  effect_avg_default = 0,
  effect_sd_default = 0.5,
  fraction_affected_default = 0.5,
  log_means = NULL,
  log_sds = NULL
)
```

## Arguments

- nfeatures:

  Number of "features" (e.g., genes) for which to simulated data.

- mode:

  Currently, EZSimulate can simulate in two modes: "standard" and
  "dynamics". The former is the default and involves simulating multiple
  conditions of standard NR-seq data. "dynamics" calls
  [`SimulateDynamics()`](https://isaacvock.github.io/EZbakR/reference/SimulateDynamics.md)
  under the hood to simulate a dynamical systems model of your choice.
  Most of the additional parameters do not apply if mode == "dynamics",
  except for those from dynamics_preset and on.

- ntreatments:

  Number of distinct treatments to simulate. This parameter is only
  relevant if `metadf` is not provided.

- nreps:

  Number of replicates of each treatment to simulate. This parameter is
  only relevant if `metadf` is not provided

- nctlreps:

  Number of -s4U replicates of each treatment to simulate. This
  parameter is only relevant if `metadf` is not provided.

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

  Total number of reads in each sample.

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

- dynamics_preset:

  Which preset model to use for simulation of dynamics. Therefore, only
  relevant if `mode` == `dynamics`. Options are:

  nuc2cyto

  :   Simplest model of nuclear and cytoplasmic RNA dynamics: 0 -\> N
      -\> C -\> 0

  preRNA

  :   Simplest model of pre-RNA and mature RNA dynamics: 0 -\> P -\> M
      -\> 0

  preRNAwithPdeg

  :   Same as preRNA, but now pre-RNA can also degrade.

  nuc2cytowithNdeg

  :   Same as nuc2cyto, but now nuclear RNA can also degrade.

  subtlseq

  :   Subcellular TimeLapse-seq model, similar to that described in
      Ietswaart et al., 2024. Simplest model discussed there, lacking
      nuclear degradation: 0 -\> CH -\> NP -\> CY -\> PL -\> 0, and CY
      can also degrade.

  nuc2cytowithpreRNA

  :   Combination of nuc2cyto and preRNA where preRNA is first
      synthesized, then either processed or exported to the cytoplasm.
      Processing can also occur in the cytoplasm, and mature nuclear RNA
      can be exported to the cytoplasm. Only mature RNA degrades.

- unassigned_name:

  String to give to reads not assigned to a given feature.

- dispersion:

  Negative binomial `size` parameter to use for simulating read counts

- lfn_sd:

  Logit(fn) replicate variability.

- treatment_effects:

  Data frame describing effects of treatment on each parameter. Should
  have five columns: "parameter_index", "treatment_index", "mean", "sd",
  and "fraction_affected". Each row corresponds to the effect the ith (i
  = treatment_index) treatment has on the jth (j = parameter_index)
  kinetic parameter. Effect sizes, on a log-scale, are drawn from a
  Normal distribution with mean and standard deviation set by the mean
  and sd columns, respectively. The number of non-zero effects is set by
  "fraction_affected", and is equal to
  `ceiling(nfeatures * fraction_affected)`. treatment_index of 1 will be
  ignored and can either be included or not.

- effect_avg_default:

  If `ntreatments` \> 1, and `treatment_effects` is not provided, this
  will be the value of `mean` for all treatments and parameters imputed
  in `treatment_effects`.

- effect_sd_default:

  If `ntreatments` \> 1, and `treatment_effects` is not provided, this
  will be the value of `sd` for all treatments and parameters imputed in
  `treatment_effects`.

- fraction_affected_default:

  If `ntreatments` \> 1, and `treatment_effects` is not provided, this
  will be the value of `fraction_affected` for all treatments and
  parameters imputed in `treatment_effects`.

- log_means:

  Vector of log-Normal logmeans from which the distribution of
  feature-specific parameters will be drawn from. Length of vector
  should be the same as max(entries in `graph`), i.e., the number of
  parameters in your specified model. If not provided, will by default
  be `c(1, seq(from = -0.3, to = -2.5, length.out = max(graph) - 1 ))`.
  `1` for the ksyn parameter (which is always denoted 1 in the preset
  `graph`) is arbitrary. Remaining parameters will make it so indices
  order parameters from fastest to slowest process.

- log_sds:

  Vector of log-Normal logsds from which the distribution of
  feature-specific parameters will be drawn from. If not provided, will
  be 0.4 for all parameters.

## Examples

``` r
# Simulate standard data
simdata_standard <- EZSimulate(30)

# Simulate dynamical systems data
simdata_ode <- EZSimulate(30,
                          mode = "dynamics",
                          ntreatments = 1,
                          label_time = c(1, 3),
                          dynamics_preset = "nuc2cyto")
```
