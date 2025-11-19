# Simulation of generalized dynamical system model.

`SimulateDynamics()` simulates any specified dynamical system of
interconverting RNA species. Its required input is similar to that of
`EstimateDynamics()`, i.e., an adjacency matrix describing the set of
species and how they are related to one another and a list of formula
relating actually assayed species to the modeled species. Currently,
`SimulateDynamics()` implements a naive heteroskedastic replicate
variability simulation and is not designed to simulate multiple
experimental conditions.

## Usage

``` r
SimulateDynamics(
  nfeatures,
  graph,
  metadf,
  log_means,
  log_sds,
  ntreatments = 1,
  treatment_effects = NULL,
  formula_list = NULL,
  unassigned_name = "__no_feature",
  seqdepth = nfeatures * 2500,
  dispersion = 100,
  lfn_sd = 0.2,
  effect_avg_default = 0,
  effect_sd_default = 0.5,
  fraction_affected_default = 0.5,
  ...
)
```

## Arguments

- nfeatures:

  Number of "features" to simulate data for. A "feature" in this case
  may contain a number of "sub-features". For example, you may want to
  simulate pre-RNA and mature RNA for a set of "genes", in which case
  the number of features is the number of genes.

- graph:

  An adjacency matrix describing the reaction diagram graph relating the
  various RNA species to one another.

- metadf:

  Data frame with two required columns (`sample` and `tl`). `sample`
  represents names given to each simulated sample. `tl` represents the
  label time for that sample. Additional columns can specify other
  features of the sample, like what subcellular compartment the sample
  is taken from. **NOTE: Not sure I am actually using these optional
  columns in any useful capacity anymore**.

- log_means:

  Vector of log-Normal logmeans from which the distribution of
  feature-specific parameters will be drawn from. Length of vector
  should be the same as max(entries in `graph`), i.e., the number of
  parameters in your specified model.

- log_sds:

  Vector of log-Normal logsds from which the distribution of
  feature-specific parameters will be drawn from.

- ntreatments:

  Number of distinct experimental treatments to simulate. By default,
  only a single "treatment" (you might refer to this as wild-type, or
  control) is simulated. Increase this if you would like to explore
  performing comparative dynamical systems modeling.

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

- formula_list:

  A list of named lists. The names of each sub-list should be the same
  as the sample names as they are found in `metadf`. Each sub-list
  should be a list of formula relating feature names that will show up
  as columns of the simulated cB to species modeled in your `graph`.
  This only needs to be specified if you want to simulate the scenario
  where some of the measured species are a sum of modeled species.

- unassigned_name:

  String to give to reads not assigned to a given feature.

- seqdepth:

  Total number of reads in each sample.

- dispersion:

  Negative binomial `size` parameter to use for simulating read counts

- lfn_sd:

  Logit(fn) replicate variability.

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

- ...:

  Parameters passed to
  [`SimulateOneRep()`](https://isaacvock.github.io/EZbakR/reference/SimulateOneRep.md).
