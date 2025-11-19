# Estimate mutation rates

Two component mixture models are fit to all data to estimate global high
and low mutation rates for all samples. Estimation of these mutation
rates are regularized through the use of weakly informative priors whose
parameters can be altered using the arguments defined below.

## Usage

``` r
EstimateMutRates(
  obj,
  populations = "all",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL
)

# S3 method for class 'EZbakRData'
EstimateMutRates(
  obj,
  populations = "all",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL
)

# S3 method for class 'EZbakRArrowData'
EstimateMutRates(
  obj,
  populations = "all",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL
)
```

## Arguments

- obj:

  An `EZbakRData` or `EZbakRArrowData` object

- populations:

  Character vector of the set of mutational populations that you want to
  infer the fractions of. For example, say your cB file contains columns
  tracking T-to-C and G-to-A

- pnew_prior_mean:

  logit-Normal mean for logit(pnew) prior.

- pnew_prior_sd:

  logit-Normal sd for logit(pnew) prior.

- pold_prior_mean:

  logit-Normal mean for logit(pold) prior.

- pold_prior_sd:

  logit-Normal sd for logit(pold) prior.

- pold_est:

  Background mutation rate estimates if you have them. Can either be a
  single number applied to all samples or a named vector of values,
  where the names should be sample names.

- pold_from_nolabel:

  Fix background mutation rate estimate to mutation rates seen in -label
  samples. By default, a single background rate is used for all samples,
  inferred from the average mutation rate across all -label samples. The
  `grouping_factors` argument can be specified to use certain -label
  samples to infer background mutation rates for certain sets of +label
  samples.

- grouping_factors:

  If `pold_from_nolabel` is TRUE, then `grouping_factors` will specify
  the sample-detail columns in the metadf that should be used to group
  -label samples by. Average mutation rates in each group of -label
  samples will be used as the background mutation rate estimate in
  +label samples with the same values for the relevant metadf columns.

## Value

`EZbakRData` object with an added `mutrates` slot containing estimated
high and low mutation rates for each mutation type modeled.

## Methods (by class)

- `EstimateMutRates(EZbakRData)`: Method for class **EZbakRData**
  Estimates mutation rates using a fully in memory object.

- `EstimateMutRates(EZbakRArrowData)`: Method for class
  **EZbakRArrowData** Estimate mutation rates using a partially on-disk
  object.

## Examples

``` r
# Simulate data to analyze
simdata <- SimulateOneRep(30)

# Create EZbakR input
metadf <- data.frame(sample = "sampleA", tl = 2)
ezbdo <- EZbakRData(simdata$cB, metadf)

# Estimate mutation rates
mutrates <- EstimateMutRates(ezbdo)
```
