# Make plots to visually assess dropout trends

Plots a measure of dropout (the ratio of -label to +label RPM) as a
function of feature fraction new, with the model fit depicted. Use this
function to qualitatively assess model fit and whether the modeling
assumptions are met.

## Usage

``` r
VisualizeDropout(
  obj,
  plot_type = c("grandR", "bakR"),
  grouping_factors = NULL,
  features = NULL,
  populations = NULL,
  fraction_design = NULL,
  repeatID = NULL,
  exactMatch = TRUE,
  n_min = 50,
  dropout_cutoff = 5,
  ...
)
```

## Arguments

- obj:

  An EZbakRFractions object, which is an EZbakRData object on which you
  have run
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).

- plot_type:

  Which type of plot to make. Options are:

  - bakR: X-axis is fraction new (a.k.a. NTR) and Y-axis is dropout (no
    label n / label n)

  - grandR: X-axis is fraction new rank (a.k.a. NTR rank) and Y-axis is
    log(dropout)

- grouping_factors:

  Which sample-detail columns in the metadf should be used to group -s4U
  samples by for calculating the average -s4U RPM? The default value of
  `NULL` will cause all sample-detail columns to be used.

- features:

  Character vector of the set of features you want to stratify reads by
  and estimate proportions of each RNA population. The default of `NULL`
  will expect there to be only one fractions table in the
  EZbakRFractions object.

- populations:

  Mutational populations that were analyzed to generate the fractions
  table to use. For example, this would be "TC" for a standard s4U-based
  nucleotide recoding experiment.

- fraction_design:

  "Design matrix" specifying which RNA populations exist in your
  samples. By default, this will be created automatically and will
  assume that all combinations of the `mutrate_populations` you have
  requested to analyze are present in your data. If this is not the case
  for your data, then you will have to create one manually. See docs for
  `EstimateFractions` (run ?EstimateFractions()) for more details.

- repeatID:

  If multiple `fractions` tables exist with the same metadata, then this
  is the numerical index by which they are distinguished.

- exactMatch:

  If TRUE, then `features` must exactly match the `features` metadata
  for a given fractions table for it to be used. Means that you cannot
  specify a subset of features by default. Set this to FALSE if you
  would like to specify a feature subset.

- n_min:

  Minimum raw number of reads to make it to plot

- dropout_cutoff:

  Maximum dropout value included in plot.

- ...:

  Parameters passed to internal `calculate_dropout()` function;

## Value

A list of `ggplot2` objects, one for each +label sample.

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

# Visualize Dropout
ezbdo <- VisualizeDropout(ezbdo)
#> Making plot for sample sample1...
#> Making plot for sample sample2...
#> Making plot for sample sample3...
#> Making plot for sample sample4...
#> Making plot for sample sample5...
#> Making plot for sample sample6...
```
