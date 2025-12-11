# Run quality control checks

`EZQC()` assesses multiple aspects of your NR-seq data and generates a
number of plots visualizing dataset-wide trends.

## Usage

``` r
EZQC(obj, ...)
```

## Arguments

- obj:

  EZbakRData or EZbakRFractions object.

- ...:

  Parameters passed to the class-specific method. If you have provided
  an EZbakRFractions object, then these can be (all play the same role
  as in
  [`EstimateKinetics()`](https://isaacvock.github.io/EZbakR/reference/EstimateKinetics.md),
  that is they get passed to
  [`EZget()`](https://isaacvock.github.io/EZbakR/reference/EZget.md) to
  find the fractions table you are interested in. See
  `?EstimateKinetics()` for details.):

  - `features`

  - `populations`

  - `fraction_design`

  If you have provided an EZbakRData object, then these can be (all same
  the same purpose as in `EstimateFractions`, so see
  `?EstimateFractions()` for details):

  - `mutrate_populations`

  - `features`

  - `filter_cols`

  - `filter_condition`

  - `remove_features`

## Value

A list of `ggplot2` objects visualizing the various aspects of your data
assessed by `EZQC()`.

## Details

`EZQC()` checks the following aspects of your NR-seq data. If you have
passed an `EZbakRData` object, then `EZQC()` checks:

- Raw mutation rates: In all sequencing reads, how many T's in the
  reference were a C in the read? The hope is that raw mutation rates
  are higher than -label controls in all +label samples. Higher raw
  mutation rates, especially when using standard label times (e.g., 2
  hours or more in mammalian systems), are typically a sign of good
  label incorporation and low labeled RNA/read dropout. If you don't
  have -label samples, know that background mutation rates are typically
  less than 0.2%, so +label raw mutation rates several times higher than
  this would be preferable.

- Mutation rates in labeled and unlabeled reads: The raw mutation rate
  counts all mutations in all reads. In a standard NR-seq experiment
  performed with a single metabolic label, there are typically two
  populations of reads:

  1.  Those from labeled RNA, having higher mutation rates due to
      chemical conversion/recoding of the metabolic label and 2) those
      from unlabeled RNA, having lower, background levels of mutations.
      `EZbakR` fits a two component mixture model to estimate the
      mutation rates in these two populations separately. A successful
      NR-seq experiment should have a labeled read mutation rate of \>
      1% and a low background mutation rate of \< 0.3%.

  2.  Read count replicate correlation: Simply the log10 read count
      correlation for replicates, as inferred from your `metadf`.

If you have passed an `EZbakRFractions` object, i.e., the output of
[`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md),
then in addition to the checks in the `EZbakRData` input case, `EZQC()`
also checks:

- Fraction labeled distribution: This is the distribution of
  feature-wise fraction labeled's (or fraction high mutation content's)
  estimated by
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
  The "ideal" is a distribution with mean around 0.5, as this maximizes
  the amount of RNA with synthesis and degradation kinetics within the
  dynamic range of the experiment. In practice, you will (and should) be
  at least a bit lower than this as longer label times risk
  physiological impacts of metabolic labeling.

- Fraction labeled replicate correlation: This is the logit(fraction
  labeled) correlation between replicates, as inferred from your
  `metadf`.

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

# Run QC
QC <- EZQC(ezbdo)
#> CHECKING RAW MUTATION RATES...
#> CHECKING INFERRED MUTATION RATES...
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the EZbakR package.
#>   Please report the issue at <https://github.com/isaacvock/EZbakR/issues>.
#> CHECKING READ COUNT CORRELATIONS...
#> log10(read counts) correlation for each pair of replicates are:
#> # A tibble: 12 × 3
#>    sample_1 sample_2 correlation
#>    <chr>    <chr>          <dbl>
#>  1 sample1  sample2        0.980
#>  2 sample1  sample3        0.979
#>  3 sample1  sample7        0.988
#>  4 sample2  sample3        0.975
#>  5 sample2  sample7        0.974
#>  6 sample3  sample7        0.974
#>  7 sample4  sample5        0.977
#>  8 sample4  sample6        0.986
#>  9 sample4  sample8        0.981
#> 10 sample5  sample6        0.978
#> 11 sample5  sample8        0.975
#> 12 sample6  sample8        0.983
#> 
#> log10(read counts) correlations are high, suggesting good reproducibility!
#> 
#> CHECKING FRACTION LABELED DISTRIBUTIONS...
#> Average fractions for each sample are:
#> # A tibble: 6 × 3
#>   sample  avg_fraction fraction_type  
#>   <chr>          <dbl> <chr>          
#> 1 sample1        0.263 fraction_highTC
#> 2 sample2        0.269 fraction_highTC
#> 3 sample3        0.266 fraction_highTC
#> 4 sample4        0.273 fraction_highTC
#> 5 sample5        0.270 fraction_highTC
#> 6 sample6        0.269 fraction_highTC
#> 
#> Labeling rates (e.g., fraction labeled for single label experiments) look good!
#> 
#> CHECKING FRACTION LABELED CORRELATIONS...
#> logit(fraction_highTC) correlation for each pair of replicates are:
#> # A tibble: 6 × 3
#>   sample_1 sample_2 correlation
#>   <chr>    <chr>          <dbl>
#> 1 sample1  sample2        0.961
#> 2 sample1  sample3        0.975
#> 3 sample2  sample3        0.955
#> 4 sample4  sample5        0.935
#> 5 sample4  sample6        0.944
#> 6 sample5  sample6        0.962
#> 
#> logit(fraction_highTC) correlations are high, suggesting good reproducibility!
#> 
```
