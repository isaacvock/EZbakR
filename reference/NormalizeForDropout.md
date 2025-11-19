# Normalize for experimental/bioinformatic dropout of labeled RNA.

Uses the strategy described
[here](https://simonlabcode.github.io/bakR/articles/Dropout.html), and
similar to that originally presented in [Berg et al.
2024](https://academic.oup.com/nar/article/52/7/e35/7612100), to
normalize for dropout. Normalizing for dropout means identifying a
reference sample with low dropout, and estimating dropout in each sample
relative to that sample.

## Usage

``` r
NormalizeForDropout(
  obj,
  normalize_across_tls = FALSE,
  grouping_factors = NULL,
  features = NULL,
  populations = NULL,
  fraction_design = NULL,
  repeatID = NULL,
  exactMatch = TRUE,
  read_cutoff = 25
)
```

## Arguments

- obj:

  An EZbakRFractions object, which is an EZbakRData object on which you
  have run
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).

- normalize_across_tls:

  If TRUE, samples from different label times will be normalized by
  finding a max inferred degradation rate constant (kdeg) sample and
  using that as a reference. Degradation kinetics with this max will be
  assumed to infer reference fraction news at different label times

- grouping_factors:

  Which sample-detail columns in the metadf should be used to group -s4U
  samples by for calculating the average -s4U RPM? The default value of
  `NULL` will cause no sample-detail columns to be used.

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

- read_cutoff:

  Minimum number of reads for a feature to be used to fit dropout model.

## Value

An `EZbakRData` object with the specified "fractions" table replaced
with a dropout corrected table.

## Details

`NormalizeForDropout()` has a number of unique advantages relative to
[`CorrectDropout()`](https://isaacvock.github.io/EZbakR/reference/CorrectDropout.md):

- `NormalizeForDropout()` doesn't require -label control data.

- `NormalizeForDropout()` compares an internally normalized quantity
  (fraction new) across samples, which has some advantages over the
  absolute dropout estimates derived from comparisons of normalized read
  counts in
  [`CorrectDropout()`](https://isaacvock.github.io/EZbakR/reference/CorrectDropout.md).

- `NormalizeForDropout()` may be used to normalize half-life estimates
  across very different biological contexts (e.g., different cell
  types).

There are also some caveats to be aware of when using
`NormalizeForDropout()`:

- Be careful using `NormalizeForDropout()` when you have multiple
  different label times. Dropout normalization requires each sample be
  compared to a reference sample with the same label time. Thus,
  normalization will be performed separately for groups of samples with
  different label times. If the extent of dropout in the references with
  different label times is different, there will still be unaccounted
  for dropout biases between some of the samples.

- `NormalizeForDropout()` effectively assumes that there are no true
  global differences in turnover kinetics of RNA. If such differences
  actually exist (e.g., half-lives in one context are on average truly
  lower than those in another), `NormalizeForDropout()` risks
  normalizing away these real differences. This is similar to how
  statistical normalization strategies implemented in differential
  expression analysis software like DESeq2 assumes that there are no
  global differences in RNA levels.

By default, all samples with same label time are normalized with respect
to a reference sample chosen from among them. If you want to further
separate the groups of samples that are normalized together, specify the
columns of your metadf by which you want to additionally group factors
in the `grouping_factors` parameter. This behavior can be changed by
setting `normalize_across_tls` to `TRUE`, which will

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

# Normalize for dropout
ezbdo <- NormalizeForDropout(ezbdo)
```
