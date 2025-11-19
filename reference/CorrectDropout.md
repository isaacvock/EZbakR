# Correct for experimental/bioinformatic dropout of labeled RNA.

Uses the strategy described
[here](https://simonlabcode.github.io/bakR/articles/Dropout.html), and
similar to that originally presented in [Berg et al.
2024](https://academic.oup.com/nar/article/52/7/e35/7612100).

## Usage

``` r
CorrectDropout(
  obj,
  strategy = c("grandR", "bakR"),
  grouping_factors = NULL,
  features = NULL,
  populations = NULL,
  fraction_design = NULL,
  repeatID = NULL,
  exactMatch = TRUE,
  read_cutoff = 25,
  dropout_cutoff = 5,
  ...
)
```

## Arguments

- obj:

  An EZbakRFractions object, which is an EZbakRData object on which you
  have run
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).

- strategy:

  Which dropout correction strategy to use. Options are:

  - grandR: Described
    [here](https://academic.oup.com/nar/article/52/7/e35/7612100). Cite
    that work and
    [grandR](https://www.nature.com/articles/s41467-023-39163-4) if
    using this strategy. Quasi-non-parametric strategy that finds an
    estimate of the dropout rate that eliminates any linear correlation
    between the newness of a transcript and the difference in +s4U and
    -s4U normalized read counts.

  - bakR: Described
    [here](https://simonlabcode.github.io/bakR/articles/Dropout.html).
    Uses a simple generative model of dropout to derive a likelihood
    function, and the dropout rate is estimated via the method of
    maximum likelihood.

  The "bakR" strategy has the advantage of being model-derived, making
  it possible to assess model fit and thus whether the simple
  assumptions of both the "bakR" and "grandR" dropout models are met.
  The "grandR" strategy has the advantage of being more robust. Thus,
  the "grandR" strategy is currently used by default.

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

- read_cutoff:

  Minimum number of reads for a feature to be used to fit the dropout
  model.

- dropout_cutoff:

  Maximum ratio of -s4U:+s4U RPMs for a feature to be used to fit the
  dropout model (i.e., simple outlier filtering cutoff).

- ...:

  Parameters passed to internal `calculate_dropout()` function; namely
  `dropout_cutoff_min`, which sets the minimum dropout value used for
  fitting the dropout model.

## Value

An `EZbakRData` object with the specified "fractions" table replaced
with a dropout corrected table.

## Details

Dropout is the disproportionate loss of labeled RNA/reads from said RNA
described independently
[here](https://academic.oup.com/nar/article/52/7/e35/7612100) and
[here](https://www.biorxiv.org/content/10.1101/2023.05.24.542133v1). It
can originate from a combination of bioinformatic (loss of high mutation
content reads due to alignment problems), technical (loss of labeled RNA
during RNA extraction), and biological (transcriptional shutoff in rare
cases caused by metabolic label toxicity) sources. `CorrectDropout()`
compares label-fed and label-free controls from the same experimental
conditions to estimate and correct for this dropout. It assumes that
there is a single number (referred to as the dropout rate, or pdo) which
describes the rate at which labeled RNA is lost (relative to unlabeled
RNA). pdo ranges from 0 (no dropout) to 1 (complete loss of all labeled
RNA), and is thus interpreted as the percentage of labeled RNA/reads
from labeled RNA disproportionately lost, relative to the equivalent
unlabeled species.

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

# Correct for dropout
ezbdo <- CorrectDropout(ezbdo)
#> Estimated rates of dropout are:
#>    sample        pdo
#> 1 sample1 0.15062004
#> 2 sample2 0.01000000
#> 3 sample3 0.09230763
#> 4 sample4 0.01000000
#> 5 sample5 0.03995805
#> 6 sample6 0.01000000
```
