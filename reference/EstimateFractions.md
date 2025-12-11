# Estimate fractions of each RNA population.

The first step of any NR-seq analysis is to figure out the fraction of
reads from each mutational population in your data. For example, if you
are performing a standard SLAM-seq or TimeLapse-seq experiment, this
means estimating the fraction of reads with high T-to-C mutation
content, and the fraction with low T-to-C mutation content. This is what
`EstimateFractions` is for.

## Usage

``` r
EstimateFractions(
  obj,
  features = "all",
  mutrate_populations = "all",
  fraction_design = NULL,
  Poisson = TRUE,
  strategy = c("standard", "hierarchical"),
  filter_cols = "all",
  filter_condition = `&`,
  remove_features = c("NA", "__no_feature"),
  split_multi_features = FALSE,
  multi_feature_cols = NULL,
  multi_feature_sep = "+",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  hier_readcutoff = 300,
  init_pnew_prior_sd = 0.8,
  pnew_prior_sd_min = 0.01,
  pnew_prior_sd_max = 0.15,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL,
  character_limit = 20,
  overwrite = TRUE
)

# S3 method for class 'EZbakRData'
EstimateFractions(
  obj,
  features = "all",
  mutrate_populations = "all",
  fraction_design = NULL,
  Poisson = TRUE,
  strategy = c("standard", "hierarchical"),
  filter_cols = "all",
  filter_condition = `&`,
  remove_features = c("NA", "__no_feature"),
  split_multi_features = FALSE,
  multi_feature_cols = NULL,
  multi_feature_sep = "+",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  hier_readcutoff = 300,
  init_pnew_prior_sd = 0.3,
  pnew_prior_sd_min = 0.01,
  pnew_prior_sd_max = 0.15,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL,
  character_limit = 20,
  overwrite = TRUE
)

# S3 method for class 'EZbakRArrowData'
EstimateFractions(
  obj,
  features = "all",
  mutrate_populations = "all",
  fraction_design = NULL,
  Poisson = TRUE,
  strategy = c("standard", "hierarchical"),
  filter_cols = "all",
  filter_condition = `&`,
  remove_features = c("NA", "__no_feature"),
  split_multi_features = FALSE,
  multi_feature_cols = NULL,
  multi_feature_sep = "+",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  hier_readcutoff = 300,
  init_pnew_prior_sd = 0.8,
  pnew_prior_sd_min = 0.01,
  pnew_prior_sd_max = 0.15,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL,
  character_limit = 20,
  overwrite = TRUE
)
```

## Arguments

- obj:

  `EZbakRData` or `EZbakRArrowData` object

- features:

  Character vector of the set of features you want to stratify reads by
  and estimate proportions of each RNA population. The default of "all"
  will use all feature columns in the `obj`'s cB file.

- mutrate_populations:

  Character vector of the set of mutational populations that you want to
  infer the rates of mutations for. By default, all mutation rates are
  estimated for all populations present in cB.

- fraction_design:

  "Design matrix" specifying which RNA populations exist in your
  samples. By default, this will be created automatically and will
  assume that all combinations of the `mutrate_populations` you have
  requested to analyze are present in your data. If this is not the case
  for your data, then you will have to create one manually.

  If you call the function `create_fraction_design(...)`, providing a
  vector of mutational population names as input, it will create a
  `fraction_design` table for you, with the assumption that every single
  possible combination of mutational populations is present in your
  data. You can then edit the `present` column as necessary to get an
  appropriate `fraction_design` for your use case. See below for details
  on the required contents of `fraction_design` and its interpretation.

  `fraction_design` must have one column per element of
  `mutrate_populations`, with these columns sharing the name of the
  `mutrate_populations`. It must also have one additional column named
  `present`. All elements of fraction_design should be booleans (`TRUE`
  or `FALSE`). It should include all possible combinations of `TRUE` and
  `FALSE` for the `mutrate_populations` columns. A `TRUE` in one of
  these columns represents a population of RNA that is expected to have
  above background mutation rates of that type. `present` will denote
  whether or not that population of RNA is expected to exist in your
  data.

  For example, assume you are doing a typical
  TimeLapse-seq/SLAM-seq/TUC-seq/etc. experiment where you have fed
  cells with s^4U and recoded any incorporated s^4U to a nucleotide that
  reverse transcriptase will read as a cytosine. That means that
  `mutrate_populations` will be "TC", since you want to estimate the
  fraction of RNA that was s^4U labeled, i.e., the fraction with high
  T-to-C mutation content. `fraction_design` will thus have two columns:
  `TC` and `present`. It will also have two rows. One of these rows must
  have a value of `TRUE` for `TC`, and the other must have a value of
  `FALSE`. The row with a value of `TRUE` for `TC` represents the
  population of reads with high T-to-C mutation content, i.e., the reads
  from RNA that were synthesized while s^4U was present. The row with a
  value of `FALSE` for `TC` reprsents the population of reads with low
  T-to-C mutation content, i.e., the reads from RNA that existed prior
  to s^4U labeling. Both of these populations exist in your data, so the
  value of the `present` column should be `TRUE` for both of these. See
  the lazily loaded `standard_fraction_design` object for an example of
  what this tibble could look like. ("lazily loaded
  `standard_fraction_design` object" means that if you run
  `print(standard_fraction_design)` after loading `EZbakR` with
  [`library(EZbakR)`](https://isaacvock.github.io/EZbakR/), then you can
  see its contents. More specifically, lazily loaded means that this
  table is not loaded into memory until you ask for it, via something
  like a [`print()`](https://rdrr.io/r/base/print.html) call.)

  As another example, consider TILAC, a NR-seq extension developed by
  the Simon lab. TILAC was originally described in [Courvan et al.,
  2022](https://pubmed.ncbi.nlm.nih.gov/36018791/). In this method, two
  populations of RNA, one from s^4U fed cells and one from s^6G fed
  cells, are pooled and prepped for sequencing together. This allows for
  internally controlled comparisons of RNA abundance without spike-ins.
  s^4U is recoded to a cytosine analog by TimeLapse chemistry (or
  similar chemistry) and s^6G is recoded to an adenine analog. Thus,
  `fraction_design` includes columns called `TC` and `GA`. A unique
  aspect of the TILAC `fraction_design` table is that one of the
  possible populations, `TC` and `GA` both `TRUE`, is denoted as not
  present (`present` = `FALSE`). This is because there is no RNA that
  was exposed to both s^4U and s^6G, thus a population of reads with
  both high T-to-C and G-to-A mutational content should not exist. To
  see an example of what a TILAC `fraction_design` table could look
  like, see the lazily loaded `tilac_fraction_design` object.

- Poisson:

  If `TRUE`, use U-content adjusted Poisson mixture modeling strategy.
  Often provides significant speed gain without sacrificing accuracy.

- strategy:

  String denoting which new read mutation rate estimation strategy to
  use. Options include:

  - standard: Estimate a single new read and old read mutation rate for
    each sample. This is done via a binomial mixture model aggregating
    over

  - hierarchical: Estimate feature-specific new read mutation rate,
    regularizing the feature-specific estimate with a sample-wide prior.
    Currently only compatible with single mutation type mixture
    modeling.

- filter_cols:

  Which feature columns should be used to filter out feature-less reads.
  The default value of "all" checks all feature columns for whether or
  not a read failed to get assigned to said feature.

- filter_condition:

  Only two possible values for this make sense: `` `&` `` and `` `|` ``.
  If set to `` `&` ``, then all features in `filter_cols` must have a
  "null" value (i.e., a value included in `remove_features`) for the row
  to get filtered out. If set to `` `|` ``, then only a single feature
  in `filter_cols` needs to have one of these "null" values to get
  filtered out.

- remove_features:

  All of the feature names that could indicate failed assignment of a
  read to a given feature. the fastq2EZbakR pipeline uses a value of
  '\_\_no_feature'.

- split_multi_features:

  If a set of reads maps ambiguously to multiple features, should data
  for such reads be copied for each feature in the ambiguous set? If
  this is `TRUE`, then `multi_feature_cols` also must be set. Examples
  where this should be set to `TRUE` includes when analyzing exonic bins
  (concept defined in original DEXSeq paper), exon-exon junctions, etc.

- multi_feature_cols:

  Character vector of columns that have the potential to include
  assignment to multiple features. Only these columns will have their
  features split if `split_multi_features` is `TRUE`.

- multi_feature_sep:

  String representing how ambiguous feature assignments are
  distinguished in the feature names. For example, the default value of
  "+" denotes that if a read maps to multiple features (call them
  featureA and featureB, for example), then the feature column will have
  a value of "featureA+featureB" for that read.

- pnew_prior_mean:

  Mean for logit(pnew) prior.

- pnew_prior_sd:

  Standard deviation for logit(pnew) prior.

- pold_prior_mean:

  Mean for logit(pold) prior.

- pold_prior_sd:

  Standard deviation for logit(pold) prior.

- hier_readcutoff:

  If `strategy` == `hierarchical`, only features with this many reads
  are used to infer the distribution of feature-specific labeled read
  mutation rates.

- init_pnew_prior_sd:

  If `strategy` == `hierarchical`, this is the initial logit(pnew) prior
  standard deviation to regularize feature-specific labeled read
  mutation rate estimates.

- pnew_prior_sd_min:

  The minimum logit(pnew) prior standard deviation when `strategy` is
  set to "hierarchcial". EZbakR will try to estimate this empirically as
  the standard deviation of initial feature-specific logit(pnew)
  estimates using high coverage features, minus the average uncertainty
  in the logit(pnew) estimates. As this difference can sometimes be
  negative, a value of `pnew_prior_sd_min` will be imputed in that case.

- pnew_prior_sd_max:

  Similar to `pnew_prior_sd_min`, but now representing the maximum
  allowed logit(pnew) prior sd.

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

- character_limit:

  Maximum number of characters for naming out fractions output. EZbakR
  will try to name this as a "\_" separated character vector of all of
  the features analyzed. If this name is greater than `character_limit`,
  then it will default to "fraction#", where "#" represents a simple
  numerical ID for the table.

- overwrite:

  If TRUE and a fractions estimate output already exists that would
  possess the same metadata (features analyzed, populations analyzed,
  and fraction_design), then it will get overwritten with the new
  output. Else, it will be saved as a separate output with the same
  name + "\_#" where "#" is a numerical ID to distinguish the similar
  outputs.

## Value

An `EZbakRFractions` object, which is just an `EZbakRData` object with a
"fractions" list element. This will include tables of fraction estimates
(e.g., fraction of reads from the high T-to-C mutation rate population
in a standard single-label s4U experiment; termed fraction_highTC in the
table) for all features in all samples.

## Details

`EstimateFractions` uses mixture modeling to estimate the fraction of
reads from each mutational population in your data, and this is done for
each feature in your data (i.e., combination of columns that specify
genomic features from which reads were derived). The set of mutational
populations in your data can be specified by providing a
`fraction_design` object, described in more depth above (also see
[`?create_fraction_design`](https://isaacvock.github.io/EZbakR/reference/create_fraction_design.md)).
There are several flavors of mixture modeling that can be performed by
`EstimateFractions`. These are as follows:

1.  The default: global mutation rate parameters (global = same for all
    reads from all features) are estimated for each sample by fitting a
    single two-component mixture model to all of the reads in that
    sample. These are used to estimate the fraction of reads from each
    feature that are from each mutational population, also using a
    two-component mixture model.

    - With `Poisson` set to `TRUE`, this is a nucleotide-content
      adjusted Poisson mixture model, which is a more efficient
      alternative to binomial mixture modeling.

    - With `Poisson` set to `FALSE`, this is a binomial mixture model.

2.  Low mutation rate from -label: if `pold_from_nolabel` is `TRUE`,
    then the background, no label mutation rates are estimated from
    -label samples. By default, a single set of background mutation
    rates are estimated for all samples, but you can change this
    behavior by setting `grouping_factors` to specify columns in your
    `metadf` by which samples should be stratified.

    - This is a great strategy to use if you have low label
      incorporation rates or if you used a fairly short label time.

3.  Hierarchical model: if `strategy == "hierarchical"`, which is
    currently only compatible for single-mutation type modeling (e.g.,
    standard T-to-C mutation modeling), then high T-to-C content
    mutation rates are estimated for each feature. The global
    sample-wide estimates are used as an informative prior to increase
    the accuracy of this process by avoiding extreme estimates.

## Methods (by class)

- `EstimateFractions(EZbakRData)`: Method for class **EZbakRData**
  Estimates fractions using an entirely in memory object.

- `EstimateFractions(EZbakRArrowData)`: Mehthod for class
  **EZbakRArrowData** This is an alternative to the fully in memory
  **EZbakRData** `EstimateFractions` method that can help with analyses
  of larger than RAM datasets. The provided "cB" is expected to be an
  on-disk Arrow Dataset. Furthermore, it is expected to be partitioned
  by the sample name, which will allow this method to read only a
  single-sample worth of data into memory at a time. This can
  significantly reduce RAM usage. Input object should be created with
  [`EZbakRArrowData()`](https://isaacvock.github.io/EZbakR/reference/EZbakRArrowData.md).

## Examples

``` r
# Simulate data to analyze
simdata <- SimulateOneRep(30)

# Create EZbakR input
metadf <- data.frame(sample = "sampleA", tl = 2)
ezbdo <- EZbakRData(simdata$cB, metadf)

# Estimate fractions
ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output
```
