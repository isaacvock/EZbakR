# Deconvolve multi-feature fractions.

Combines the output of `EstimateFractions` with feature quantification
performed by an outside tool (e.g., RSEM, kallisto, salmon, etc.) to
infer fraction news for features that reads cannot always be assigned to
unambiguously. This is a generalization of `EstimateIsoformFractions`
which performs this deconvolution for transcript isoforms.

## Usage

``` r
DeconvolveFractions(
  obj,
  feature_type = c("gene", "isoform"),
  features = NULL,
  populations = NULL,
  fraction_design = NULL,
  repeatID = NULL,
  exactMatch = TRUE,
  fraction_name = NULL,
  quant_name = NULL,
  gene_to_transcript = NULL,
  overwrite = TRUE,
  TPM_min = 1,
  count_min = 10
)
```

## Arguments

- obj:

  An `EZbakRData` object

- feature_type:

  Either `"gene"` (if deconvolving gene-level fraction news, i.e., for
  resolving fusion gene and component gene fraction news) or `"isoform"`
  (if deconvolving transcript isoform fraction news).

- features:

  Character vector of the set of features you want to stratify reads by
  and estimate proportions of each RNA population. The default of "all"
  will use all feature columns in the `obj`'s cB.

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

  If TRUE, then `features` and `populations` have to exactly match those
  for a given fractions table for that table to be used. Means that you
  can't specify a subset of features or populations by default, since
  this is TRUE by default.

- fraction_name:

  Name of fraction estimate table to use. Should be stored in the
  `obj$fractions` list under this name. Can also rely on specifying
  `features` and/or `populations` and having
  [`EZget()`](https://isaacvock.github.io/EZbakR/reference/EZget.md)
  find it.

- quant_name:

  Name of transcript isoform quantification table to use. Should be
  stored in the obj\$readcounts list under this name. Use
  [`ImportIsoformQuant()`](https://isaacvock.github.io/EZbakR/reference/ImportIsoformQuant.md)
  to create this table. If `quant_name` is `NULL`, it will search for
  tables containing the string "isoform_quant" in their name, as that is
  the naming convention used by
  [`ImportIsoformQuant()`](https://isaacvock.github.io/EZbakR/reference/ImportIsoformQuant.md).
  If more than one such table exists, an error will be thrown and you
  will have to specify the exact name in `quant_name`.

- gene_to_transcript:

  Table with columns `transcript_id` and all feature related columns
  that appear in the relevant fractions table. This is only relevant as
  a hack to to deal with the case where STAR includes in its
  transcriptome alignment transcripts on the opposite strand from where
  the RNA actually originated. This table will be used to filter out
  such transcript-feature combinations that should not exist.

- overwrite:

  If TRUE and a fractions estimate output already exists that would
  possess the same metadata (features analyzed, populations analyzed,
  and fraction_design), then it will get overwritten with the new
  output. Else, it will be saved as a separate output with the same
  name + "\_#" where "#" is a numerical ID to distinguish the similar
  outputs.

- TPM_min:

  Minimum TPM for a transcript to be kept in analysis.

- count_min:

  Minimum expected_count for a transcript to be kept in analysis.

## Value

An `EZbakRData` object with an additional table under the "fractions"
list. Has the same form as the output of
[`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md),
and will have the feature column "transcript_id".

## Details

`DeconvolveFractions` expects as input a "fractions" table with
estimates for fraction news of at least one convolved feature set. A
convolved feature set is one where some reads cannot be unambiguously
assigned to one instance of that feature type. For example, it is often
impossible to assign short reads to a single transcript isoform. Thus,
something like the "TEC" assignment provided by fastq2EZbakR is an
instance of a convolved feature set, as it is an assignment of reads to
transcript isoforms with which they are compatible. Another example is
assignment to the exonic regions of genes, for fusion genes (where a
read may be consistent with both the putative fusion gene as well as one
of the fusion components).

`DeconvolveFractions` deconvolves fraction news by fitting a linear
mixing model to the convolved fraction new estimates + feature abundance
estimates. In other words, each convolved fraction new (data) is modeled
as a weighted average of single feature fraction news (parameters to
estimate), with the weights determined by the relative abundances of the
features in the convolved set (data). The convolved fraction new is
modeled as a beta distribution with mean equal to the weighted feature
fraction new average and variance related to the number of reads in the
convolved feature set.

Features with estimated TPMs less than `TPM_min` (1 by default) or an
estimated number of read counts less than `count_min` (10 by default)
are removed from convolved feature sets and will not have their fraction
news estimated.

## Examples

``` r
# Load dependencies
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# Simulates a single sample worth of data
simdata_iso <- SimulateIsoforms(nfeatures = 30)

# We have to manually create the metadf in this case
metadf <- tibble(sample = 'sampleA',
                     tl = 4,
                     condition = 'A')

ezbdo <- EZbakRData(simdata_iso$cB,
                    metadf)

ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output

### Hack in the true, simulated isoform levels
reads <- simdata_iso$ground_truth %>%
  dplyr::select(transcript_id, true_count, true_TPM) %>%
  dplyr::mutate(sample = 'sampleA',
                effective_length = 10000) %>%
  dplyr::rename(expected_count = true_count,
                TPM = true_TPM)

# Name of table needs to have "isoform_quant" in it
ezbdo[['readcounts']][['simulated_isoform_quant']] <- reads

### Perform deconvolution
ezbdo <- DeconvolveFractions(ezbdo, feature_type = "isoform")
#> Analyzing sample sampleA...
```
