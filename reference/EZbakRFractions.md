# `EZbakRFractions` helper function for users

`EZbakRFractions` creates an object of class `EZbakRFractions` and
checks the validity of the provided input.

## Usage

``` r
EZbakRFractions(fractions, metadf, name = NULL, character_limit = 20)
```

## Arguments

- fractions:

  Data frame with the following columns:

  - sample: Name given to particular sample from which data was
    collected.

  - estimates of population fractions: These columns refer to the
    estimate for the fraction of reads coming from a particular
    mutational population. For example, in a standard NR-seq experiment,
    you should have one column named `fraction_highTC`. This refers to
    the fraction of RNA inferred to have a high T-to-C mutation rate
    (e.g., the newly synthesized RNA in a pulse-labeling NR-seq
    experiment). If you estimated the fractions of more than 2 mutation
    types (e.g., T-to-C and G-to-A), then you need to explicitly list
    all fractions of interest estimated. For example, in a TILAC
    experiment, this would be `fraction_highTC_lowGA`,
    `fraction_lowTC_highGA`, and `fraction_lowTC_lowGA`.

  - n: Number of reads assigned to a given feature in a given sample.

  - features: Any columns that cannot be interpreted as an estimate of
    population fractions (and that aren't named `sample` or `n`) will be
    interpreted as an ID for a genomic "feature" from which a read
    originated. Common examples of features and typical column names for
    said features include:

    - Genes; common column names: gene, gene_id, gene_name, GF

    - Genes-exonic; common column names: gene_exon, gene_id_exon,
      gene_name_exon, XF

    - Transcripts; common column names: transcripts, TF

    - Exonic bins; common column names: exonic_bins, EF, EB

    - Exons; common column names: exons, exon_ids

    In some cases, a read will often map to multiple features (e.g.,
    exons). Many functions in bakR expect each of the feature IDs in
    these cases to be separated by `+`. For example, if a read overlaps
    with two exons, with IDs exon_1 and exon_2, then the corresponding
    entry in a column of exonic assignments would be "exon_1+exon_2".
    The default expectation can be overwritten though and is thus not
    strictly enforced.

  - n: Number of reads with identical values for all other columns.

- metadf:

  Data frame detailing various aspects of each of the samples included
  in the fractions data frame. This includes:

  - `sample`: The sample ID, which should correspond to a sample ID in
    the provided fractions data frame.

  - `tl`: Metabolic label time. There are several edge cases to be aware
    of:

    - If more than one metabolic label was used in the set of samples
      described by the metadf (e.g., s4U and s6G were used), then the
      `tl` column should be replaced by `tl_<muttype>`, where
      `<muttype>` represents the corresponding mutation type referenced
      in the fractions that the label whose incubation time will be
      listed in this column. For example, if feeding with s4U in some
      samples and s6G in others, then performing standard nucleotide
      recoding chemistry, you will include `tl_TC` and `tl_GA` columns
      corresponding to the s4U and s6G label times, respectively.

    - If a pulse-chase experimental design was used (!!this is strongly
      discouraged unless you have a legitimate reason to prefer this
      design to a pulse-label design!!), then you should have columns
      named `tpulse` and `tchase`, corresponding to the pulse and chase
      times respectively. The same `_<muttype>` convention should be
      used in the case of multi-label pulse-chase designs.

  - sample characteristics: The remaining columns can be named whatever
    you like and should include distinguishing features of groups of
    samples. Common columns might include:

    - `treatment`: The experimental treatment applied to a set of
      samples. This could represent things like genetic knockouts or
      knockdowns, drug treatments, etc.

    - `batch`: An ID for sets of samples that were collected and/or
      processed together. Useful for regressing out technical batch
      effects

- name:

  Optional; name to give to fractions table.

- character_limit:

  Maximum number of characters for naming out fractions output. EZbakR
  will try to name this as a "\_" separated character vector of all of
  the features analyzed. If this name is greater than `character_limit`,
  then it will default to "fraction#", where "#" represents a simple
  numerical ID for the table.

## Value

An EZbakRFractions object. This is simply a list of the provide
`fractions` and `metadf` with class `EZbakRFractions`

## Examples

``` r
# Simulate data
simdata <- EZSimulate(30)

# Get fractions table by estimating (for demonstration)
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output
fxns <- EZget(ezbdo, type = "fractions")

# Create EZbakRFractions object
ezbfo <- EZbakRFractions(fxns, simdata$metadf)
```
