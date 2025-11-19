# `EZbakRKinetics` helper function for users

`EZbakRKinetics` creates an object of class `EZbakRKinetics` and checks
the validity of the provided input.

## Usage

``` r
EZbakRKinetics(kinetics, metadf, features, name = NULL, character_limit = 20)
```

## Arguments

- kinetics:

  Data frame with the following columns:

  - sample: Name given to particular sample from which data was
    collected.

  - features: Any columns that cannot be interpreted as a mutation count
    or base nucleotide count (and that aren't named `sample` or `n`)
    will be interpreted as an ID for a genomic "feature" from which a
    read originated. Common examples of features and typical column
    names for said features include:

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

  - kinetic parameter estimates: These can be named whatever you would
    like as long as they do not start with the string "se\_". This
    should be reserved for kinetic parameter uncertainties, if provided.

  - kinetic parameter uncertainties: Uncertainty in your kinetic
    parameter estimates. These should be named "se\_" followed by the
    kinetic parameter as its name appears in the relevant column name of
    the `kinetics` table.

- metadf:

  Data frame detailing various aspects of each of the samples included
  in the kinetics data frame. This includes:

  - `sample`: The sample ID, which should correspond to a sample ID in
    the provided kinetics data frame.

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
      times respectively. The same \_ convention should be used in the
      case of multi-label pulse-chase designs.

  - sample characteristics: The remaining columns can be named whatever
    you like and should include distinguishing features of groups of
    samples. Common columns might include:

    - `treatment`: The experimental treatment applied to a set of
      samples. This could represent things like genetic knockouts or
      knockdowns, drug treatments, etc.

    - `batch`: An ID for sets of samples that were collected and/or
      processed together. Useful for regressing out technical batch
      effects

  - `assay`: This optional column should include a string that describes
    the type of experiment that was done so as to influence how EZbakR
    analyzes and interprets the data from those samples. Possible values
    for `assay` currently include:

    - standard: Refers to the "standard" nucleotide recoding RNA-seq
      methods (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.), in which
      cells are fed with a single metabolic label, RNA is extracted and
      sequenced, and mutations of a particular type are counted

    - STL: Refers to Start-TimeLapse-seq, a method combining Start-seq
      (developed by Karen Adelman's lab) with TimeLapse-seq. Used to
      infer the kinetics of transcription initiation and
      promoter-proximal pause site departure.

    - TT: Refers to Transient-Transcriptome NR-seq, a method combining
      TT-seq (developed by Patrick Cramer's lab) with NR-seq. TT-seq
      involves biochemically enriching for labeled RNA. By combining
      this method with nucleotide recoding chemistry (as was first done
      by the Simon lab with TT-TimeLapse-seq and has since been done
      with SLAM chemistry, often referred to as TTchem-seq), it is
      possible to bioinformatically filter out reads coming from
      unlabeled RNA background.

    - TILAC: Refers to TILAC, a method developed by the Simon lab to
      achieve spike-in free normalization of RNA-seq data through the
      use of a dual labeling approach inspired by the proteomic method
      SILAC.

    - subcellular: Refers to techniques such as subcellular
      TimeLapse-seq (developed by Stirling Churchman's lab) which
      combine subcellular fractionation with NR-seq to infer additional
      kinetic parameters.

    - sc: Refers to single-cell RNA-seq implementations of NR-seq.

- features:

  Features tracked in `kinetics` data frame. Needs to be specified
  explicitly as it cannot be automatically inferred.

- name:

  Optional; name to give to fractions table.

- character_limit:

  If name is chosen automatically, limit on the number of characters in
  said name. If default name yields a string longer than this, then
  kinetics table will be named `kinetics1`

## Value

An EZbakRKinetics object. This is simply a list of the provided
`kinetics` and `metadf` with class `EZbakRKinetics`

## Examples

``` r
# Simulate data
simdata <- EZSimulate(30)

# Get kinetics table by estimating (for demonstration)
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output
ezbdo <- EstimateKinetics(ezbdo)
kinetics <- EZget(ezbdo, type = "kinetics")

# Create EZbakRKinetics object
ezbko <- EZbakRKinetics(kinetics, simdata$metadf, features = "feature")
#> Kinetic parameter estimates missing corresponding uncertainty
#>           quantification will have an imputed uncertainty of 0.
```
