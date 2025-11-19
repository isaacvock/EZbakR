# `EZbakRData` object helper function for users

`EZbakRData` creates an object of class `EZbakRData` and checks the
validity of the provided input.

## Usage

``` r
EZbakRData(cB, metadf)
```

## Arguments

- cB:

  Data frame with the following columns:

  - sample: Name given to particular sample from which data was
    collected.

  - mutational counts: Integers corresponding to the number of a
    particular mutation seen in a sequencing read. The following column
    names are allowed:

    - TC: Number of Thymine-to-Cytosine mutations

    - TA: Number of Thymine-to-Adenine mutations

    - TG: Number of Thymine-to-Guanine mutations

    - CT: Number of Cytosine-to-Thymine mutations

    - CA: Number of Cytosine-to-Adenine mutations

    - CG: Number of Cytosine-to-Guanine mutations

    - CU: Number of Cytosine-to-Uridine mutations

    - AT: Number of Adenine-to-Thymine mutations

    - AC: Number of Adenine-to-Cytosine mutations

    - AG: Number of Adenine-to-Guanine mutations

    - AU: Number of Adenine-to-Uridine mutations

    - GT: Number of Guanine-to-Thymine mutations

    - GC: Number of Guanine-to-Cytosine mutations

    - GA: Number of Guanine-to-Adenine mutations

    - GU: Number of Guanine-to-Uridine mutations

    - TN: Number of Thymine-to-Adenine/Cytosine/Guanine mutations

    - CN: Number of Cytosine-to-Adenine/Thymine/Guanine/Uridine
      mutations

    - AN: Number of Adenine-to-Thymine/Cytosine/Guanine/Uridine
      mutations

    - GN: Number of Guanine-to-Adenine/Cytosine/Thymine/Uridine
      mutations

    - UN: Number of Uridine-to-Adenine/Cytosine/Guanine mutations

    - NT: Number of Adenine/Cytosine/Guanine-to-Thymine mutations

    - NC: Number of Adenine/Thymine/Guanine/Uridine-to-Cytosine
      mutations

    - NtoA: Number of Thymine/Cytosine/Guanine/Uridine-to-Adenine
      mutations. (Naming convention changed because NA taken)

    - NU: Number of Cytosine/Guanine/Adenine-to-Uridine mutations.

    - NN: Number of any kind of mutation

  - base nucleotide count: Integers corresponding to the number of
    instances of a particular type of nucleotide whose mutations are
    tracked in a corresponding mutation count column. The following
    column names are allowed:

    - nT: Number of Thymines

    - nG: Number of Guanines

    - nA: Number of Adenines

    - nC: Number of Cytosines

    - nU: Number of Uridines

    - nN: number of any kind of nucleotide

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

- metadf:

  Data frame detailing various aspects of each of the samples included
  in the cB. This includes:

  - `sample`: The sample ID, which should correspond to a sample ID in
    the provided cB.

  - `tl`: Metabolic label time. There are several edge cases to be aware
    of:

    - If more than one metabolic label was used in the set of samples
      described by the metadf (e.g., s4U and s6G were used), then the
      `tl` column should be replaced by `tl_<muttype>`, where
      `<muttype>` represents the corresponding mutation type count
      column in the cB that the label whose incubation time will be
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

## Value

An EZbakRData object. This is simply a list of the provide `cB` and
`metadf` with class `EZbakRData`

## Examples

``` r
# Simulate data
simdata <- EZSimulate(30)

# Create EZbakRData object
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
```
