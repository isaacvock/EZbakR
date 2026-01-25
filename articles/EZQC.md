# QCing your NR-seq data with EZbakR: EZQC()

## Introduction

Quality control checks are an important part of any high throughput
dataset analysis, and analyzing NR-seq data is no different. Therefore,
EZbakR provides a function,
[`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.md) to help
you identify potential problems in your NR-seq data. In this section, we
will discuss how to run
[`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.md) and
what it looks for. I will also provide some suggestions for how to best
design and analyze NR-seq experiments. Let’s start by loading EZbakR,
which we will use throughout this vignette:

``` r
library(EZbakR)
```

**NOTE**:
[`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.md) is
EZbakR’s instantiation of bakR’s `QC_checks()`, though differs in a
number of key ways. Thus, don’t expect its output or behavior to exactly
mimic that of `QC_checks()`. That being said, much of the discussion of
QCing NR-seq data present in [bakR’s
vignette](https://simonlabcode.github.io/bakR/articles/Troubleshooting.html)
is still relevant for interpreting EZbakR’s output.

## EZQC

[`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.md) can
take two different inputs:

1.  An `EZbakRData` object. This can be created with
    [`EZbakRData()`](https://isaacvock.github.io/EZbakR/reference/EZbakRData.md).
2.  An `EZbakRFractions` object. This is an `EZbakRData` object on which
    you have run
    [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).

An example of both of these are shown below:

``` r
simdata <- EZSimulate(250)

ezbdo <- EZbakRData(simdata$cB,
                    simdata$metadf)

### Input: EZbakRData object
qc <- EZQC(ezbdo)
#> CHECKING RAW MUTATION RATES...
#> CHECKING INFERRED MUTATION RATES...
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the EZbakR package.
#>   Please report the issue at <https://github.com/isaacvock/EZbakR/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> CHECKING READ COUNT CORRELATIONS...
#> log10(read counts) correlation for each pair of replicates are:
#> # A tibble: 12 × 3
#>    sample_1 sample_2 correlation
#>    <chr>    <chr>          <dbl>
#>  1 sample1  sample2        0.984
#>  2 sample1  sample3        0.983
#>  3 sample1  sample7        0.983
#>  4 sample2  sample3        0.985
#>  5 sample2  sample7        0.986
#>  6 sample3  sample7        0.985
#>  7 sample4  sample5        0.982
#>  8 sample4  sample6        0.984
#>  9 sample4  sample8        0.986
#> 10 sample5  sample6        0.983
#> 11 sample5  sample8        0.985
#> 12 sample6  sample8        0.986
#> 
#> log10(read counts) correlations are high, suggesting good reproducibility!
#> 


### Input: EZbakRFractions object
ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output
qc_fn <- EZQC(ezbdo)
#> CHECKING RAW MUTATION RATES...
#> CHECKING INFERRED MUTATION RATES...
#> CHECKING READ COUNT CORRELATIONS...
#> log10(read counts) correlation for each pair of replicates are:
#> # A tibble: 12 × 3
#>    sample_1 sample_2 correlation
#>    <chr>    <chr>          <dbl>
#>  1 sample1  sample2        0.984
#>  2 sample1  sample3        0.983
#>  3 sample1  sample7        0.983
#>  4 sample2  sample3        0.985
#>  5 sample2  sample7        0.986
#>  6 sample3  sample7        0.985
#>  7 sample4  sample5        0.982
#>  8 sample4  sample6        0.984
#>  9 sample4  sample8        0.986
#> 10 sample5  sample6        0.983
#> 11 sample5  sample8        0.985
#> 12 sample6  sample8        0.986
#> 
#> log10(read counts) correlations are high, suggesting good reproducibility!
#> 
#> CHECKING FRACTION LABELED DISTRIBUTIONS...
#> Average fractions for each sample are:
#> # A tibble: 6 × 3
#>   sample  avg_fraction fraction_type  
#>   <chr>          <dbl> <chr>          
#> 1 sample1        0.303 fraction_highTC
#> 2 sample2        0.307 fraction_highTC
#> 3 sample3        0.303 fraction_highTC
#> 4 sample4        0.308 fraction_highTC
#> 5 sample5        0.309 fraction_highTC
#> 6 sample6        0.307 fraction_highTC
#> 
#> Labeling rates (e.g., fraction labeled for single label experiments) look good!
#> 
#> CHECKING FRACTION LABELED CORRELATIONS...
#> logit(fraction_highTC) correlation for each pair of replicates are:
#> # A tibble: 6 × 3
#>   sample_1 sample_2 correlation
#>   <chr>    <chr>          <dbl>
#> 1 sample1  sample2        0.964
#> 2 sample1  sample3        0.957
#> 3 sample2  sample3        0.969
#> 4 sample4  sample5        0.949
#> 5 sample4  sample6        0.969
#> 6 sample5  sample6        0.962
#> 
#> logit(fraction_highTC) correlations are high, suggesting good reproducibility!
#> 
```

In the first case (`EZbakRData` input), the following are checked:

1.  Correlation of read counts between replicates
2.  Raw mutation rates
3.  Labeled and unlabeled read mutation rates

In the second case (`EZbakRFractions` input), all of the same things are
checked, with the addition of:

1.  Correlation of fraction estimates (e.g., fraction labeled in single
    label experiments) between replicates.
2.  Distribution of fraction estimates.

Replicate correlation is an intuitive QC metric that ensures replicates
agree with each other well. The other three metrics are NR-seq specific
metrics that assess the extent to which metabolic label was readily
incorporated into nascent RNA, and the appropriateness of the metabolic
label feed time used (i.e., how long cells were fed with the label,
sometimes referred to as the label time). In general, you want to see:

1.  High labeled read mutation rates (referred to as pnew within
    EZbakR).
2.  Low unlabeled read mutation rates (referred to as pold within
    EZbakR).
3.  Roughly equal proportions of all mutational populations (e.g.,
    labeled and unlabeled reads in a single label experiment).

### EZQC output

Inside of the output objects (`qc` and `qc_fn` in the above code), you
will find a number of plots. The output is a named list, with one or
more plot present in each element of this list. The named elements and
their contents are as follows:

1.  `Raw_mutrates`: If you have a single mutation type in your cB, this
    will be a barplot of raw mutation rates of that type in all samples.
    If you have multiple mutation types, then this will be a named list,
    with one element per mutation type, with each element being this
    same barplot.
    - The raw mutation rate = (total \# of mutations) / (total \# of
      mutable nucleotides)
    - Raw mutation rate should be high in samples fed the relevant
      label, and low in unfed samples.
2.  `Inferred_mutrates`: If you have a single mutation type in your cB,
    this will be a barplot of the inferred labeled and unlabeled read
    mutation rates in all samples. If you have multiple mutation types,
    then it will be a named list of similar barplots for each mutation
    type, as in `Raw_mutrates`.
    - The inferred labeled mutation rate is the rate of a mutation of a
      particular type seen in reads from RNA synthesized in the presence
      of the relevant metabolic label. The inferred unlabeled mutation
      rate is this same mutation rate, but in reads from RNA synthesized
      prior to labeling. These are sometimes referred to as the
      conversion and background rates, respectively.
    - The plot includes a subtitle with some guidelines for
      interpretation. Inferred labeled mutation rates is preferably much
      higher than the unlabeled mutation rates in all samples.
3.  `Readcount_corr`: This is a list, whose elements are collections of
    read count correlation plots for a given group of replicates. Groups
    of replicates are determined by the provided metadf in your
    `EZbakRData` object. +label and -label replicates of a given
    treatment condition are grouped together. Different label times for
    the same treatment condition are also grouped together here.
    - The plots are named according to the samples compared (e.g.,
      “sample1_vs_sample2”).
4.  `Fraction_labeled_dist`: This is a list of density plots for each
    label fed sample. In the case of a single label experiment, these
    will be the distribution of estimated fraction labeled’s for each
    sample. If you have multiple labels/mutation types, then this will
    be a list with density plots for each estimated fraction.
    - Ideally, these will peak around 1 / (# of populations). For a
      standard, single label experiment, this would be 0.5 (# of
      populations = 2; labeled and unlabeled).
    - This is only present if you provided a `EZbakRFractions` object as
      input!
5.  `Fraction_labeled_corr`: Same as `Readcount_corr`, but this time
    correlating fraction estimates provided by
    [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
    Unlike in `Readcount_corr`, this excludes -label samples, and
    considers distinct label times distinct replicates.
    - This is only present if you provided a `EZbakRFractions` object as
      input!

### What to do if the QC looks iffy?

**Potential problem \#1**: Labeled read mutation rates are low

Possible solutions:

1.  Better read alignment settings to ensure recovery of high mutation
    content reads.
    [fastq2EZbakR’s](https://github.com/isaacvock/fastq2EZbakR) default
    [config](https://github.com/isaacvock/fastq2EZbakR/blob/main/config/config.yaml)
    includes STAR settings used by multiple labs to help with this.
    Tools like [grandRescue](https://github.com/erhard-lab/grandRescue)
    developed by the Erhard lab may also help.
2.  If you ran fastq2EZbakR, did you properly set the strandedness of
    your library? You can check some of the intermediate files produced
    by fastq2EZbakR to test this hypothesis; see
    [here](https://fastq2ezbakr.readthedocs.io/en/latest/faq/) for
    details.
3.  If using a short label time, or if incorporation rates are low,
    EZbakR may just be having a tough time estimating the labeled read
    mutation rate. If you have -label data, you can use this to infer a
    single background mutation rate for all samples, and then estimate
    the labeled read mutation rate, fixing this unlabeled read mutation
    rate. This can be done by setting `pold_from_no_label=TRUE` in your
    call to
    [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
    In fact, I even like defaulting to this setting when -label samples
    are available.
4.  You may need to optimize the concentration of metabolic label being
    fed to cells. A simple TAMRA-labeling dot blot assay can be used to
    broadly access the incorporation of metabolic label in your
    particular system.

**Potential problem \#2**: Labeled read mutation rates are acceptable,
but raw mutation rates are low. In this case, fraction estimate
replicate correlation may also be low.

Possible solutions:

1.  You may need to use a longer label time to ensure better
    representation of labeled RNA in your samples.
2.  Either way, using `pold_from_no_label=TRUE` in your
    [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md)
    call, as described above, is likely a good idea to ensure accurate
    labeled read mutation rate inference.

**Potential problem \#3**: Poor correlation of read counts across
replicates, especially between -label and +label, or different label
time +label, samples.

This may be a sign of dropout, described
[here](https://pubmed.ncbi.nlm.nih.gov/38381903/), ,
[here](https://pubmed.ncbi.nlm.nih.gov/37292657/), and
[here](https://simonlabcode.github.io/bakR/articles/Dropout.html). This
is when reads from highly labeled RNA are underrepresented in +label
data. You can try:

1.  Decreasing metabolic label concentration.
2.  Decreasing metabolic label feed time.
3.  Using the RNA handling protocol suggested
    [here](https://pubmed.ncbi.nlm.nih.gov/37292657/).
4.  Better read alignment settings to avoid loss of high mutation
    content reads. See the potential solutions to problem \#1 described
    above for details.
5.  Running EZbakR’s
    [`CorrectDropout()`](https://isaacvock.github.io/EZbakR/reference/CorrectDropout.md)
    to try and bioinformatically correct for the biases induced by
    dropout. This requires you to have -label data from all distinct
    biological conditions.

Suggestions 1 and 2 can help no matter the cause of dropout (e.g.,
adverse effects of the label, RNA handling problems, read alignment
problems, etc.). The other three suggestions are more suited in the case
when dropout does not seem to have a biological origin. A general,
proven strategy to distinguish the cause of dropout does not exist, but
you may want to try assessing trends in sequencing tracks (e.g., look
dropoff in coverage near the 3’ end of transcripts upon increased
labeling) or performing differential expression and GO analysis of + and
-label samples to assess potential upregulation of transcriptional
repression stress pathways.

## Characteristics of a good NR-seq experiment:

Below I will list and discuss what, in my opinion, makes an ideal NR-seq
dataset:

1.  Multiple +label (e.g., s4U) replicates for each biological condition
    being assessed. Nothing radical here, replicates are important in
    all settings.
2.  At least one, though preferably multiple, replicates of -label
    (e.g., regular RNA-seq) controls in each biological condition being
    assessed. This is a crucial control to ensure that the label(s) used
    do not have a significant impact on the treated cells. These can
    also help correct for certain instances of bias (i.e.,
    [dropout](https://simonlabcode.github.io/bakR/articles/Dropout.html))
    present in some NR-seq datasets.
3.  A pulse-label (label for certain amount of time and extract RNA)
    rather than a pulse-chase (label for certain amount of time, then
    chase with regular nucleotide for a certain amount of time, then
    extract RNA) design. I’ve ranted about this
    [elsewhere](https://simonlabcode.github.io/bakR/articles/Troubleshooting.html)
    (bottom of linked vignette), but in short, a pulse-label design is
    almost always optimal relative to a pulse-chase design.
    Pulse-chases:
    - suffer from longer exposure of cells to metabolic label, thus
      incurring a greater risk for physiological impact.
    - require comparisons across samples to infer kinetic parameters,
      yielding noisier estimates.
    - Are almost never necessary as the dynamics of the unlabeled RNA in
      a pulse-label experiment are identical to that of the labeled RNA
      in a pulse-chase experiment.
    - Rare exceptions to these rules exist, but unless you have a really
      good reason to use a pulse-chase design, pulse-labels should be
      your default.
4.  A metabolic label feed time length around that of the average
    half-life of the RNAs of interest. Shorter is better than longer to
    limit any adverse effects of the metabolic label.
    - In humans and mice, this means about 4 hours, if there isn’t a
      specific species of RNA you are interested in probing. 2 hour
      label times are typical and tend to work great in these systems.
    - In yeast this is much shorter, around 30 minutes.
5.  An optimized concentration of metabolic label.
    - Typical values range from 100 $\mu$M to 1 mM, but simple, low
      throughput assays (e.g., dot blots) should be used to ensure
      sufficient incorporation of metabolic label at a given
      concentration. Different cell lines can exhibit wildly different
      levels of label uptake, and thus it would be best to assess this
      in each unique biological context you plan to do NR-seq in.
    - Lower concentrations are better to avoid adverse effects, but make
      sure that incorporation at lower concentrations is still robust
      before committing.
