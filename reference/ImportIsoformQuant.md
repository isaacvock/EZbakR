# Import transcript isoform quantification into EZbakRData object

A convenient wrapper to
[`tximport()`](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
for importing isoform quantification data into an EZbakRData object. You
need to run this before running `EstimateIsoformFractions`.

## Usage

``` r
ImportIsoformQuant(
  obj,
  files,
  quant_tool = c("none", "salmon", "sailfish", "alevin", "piscem", "kallisto", "rsem",
    "stringtie"),
  txIn = TRUE,
  ...
)
```

## Arguments

- obj:

  An `EZbakRData` object.

- files:

  A named vector of paths to all transcript quantification files that
  you would like to import. This will be passed as the first argument of
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html)
  (also named `files`). The names of this vector should be the same as
  the sample names as they appear in the metadf of the `EZbakRData`
  object.

- quant_tool:

  String denoting the type of software used to generate the abundances.
  Will get passed to the `type` argument of
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html).
  As described in the documentation for `tximport` 'Options are
  "salmon", "sailfish", "alevin", "piscem", "kallisto", "rsem",
  "stringtie", or "none". This argument is used to autofill the
  arguments below (geneIdCol, etc.) "none" means that the user will
  specify these columns. Be aware that specifying type other than "none"
  will ignore the arguments below (geneIdCol, etc.)'. Referenced
  'arguments below' can be specified as part of `...`.

- txIn:

  Whether or now you are providing isoform level quantification files.
  Alternative (`txIn = FALSE`) is gene-level quantification. In
  `ImportIsoformQuant`, `txIn` gets passed to BOTH the `txIn` and
  `txOut` parameters in `tximport()`.

- ...:

  Additional arguments to be passed to
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html).
  Especially relevant if you set `quant_tool` to "none".

## Value

An `EZbakRData` object with an additional element in the `readcounts`
list named "isform_quant\_\<quant_tool\>". It contains TPM,
expected_count, and effective length information for each transcript_id
and each sample.

## Examples

``` r
# Dependencies for example
library(dplyr)
library(data.table)
#> 
#> Attaching package: ‘data.table’
#> The following objects are masked from ‘package:dplyr’:
#> 
#>     between, first, last

# Simulate and analyze data
simdata <- EZSimulate(30)
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
ezbdo <- EstimateFractions(ezbdo)
#> Estimating mutation rates
#> Summarizing data for feature(s) of interest
#> Averaging out the nucleotide counts for improved efficiency
#> Estimating fractions
#> Processing output

# Hack to generate example quantification files
savedir <- tempdir()
rsem_data <- tibble(
  transcript_id = paste0("tscript_feature", 1:30),
  gene_id = paste0("feature", 1:30),
  length = 1000,
  effective_length = 1000,
  expected_count = 1000,
  TPM = 10,
  FPKM = 10,
  IsoPct = 1
)

fwrite(rsem_data, file.path(savedir, "Sample_1.isoforms.results"), sep = '\t')

files <- file.path(savedir,"Sample_1.isoforms.results")
names(files) <- "Sample_1"

# Read in file
ezbdo <- ImportIsoformQuant(ezbdo, files, quant_tool = "rsem")
#> reading in files with read.delim (install 'readr' package for speed up)
#> 1 
#> 


```
