# Easily get EZbakR table of estimates of interest

`EZget()` returns a table of interest from your `EZbakRData` object. It
is meant to make it easier to find and access certain analyses, as a
single `EZbakRData` object may include analyses of different features,
kinetic parameters, dynamical systems models, comparisons, etc.

## Usage

``` r
EZget(
  obj,
  type = c("fractions", "kinetics", "readcounts", "averages", "comparisons", "dynamics"),
  features = NULL,
  populations = NULL,
  fraction_design = NULL,
  isoforms = NULL,
  kstrat = NULL,
  parameter = NULL,
  counttype = NULL,
  design_factor = NULL,
  dynamics_design_factors = NULL,
  scale_factors = NULL,
  cstrat = NULL,
  feature_lengths = NULL,
  experimental = NULL,
  reference = NULL,
  param_name = NULL,
  param_function = NULL,
  formula_mean = NULL,
  sd_grouping_factors = NULL,
  fit_params = NULL,
  repeatID = NULL,
  sub_features = NULL,
  grouping_features = NULL,
  sample_feature = NULL,
  modeled_to_measured = NULL,
  graph = NULL,
  normalize_by_median = NULL,
  deconvolved = NULL,
  returnNameOnly = FALSE,
  exactMatch = FALSE,
  alwaysCheck = FALSE
)
```

## Arguments

- obj:

  EZbakRData object

- type:

  The class of EZbakR outputs would you like to search through.
  Equivalent to the name of the list in the EZbakRData object that
  contains the tables of interest.

- features:

  Features that must be present in the table of interest. If
  `exactMatch` is TRUE, then these features must also be the only
  features present in the table.

- populations:

  Only relevant if `type` == "fractions". Mutational populations that
  must have been analyzed to generate the table of interest.

- fraction_design:

  Only relevant if `type` == "fractions". Fraction design table used to
  generate the table of interest.

- isoforms:

  If the relevant table is the result of isoform deconvolution

- kstrat:

  Only relevant if `type` == "kinetics". Short for "kinetics strategy";
  the strategy used to infer kinetic parameters.

- parameter:

  Only relevant if `type` == "averages" or "comparisons". Which
  parameter was being averaged or compared?

- counttype:

  String denoting what type of read count information you are looking
  for. Current options are "TMM_normalized", "transcript", and "matrix".
  TO-DO: Not sure this is being used in any way currently...

- design_factor:

  design_factor specified in relevant run of
  [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md).
  Therefore, only relevant if type == "comparisons".

- dynamics_design_factors:

  design_factors included in final
  [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md)
  output. Therefore, only relevant if type == "dynamics".

- scale_factors:

  Sample group scale factors used in
  [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md).
  Therefore, only relevant if type == "dynamics"

- cstrat:

  Strategy used for comparative analyses. Can be:

  - contrast: If two parameters were compared via specifying the
    reference and experimental levels in
    [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md)
    (for type == "averages").

  - single_param: If a single parameter was passed to
    [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md)
    via its `param_name` option.

  - dynamics: If output of
    [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md)
    was passed to
    [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md)

  - function: If function of multiple parameters was passed to
    `CompareParameter()` via its `param_function` option.

- feature_lengths:

  Table of feature lengths used for length normalization.

- experimental:

  Experimental condition specified in relevant run of
  [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md).
  Therefore, only relevant if type == "comparisons".

- reference:

  Reference condition specified in relevant run of
  [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md).
  Therefore, only relevant if type == "comparisons".

- param_name:

  Parameter name specified in relevant run of
  [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md).
  Therefore, only relevant if type == "comparisons"

- param_function:

  Function of parameters specified in relevant run of
  [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md).
  Therefore, only relevant if type == "comparisons".

- formula_mean:

  An R formula object specifying how the `parameter` of interest depends
  on the sample characteristics specified in `obj`'s metadf. Therefore,
  only relevant if type == "averages".

- sd_grouping_factors:

  What metadf columns should data be grouped by when estimating standard
  deviations across replicates? Therefore, only relevant if type ==
  "averages".

- fit_params:

  Character vector of names of parameters in linear model fit.
  Therefore, only relevant if type == "averages".

- repeatID:

  Numerical ID for duplicate objects with same metadata.

- sub_features:

  Only relevant if type == "dynamics". Feature columns that
  distinguished between the different measured species when running
  [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md).

- grouping_features:

  Only relevant if type == "dynamics. Features that were the overarching
  feature assignments by which `sub_features` were grouped when running
  [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md).

- sample_feature:

  Only relevant if type == "dynamics". Name of the metadf column that
  distinguished the different classes of samples when running
  [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md).

- modeled_to_measured:

  Only relevant if type == "dynamics". Specifies the relationship
  between `sub_features`, `sample_feature` (if specified), and the
  species in `graph`.

- graph:

  Only relevant if type == "dynamics". NxN adjacency matrix, where N
  represents the number of species modeled when running
  [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.md).

- normalize_by_median:

  Whether median difference was subtracted from estimated kinetic
  parameter difference. Thus, only relevant if type == "comparisons".

- deconvolved:

  Only relevant if type == "fractions". Boolean that is TRUE if
  fractions table is result of performing multi-feature deconvolution
  with
  [`DeconvolveFractions()`](https://isaacvock.github.io/EZbakR/reference/DeconvolveFractions.md).

- returnNameOnly:

  If TRUE, then only the names of tables that passed your search
  criteria will be returned. Else, the single table passing your search
  criteria will be returned. If there is more than one table that passes
  your search criteria and `returnNameOnly` == `FALSE`, an error will be
  thrown.

- exactMatch:

  If TRUE, then for `features` and `populations`, which can be vectors,
  ensure that provided vectors of features and populations exactly match
  the relevant metadata vectors.

- alwaysCheck:

  If TRUE, then even if there is only a single table for the `type` of
  interest, still run all checks against queries.

## Value

Table of interest from the relevant `EZbakRdata` list (set by the type
parameter).

## Details

The input to `EZget()` is 1) the type of table you want to get
("fractions", "kinetics", "averages", "comparisons", or "dynamics") and
2) the metadata necessary to uniquely specify the table of interest.
Above, every available piece of metadata that can be specified for this
purpose is documented. You only need to specify the minimum information
necessary. For example, if you would like to get a "fractions" table
from an analysis of exon bins (feature == "exon_bins", and potentially
other overarching features like "XF", "GF", or "rname"), and none of
your other "fractions" tables includes exon_bins as a feature, then you
can get this table with
`EZget(ezbdo, type = "fractions", features = "exon_bins")`, where
`ezbdo` is your `EZbakRData` object.

As another example, imagine you want to get a "kinetics" table from an
analysis of gene-wise kinetic parameters (e.g., features == "XF"). You
may have multiple "kinetics" tables, all with "XF" as at least one of
their features. If all of the other tables have additional features
though, then you can tell `EZget()` that "XF" is the only feature
present in your table of interest by setting `exactMatch` to `TRUE`,
which tells `EZget()` that the metadata you specify should exactly match
the relevant metadata for the table of interest. So the call in this
case would look like
`EZget(ezbdo, type = "fractions", features = "XF", exactMatch = TRUE)`.

`EZget()` is used internally in almost every single EZbakR function to
specify the input table for each analysis. Thus, the usage and metadata
described here also applies to all functions that require you to specify
which table you want to use (e.g.,
[`EstimateKinetics()`](https://isaacvock.github.io/EZbakR/reference/EstimateKinetics.md),
[`AverageAndRegularize()`](https://isaacvock.github.io/EZbakR/reference/AverageAndRegularize.md),
[`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.md),
etc.).

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

# Estimate Kinetics
ezbdo <- EstimateKinetics(ezbdo)

# Average log(kdeg) estimates across replicate
ezbdo <- AverageAndRegularize(ezbdo)
#> Fitting linear model
#> Estimating coverage vs. variance trend
#> Regularizing variance estimates

#' # Average log(ksyn) estimates across replicate
ezbdo <- AverageAndRegularize(ezbdo, parameter = "log_ksyn")
#> Fitting linear model
#> Estimating coverage vs. variance trend
#> Regularizing variance estimates

# Compare log(kdeg) across conditions
ezbdo <- CompareParameters(
ezbdo,
design_factor = "treatment",
reference = "treatment1",
experimental = "treatment2"
)

# Get the one and only fractions object
fxns <- EZget(ezbdo, type = "fractions")

# Get the log(ksyn) averages table
ksyn_avg <- EZget(ezbdo, type = "averages", parameter = "log_ksyn")
```
