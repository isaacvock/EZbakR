# Get normalized read counts from either a cB table or `EZbakRFractions` object.

Uses TMM normalization strategy, similar to that used by DESeq2 and
edgeR.

## Usage

``` r
get_normalized_read_counts(
  obj,
  features_to_analyze,
  fractions_name = NULL,
  feature_lengths = NULL,
  scale_factors = NULL
)

# S3 method for class 'EZbakRFractions'
get_normalized_read_counts(
  obj,
  features_to_analyze,
  fractions_name = NULL,
  feature_lengths = NULL,
  scale_factors = NULL
)

# S3 method for class 'EZbakRData'
get_normalized_read_counts(
  obj,
  features_to_analyze,
  fractions_name = NULL,
  feature_lengths = NULL,
  scale_factors = NULL
)
```

## Arguments

- obj:

  An `EZbakRData` or `EZbakRFractions` object.

- features_to_analyze:

  Features in relevant table

- fractions_name:

  Name of fractions table to use

- feature_lengths:

  Table of effective lengths for each feature combination in your data.
  For example, if your analysis includes features named GF and XF, this
  should be a data frame with columns GF, XF, and length.

- scale_factors:

  Dataframe with two columns, one being "sample" (sample names) and the
  other being "scale_factor" (value to divide read counts by to
  normalize them)

## Value

Data table of normalized read counts.

## Methods (by class)

- `get_normalized_read_counts(EZbakRFractions)`: Method for class
  **EZbakRFractions** Get normalized read counts from fractions table.

- `get_normalized_read_counts(EZbakRData)`: Method for class
  **EZbakRData** Get normalized read counts from a cB table.

## Examples

``` r
# Simulate data
simdata <- EZSimulate(30)

# Create EZbakRData object
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)

# Get normalized read counts
reads <- get_normalized_read_counts(ezbdo, features_to_analyze = "feature")
```
