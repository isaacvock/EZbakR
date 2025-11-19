# Run quality control checks

Run quality control checks

## Usage

``` r
# S3 method for class 'EZbakRArrowData'
EZQC(
  obj,
  mutrate_populations = "all",
  features = "all",
  filter_cols = "all",
  filter_condition = `&`,
  remove_features = c("NA", "__no_feature"),
  ...
)
```

## Arguments

- obj:

  An EZbakRData object

- mutrate_populations:

  Same as in
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
  See `?EstimateFractions()` for details.

- features:

  Same as in
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
  See `?EstimateFractions()` for details.

- filter_cols:

  Same as in
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
  See `?EstimateFractions()` for details.

- filter_condition:

  Same as in
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
  See `?EstimateFractions()` for details.

- remove_features:

  Same as in
  [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.md).
  See `?EstimateFractions()` for details.

- ...:

  Additional arguments. Currently goes unused.

## Value

A list of `ggplot2` objects visualizing the various aspects of your data
assessed by
[`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.md).
