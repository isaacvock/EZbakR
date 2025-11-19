# Run quality control checks

Run quality control checks

## Usage

``` r
# S3 method for class 'EZbakRFractions'
EZQC(obj, features = NULL, populations = NULL, fraction_design = NULL, ...)
```

## Arguments

- obj:

  EZbakRFractions object, which is an EZbakRData object on which
  `EstimateFractions` has been run.

- features:

  Set of features analyzed in the fractions table you are interested
  QCing. This gets passed to
  [`EZget()`](https://isaacvock.github.io/EZbakR/reference/EZget.md) to
  help find this table.

- populations:

  Set of mutation types analyzed in the fractions table you are
  interested in QCing. This gets passed to
  [`EZget()`](https://isaacvock.github.io/EZbakR/reference/EZget.md) to
  help find this table.

- fraction_design:

  The fraction "design matrix" specified to get the fractions table you
  are interested in QCing. This gets passed to
  [`EZget()`](https://isaacvock.github.io/EZbakR/reference/EZget.md) to
  help find this table.

- ...:

  Additional arguments. Currently goes unused.

## Value

A list of `ggplot2` objects visualizing the various aspects of your data
assessed by
[`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.md).
