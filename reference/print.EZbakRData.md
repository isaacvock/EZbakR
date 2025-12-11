# Print method for `EZbakRData` objects

Print method for `EZbakRData` objects

## Usage

``` r
# S3 method for class 'EZbakRData'
print(x, max_name_chars = 60, ...)
```

## Arguments

- x:

  An `EZbakRData` object.

- max_name_chars:

  Maximum number of characters to print on each line

- ...:

  Ignored

## Value

The input `EZbakRData` object, invisibly

## Examples

``` r
# Simulate data to analyze
simdata <- SimulateOneRep(30)

# Create EZbakR input
metadf <- data.frame(sample = "sampleA", tl = 2)
ezbdo <- EZbakRData(simdata$cB, metadf)

# Print
print(ezbdo)
#> EZbakRData object with 0 analyses of 1 samples
```
