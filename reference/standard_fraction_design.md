# Standard `fraction_design` table for `EstimateFractions`

An example `fraction_design` table for a standard NR-seq experiment with
s^4U labeling. This table tells `EstimateFractions` that there are two
populations of reads, one with high T-to-C mutation content and one with
low T-to-C mutation content

## Usage

``` r
standard_fraction_design
```

## Format

### `standard_fraction_design`

A tibble with 2 rows and 2 columns:

- TC:

  Boolean denoting if population represented by that row has high T-to-C
  mutational content

- present:

  Boolean denoting if population represented by that row is expected to
  be present in this dataset
