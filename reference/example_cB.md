# Example cB table

An example `cB` table used to create an `EZbakRData` object. This cB
table is a subset of a cB file from the DCP2 KO dataset published in Luo
et al., 2020. The original file is large (69 MB), so the example cB file
has been downsampled and contains only a subset of reads from chromosome
21.

## Usage

``` r
example_cB
```

## Format

### `example_cB`

A tibble with 10,000 rows and 7 columns:

- sample:

  Sample name

- rname:

  Chromosome name

- GF:

  Gene name for reads aligning to any region of a gene

- XF:

  Gene name for reads aligning to exclusively exonic regions of a gene

- TC:

  Number of T-to-C mutations

- nT:

  Number of Ts

- n:

  Number of reads with same value for the first 6 columns

## References

Luo et al. (2020) Biochemistry. 59(42), 4121-4142
