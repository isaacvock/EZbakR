# `EZbakRArrowData` object constructor for internal use

`new_EZbakRArrowData` efficiently creates an object of class
`EZbakRArrowData`. It does not perform any rigorous checks of the
legitimacy of this object.

## Usage

``` r
new_EZbakRArrowData(cBds, metadf)
```

## Arguments

- cBds:

  Arrow dataset tracking the sample ID, mutational and nucleotide
  content, and feature assignment of sequencing reads.

- metadf:

  Data frame tracking features of each of the samples included in
  `cBds`.
