# `EZbakRDataobject` constructor for internal use

`new_EZbakRData` efficiently creates an object of class `EZbakRData`. It
does not perform any rigorous checks of the legitimacy of this object.

## Usage

``` r
new_EZbakRData(cB, metadf)
```

## Arguments

- cB:

  Data frame tracking the sample ID, mutational and nucleotide content,
  and feature assignment of sequencing reads.

- metadf:

  Data frame tracking features of each of the samples included in `cB`.
