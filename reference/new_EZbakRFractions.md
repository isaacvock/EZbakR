# `EZbakRFractions` object constructor

`new_EZbakRFractions` efficiently creates an object of class
`EZbakRFractions`. It does not perform any rigorous checks of the
legitimacy of this object.

## Usage

``` r
new_EZbakRFractions(fractions, metadf, name = NULL, character_limit = 20)
```

## Arguments

- fractions:

  Data frame containing information about the fraction of reads from
  each mutational population of interest.

- metadf:

  Data frame reporting aspects of each of the samples included

- name:

  Optional; name to give to fractions table.

- character_limit:

  Maximum number of characters for naming out fractions output. EZbakR
  will try to name this as a "\_" separated character vector of all of
  the features analyzed. If this name is greater than `character_limit`,
  then it will default to "fraction#", where "#" represents a simple
  numerical ID for the table. in `fractions`
