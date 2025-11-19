# `EZbakRKinetics` object constructor

`new_EZbakRKinetics` efficiently creates an object of class
`EZbakRKinetics`. It does not perform any rigorous checks of the
legitimacy of this object.

## Usage

``` r
new_EZbakRKinetics(
  kinetics,
  features,
  metadf,
  name = NULL,
  character_limit = 20
)
```

## Arguments

- kinetics:

  Data frame containing information about the kinetic parameters of
  interest for each set of features tracked.

- features:

  Features tracked in `kinetics` data frame. Needs to be specified
  explicitly as it cannot be automatically inferred.

- metadf:

  Data frame describing each of the samples included

- name:

  Optional; name to give to fractions table.

- character_limit:

  Maximum number of characters for naming out fractions output. EZbakR
  will try to name this as a "\_" separated character vector of all of
  the features analyzed. If this name is greater than `character_limit`,
  then it will default to "fraction#", where "#" represents a simple
  numerical ID for the table. in `kinetics`
