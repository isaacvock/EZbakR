# Example ODE model graphs and formulas

A list of example "graphs" and specie formulas that can be passed to
`EZDynamics` or `SimulateDynamics`, and that are used under the hood by
`EZSimulate` to facilitate simulations of ODE models.

## Usage

``` r
ode_models
```

## Format

### `ode_models`

A list with 6 elements

- nuc2cyto:

  Simplest model of nuclear and cytoplasmic RNA dynamics: 0 -\> N -\> C
  -\> 0

- preRNA:

  Simplest model of pre-RNA and mature RNA dynamics: 0 -\> P -\> M -\> 0

- preRNAwithPdeg:

  Same as preRNA, but now pre-RNA can also degrade.

- nuc2cytowithNdeg:

  Same as nuc2cyto, but now nuclear RNA can also degrade.

- subtlseq:

  Subcellular TimeLapse-seq model, similar to that described in
  Ietswaart et al., 2024. Simplest model discussed there, lacking
  nuclear degradation: 0 -\> CH -\> NP -\> CY -\> PL -\> 0, and CY can
  also degrade.

- nuc2cytowithpreRNA:

  Combination of nuc2cyto and preRNA where preRNA is first synthesized,
  then either processed or exported to the cytoplasm. Processing can
  also occur in the cytoplasm, and mature nuclear RNA can be exported to
  the cytoplasm. Only mature RNA degrades.

Each element of list has two items

- graph:

  Matrix representation of ODE system graph.

- formulas:

  Formula objects relating measured species to modeled species.

## References

Ietswaart et al. (2024) Molecular Cell. 84(14), 2765-2784.
