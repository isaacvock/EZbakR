# Vectorized simulation of one replicate of multi-label NR-seq data

Generalizes SimulateOneRep() to simulate any combination of mutation
types. Currently, no kinetic model is used to relate certain parameters
to the fractions of reads belonging to each simulated mutational
population. Instead these fractions are drawn from a Dirichlet
distribution with gene-specific parameters.

## Usage

``` r
VectSimulateMultiLabel(
  nfeatures,
  populations = c("TC"),
  fraction_design = create_fraction_design(populations),
  fractions_matrix = NULL,
  read_vect = NULL,
  sample_name = "sampleA",
  feature_prefix = "Gene",
  kdeg_vect = NULL,
  ksyn_vect = NULL,
  logkdeg_mean = -1.9,
  logkdeg_sd = 0.7,
  logksyn_mean = 2.3,
  logksyn_sd = 0.7,
  phighs = stats::setNames(rep(0.05, times = length(populations)), populations),
  plows = stats::setNames(rep(0.002, times = length(populations)), populations),
  seqdepth = nfeatures * 2500,
  readlength = 200,
  alpha_min = 3,
  alpha_max = 6,
  Ucont = 0.25,
  Acont = 0.25,
  Gcont = 0.25,
  Ccont = 0.25
)
```

## Arguments

- nfeatures:

  Number of "features" (e.g., genes) to simulate data for

- populations:

  Vector of mutation populations you want to simulate.

- fraction_design:

  Fraction design matrix, specifying which potential mutational
  populations should actually exist. See ?EstimateFractions for more
  details.

- fractions_matrix:

  Matrix of fractions of each mutational population to simulate. If not
  provided, this will be simulated. One row for each feature, one column
  for each mutational population, rows should sum to 1.

- read_vect:

  Vector of length = `nfeatures`; specifies the number of reads to be
  simulated for each feature. If this is not provided, the number of
  reads simulated is equal to
  `round(seqdepth * (ksyn_i/kdeg_i)/sum(ksyn/kdeg))`. In other words,
  the normalized steady-state abundance of a feature is multiplied by
  the total number of reads to be simulated and rounded to the nearest
  integer.

- sample_name:

  Character vector to assign to `sample` column of output simulated data
  table (the cB table).

- feature_prefix:

  Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows
  up in the `feature` column of the output simulated data table.

- kdeg_vect:

  Vector of length = `nfeatures`; specifies the degradation rate
  constant to use for each feature's simulation. If this is not provided
  and `fn_vect` is, then `kdeg_vect = -log(1 - fn_vect)/label_time`. If
  both `kdeg_vect` and `fn_vect` are not provided, each feature's
  `kdeg_vect` value is drawn from a log-normal distrubition with meanlog
  = `logkdeg_mean` and sdlog = `logkdeg_sd`. `kdeg_vect` is actually
  only simulated in the case where `read_vect` is also not provided, as
  it will be used to simulate read counts as described above.

- ksyn_vect:

  Vector of length = `nfeatures`; specifies the synthesis rate constant
  to use for each feature's simulation. If this is not provided, and
  `read_vect` is also not provided, then each feature's `ksyn_vect`
  value is drawn from a log-normal distribution with meanlog =
  `logksyn_mean` and sdlog = `logksyn_sd`. ksyn's do not need to be
  simulated if `read_vect` is provided, as they only influence read
  counts.

- logkdeg_mean:

  If necessary, meanlog of a log-normal distribution from which kdegs
  are simulated

- logkdeg_sd:

  If necessary, sdlog of a log-normal distribution from which kdegs are
  simulated

- logksyn_mean:

  If necessary, meanlog of a log-normal distribution from which ksyns
  are simulated

- logksyn_sd:

  If necessary, sdlog of a log-normal distribution from which ksyns are
  simulated

- phighs:

  Vector of probabilities of mutation rates in labeled reads of each
  type denoted in `populations`. Should be a named vector, with names
  being the corresponding `population`.

- plows:

  Vector of probabilities of mutation rates in unlabeled reads of each
  type denoted in `populations`. Should be a named vector, with names
  being the corresponding `population`.

- seqdepth:

  Only relevant if `read_vect` is not provided; in that case, this is
  the total number of reads to simulate.

- readlength:

  Length of simulated reads. In this simple simulation, all reads are
  simulated as being exactly this length.

- alpha_min:

  Minimum possible value of alpha element of Dirichlet random variable

- alpha_max:

  Maximum possible value of alpha element of Dirichlet random variable

- Ucont:

  Probability that a nucleotide in a simulated read is a U.

- Acont:

  Probability that a nucleotide in a simulated read is an A.

- Gcont:

  Probability that a nucleotide in a simulated read is a G.

- Ccont:

  Probability that a nucleotide in a simulated read is a C.

## Value

List with two elements:

- cB: Tibble that can be passed as the `cB` arg to
  [`EZbakRData()`](https://isaacvock.github.io/EZbakR/reference/EZbakRData.md).

- ground_truth: Tibble containing simulated ground truth.

## Examples

``` r
simdata <- VectSimulateMultiLabel(30)
```
