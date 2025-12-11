# Simulation of transcript isoform kinetic parameters.

`SimulateIsoforms()` performs a simple simulation of isoform-specific
kinetic parameters to showcase and test
[`EstimateIsoformFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateIsoformFractions.md).
It assumes that there are a set of reads (fraction of total set by
`funique` parameter) which map uniquely to a given isoform, while the
rest are ambiguous to all isoforms from that gene. Mutational content of
these reads are simulated as in
[`SimulateOneRep()`](https://isaacvock.github.io/EZbakR/reference/SimulateOneRep.md).

## Usage

``` r
SimulateIsoforms(
  nfeatures,
  nt = NULL,
  seqdepth = nfeatures * 2500,
  label_time = 4,
  sample_name = "sampleA",
  feature_prefix = "Gene",
  pnew = 0.1,
  pold = 0.002,
  funique = 0.2,
  readlength = 200,
  Ucont = 0.25,
  avg_numiso = 2,
  psynthdiff = 0.5,
  logkdeg_mean = -1.9,
  logkdeg_sd = 0.7,
  logksyn_mean = 2.3,
  logksyn_sd = 0.7
)
```

## Arguments

- nfeatures:

  Number of "features" to simulate data for. Each feature will have a
  simulated number of transcript isoforms

- nt:

  (Optional), can provide a vector of the number of isoforms you would
  like to simulate for each of the `nfeatures` features. Vector can
  either be length 1, in which case that many isoforms will be simulated
  for all features, or length equal to `nfeatures`.

- seqdepth:

  Total number of sequencing reads to simulate

- label_time:

  Length of s^4^U feed to simulate.

- sample_name:

  Character vector to assign to `sample` column of output simulated data
  table (the cB table).

- feature_prefix:

  Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows
  up in the `feature` column of the output simulated data table.

- pnew:

  Probability that a T is mutated to a C if a read is new.

- pold:

  Probability that a T is mutated to a C if a read is old.

- funique:

  Fraction of reads that uniquely "map" to a single isoform.

- readlength:

  Length of simulated reads. In this simple simulation, all reads are
  simulated as being exactly this length.

- Ucont:

  Percentage of nucleotides simulated to be U's.

- avg_numiso:

  Average number of isoforms for each feature. Feature-specific isoform
  counts are drawn from a Poisson distribution with this average. NOTE:
  to insure that all features have multiple isoforms, the simulated
  number of isoforms drawn from a Poisson distribution is incremented
  by 2. Thus, the actual average number of isoforms from each feature is
  `avg_numiso` + 2.

- psynthdiff:

  Percentage of genes for which all isoform abundance differences are
  synthesis driven. If not synthesis driven, then isoform abundance
  differences will be driven by differences in isoform kdegs.

- logkdeg_mean:

  meanlog of a log-normal distribution from which kdegs are simulated

- logkdeg_sd:

  sdlog of a log-normal distribution from which kdegs are simulated

- logksyn_mean:

  meanlog of a log-normal distribution from which ksyns are simulated

- logksyn_sd:

  sdlog of a log-normal distribution from which ksyns are simulated

## Value

List with two elements:

- cB: Tibble that can be passed as the `cB` arg to
  [`EZbakRData()`](https://isaacvock.github.io/EZbakR/reference/EZbakRData.md).

- ground_truth: Tibble containing simulated ground truth.

## Examples

``` r
simdata <- SimulateIsoforms(30)
```
