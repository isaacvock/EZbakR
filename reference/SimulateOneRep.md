# Simulate a single replicate of NR-seq data

In `SimulateOneRep`, users have the option to either provide vectors of
feature-specific read counts, fraction news, kdegs, and ksyns for the
simulation, or to have those drawn from relevant distributions whose
properties can be tuned by the various optional parameters of
`SimulateOneRep`. The number of mutable nucleotides (nT) in a read is
drawn from a binomial distribution with `readlength` trials and a
probability of "success" equal to `Ucont`. A read's status as new or old
is drawn from a Bernoulli distribution with probability of "success"
equal to the feature's fraction new. If a read is new, the number of
mutations in the read is drawn from a binomial distribution with
probability of mutation equal to pnew. If a read is old, the number of
mutations is instead drawn from a binomial distribution with probability
of mutation equal to pold.

## Usage

``` r
SimulateOneRep(
  nfeatures,
  read_vect = NULL,
  label_time = 2,
  sample_name = "sampleA",
  feature_prefix = "Gene",
  fn_vect = NULL,
  kdeg_vect = NULL,
  ksyn_vect = NULL,
  pnew = 0.05,
  pold = 0.002,
  logkdeg_mean = -1.9,
  logkdeg_sd = 0.7,
  logksyn_mean = 2.3,
  logksyn_sd = 0.7,
  seqdepth = nfeatures * 2500,
  readlength = 200,
  Ucont_alpha = 25,
  Ucont_beta = 75,
  feature_pnew = FALSE,
  pnew_kdeg_corr = FALSE,
  logit_pnew_mean = -2.5,
  logit_pnew_sd = 0.1
)
```

## Arguments

- nfeatures:

  Number of "features" (e.g., genes) to simulate data for

- read_vect:

  Vector of length = `nfeatures`; specifies the number of reads to be
  simulated for each feature. If this is not provided, the number of
  reads simulated is equal to
  `round(seqdepth * (ksyn_i/kdeg_i)/sum(ksyn/kdeg))`. In other words,
  the normalized steady-state abundance of a feature is multiplied by
  the total number of reads to be simulated and rounded to the nearest
  integer.

- label_time:

  Length of s^4^U feed to simulate.

- sample_name:

  Character vector to assign to `sample` column of output simulated data
  table (the cB table).

- feature_prefix:

  Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows
  up in the `feature` column of the output simulated data table.

- fn_vect:

  Vector of length = `nfeatures`; specifies the fraction new to use for
  each feature's simulation. If this is not provided and `kdeg_vect` is,
  then `fn_vect = 1 - exp(-kdeg_vect*label_time)`. If both `fn_vect` and
  `kdeg_vect` are not provided, then kdegs are simulated from a joint
  distribution as described below and converted to a `fn_vect` as when
  `kdeg_vect` is user-provided.

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

- pnew:

  Probability that a T is mutated to a C if a read is new.

- pold:

  Probability that a T is mutated to a C if a read is old.

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

- seqdepth:

  Only relevant if `read_vect` is not provided; in that case, this is
  the total number of reads to simulate.

- readlength:

  Length of simulated reads. In this simple simulation, all reads are
  simulated as being exactly this length.

- Ucont_alpha:

  Probability that a nucleotide in a simulated read from a given feature
  is a U is drawn from a beta distribution with shape1 = `Ucont_alpha`.

- Ucont_beta:

  Probability that a nucleotide in a simulated read from a given feature
  is a U is drawn from a beta distribution with shape2 = `Ucont_beta`.

- feature_pnew:

  Boolean; if TRUE, simulate a different pnew for each feature

- pnew_kdeg_corr:

  Boolean; only relevant if `feature_pnew` is TRUE. If so, then setting
  `pnew_kdeg_corr` to TRUE will ensure that higher kdeg transcripts have
  a higher pnew.

- logit_pnew_mean:

  If `feature_pnew` is TRUE, then the logit(pnew) for each feature will
  be drawn from a normal distribution with this mean.

- logit_pnew_sd:

  If `feature_pnew` is TRUE, then the logit(pnew) for each feature will
  be drawn from a normal distribution with this standard deviation.

## Examples

``` r
simdata <- SimulateOneRep(30)
```
