# Simple simulation function

Simple simulation function

## Usage

``` r
SimpleSim(
  nreads = 1000,
  fn = 0.5,
  pnew = 0.05,
  pold = 0.001,
  rlen = 100,
  Ucont = 0.25
)
```

## Arguments

- nreads:

  Number of reads to simulate

- fn:

  Fraction of reads that are new in simulation. Whether a read will be
  new will be determined by a draw from a Bernoulli(fn) distribution.

- pnew:

  T-to-C mutation rate in new reads

- pold:

  T-to-C mutation rate in old reads

- rlen:

  Length of simulated reads

- Ucont:

  Fraction of nucleotides in simulated reads that are Ts (U in RNA)

## Examples

``` r
# Simulate 1 gene worth of data data
simdata <- SimpleSim()
```
