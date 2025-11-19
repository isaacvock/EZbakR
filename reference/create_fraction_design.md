# Generate a `fraction_design` table for `EstimateFractions`.

A `fraction_design` table denotes what populations of labeled/unlabeled
RNA are present in your data. A `fraction_design` table as one column
for each mutation type (e.g., TC) present in your cB file, and one
column named "present". Each entry is either `TRUE` or `FALSE`. The rows
include all possible combinations of `TRUE` and `FALSE` for all mutation
types columns. A value of `TRUE` in a mutation type column represents a
population of reads that have high amounts (on average) of that mutation
type. For example, if your `fraction_design` table has mutation type
columns "TC" and "GA", the row with TC == `TRUE` and GA == `FALSE`
represents a population of reads with high T-to-C mutation content and
low G-to-A mutation content. In other words, these are reads from RNA
synthesized in the presence of s4U but not s6G. If such a population
exists in your data, the "present" column for that row should have a
value of `TRUE`.

## Usage

``` r
create_fraction_design(mutrate_populations)
```

## Arguments

- mutrate_populations:

  Character vector of the set of mutational populations present in your
  data. For example, s4U fed data with standard nucleotide recoding
  chemistry (e.g., TimeLapse, SLAM, TUC, AMUC, etc.) would have a
  `mutrate_populations` of c("TC"). Dual labeling experiments with s4U
  and s6G feeds would have a `mutrate_populations` of c("TC", "GA").

## Value

A `fraction_design` table that assumes that every possible combination
of mutational populations listed in `mutrate_populations` are present in
your data. The `present` column can be modified if this assumption is
incorrect. This default is chosen as it will in theory work for all
analyses, it may just be unnecessarily inefficient and estimate the
abundance of populations that don't exist.

## Examples

``` r
# Standard, single-label NR-seq
fd <- create_fraction_design(c("TC"))

# Dual-label NR-seq
fd2 <- create_fraction_design(c("TC", "GA"))

# Adjust dual-label output for TILAC
fd2$present <- ifelse(fd2$TC & fd2$GA, FALSE, fd2$present)

```
