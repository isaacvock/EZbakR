## Resubmission

This is a resubmission. In this version I have:

* Reduced example runtimes and CPU usage in examples flagged by CRAN incoming checks.
* Confirmed Bioconductor repository information for Bioconductor dependencies.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a patch release to address a newly exposed R CMD check NOTE:

  check_plabeled: no visible binding for global variable 'type'

The NOTE was identified during reverse-dependency checks for an upcoming `arrow` release.
EZbakR does not use `arrow::type()` directly; this patch declares `type` as a data-column symbol used in tidyverse-style code. 
There are no intentional user-facing API changes.