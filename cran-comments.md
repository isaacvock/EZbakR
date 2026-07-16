## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a patch release to address a newly exposed R CMD check NOTE:

  check_plabeled: no visible binding for global variable 'type'

The NOTE was identified during reverse-dependency checks for an upcoming `arrow` release.
EZbakR does not use `arrow::type()` directly; this patch declares `type` as a data-column symbol used in tidyverse-style code. 
There are no intentional user-facing API changes.