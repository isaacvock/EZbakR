#' Standard `fraction_design` table for `EstimateFractions`
#'
#' An example `fraction_design` table for a standard NR-seq experiment with
#' s^4U labeling. This table tells `EstimateFractions` that there are two
#' populations of reads, one with high T-to-C mutation content and one
#' with low T-to-C mutation content
#'
#' @format ## `standard_fraction_design`
#' A tibble with 2 rows and 2 columns:
#' \describe{
#'   \item{TC}{Boolean denoting if population represented by that row has high T-to-C mutational content}
#'   \item{present}{Boolean denoting if population represented by that row is expected to be present in this dataset}
#' }
"standard_fraction_design"

#' TILAC `fraction_design` table for `EstimateFractions`
#'
#' An example `fraction_design` table for a TILAC experiment. TILAC was originally
#' described in [Courvan et al., 2022](https://academic.oup.com/nar/article/50/19/e110/6677324). In this
#' method, two populations of RNA, one from s^4U fed cells and one from s^6G fed cells, are pooled
#' and prepped for sequencing together. This allows for internally controlled comparisons of RNA
#' abundance without spike-ins. s^4U is recoded to a cytosine analog by TimeLapse chemistry
#' (or similar chemistry) and s^6G is recoded to an adenine analog. Thus, `fraction_design` includes
#' columns called `TC` and `GA`. A unique aspect of the TILAC `fraction_design` table is that
#' one of the possible populations, `TC` and `GA` both `TRUE`, is denoted as not present (`present` = `FALSE`).
#' This is because there is no RNA was exposed to both s^4U and s^6G, thus a population of reads
#' with both high T-to-C and G-to-A mutational content should not exist.
#'
#'
#' @format ## `tilac_fraction_design`
#' A tibble with 4 rows and 3 columns:
#' \describe{
#'   \item{TC}{Boolean denoting if population represented by that row has high T-to-C mutational content}
#'   \item{GA}{Boolean denoting if population represented by that row has high G-to-A mutational content}
#'   \item{present}{Boolean denoting if population represented by that row is expected to be present in this dataset}
#' }
"tilac_fraction_design"
