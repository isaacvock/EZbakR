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


#' Example cB table
#'
#' An example `cB` table used to create an `EZbakRData` object. This cB table is a
#' subset of a cB file from the DCP2 KO dataset published in Luo et al., 2020.
#' The original file is large (69 MB), so the example cB file has been
#' downsampled and contains only a subset of reads from chromosome 21.
#'
#' @format ## `example_cB`
#' A tibble with 10,000 rows and 7 columns:
#' \describe{
#'   \item{sample}{Sample name}
#'   \item{rname}{Chromosome name}
#'   \item{GF}{Gene name for reads aligning to any region of a gene}
#'   \item{XF}{Gene name for reads aligning to exclusively exonic regions of a gene}
#'   \item{TC}{Number of T-to-C mutations}
#'   \item{nT}{Number of Ts}
#'   \item{n}{Number of reads with same value for the first 6 columns}
#' }
#' @references Luo et al. (2020) Biochemistry. 59(42), 4121-4142
"example_cB"


#' Example metadf
#'
#' An example `metadf` table used to create an `EZbakRData` object. This metadf
#' describes the DCP2 KO dataset published in Luo et al., 2020.
#' @format ## `example_metadf`
#' A tibble with 6 rows and 3 columns
#' \describe{
#'   \item{sample}{Sample name}
#'   \item{tl}{Metabolic label feed time}
#'   \item{genotype}{Whether sample was collected from WT or DCP2 KO cells}
#' }
#' @references Luo et al. (2020) Biochemistry. 59(42), 4121-4142
"example_metadf"


#' Example ODE model graphs and formulas
#'
#' A list of example "graphs" and specie formulas that can be passed to
#' `EZDynamics` or `SimulateDynamics`, and that are used under the hood
#' by `EZSimulate` to facilitate simulations of ODE models.
#' @format ## `ode_models`
#' A list with 6 elements
#' \describe{
#'  \item{nuc2cyto}{Simplest model of nuclear and cytoplasmic RNA dynamics: 0 -> N -> C -> 0}
#'  \item{preRNA}{Simplest model of pre-RNA and mature RNA dynamics: 0 -> P -> M -> 0}
#'  \item{preRNAwithPdeg}{Same as preRNA, but now pre-RNA can also degrade.}
#'  \item{nuc2cytowithNdeg}{Same as nuc2cyto, but now nuclear RNA can also degrade.}
#'  \item{subtlseq}{Subcellular TimeLapse-seq model, similar to that described in Ietswaart et al., 2024.
#'  Simplest model discussed there, lacking nuclear degradation: 0 -> CH -> NP -> CY -> PL -> 0, and CY can
#'  also degrade.}
#'  \item{nuc2cytowithpreRNA}{Combination of nuc2cyto and preRNA where preRNA is first synthesized,
#'  then either processed or exported to the cytoplasm. Processing can also occur in the cytoplasm, and
#'  mature nuclear RNA can be exported to the cytoplasm. Only  mature RNA degrades.}
#' }
#' Each element of list has two items
#' \describe{
#'  \item{graph}{Matrix representation of ODE system graph.}
#'  \item{formulas}{Formula objects relating measured species to modeled species.}
#' }
#' @references Ietswaart et al. (2024) Molecular Cell. 84(14), 2765-2784.
"ode_models"
