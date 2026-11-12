#' Generate a `fraction_design` table for `EstimateFractions`.
#'
#' A `fraction_design` table denotes what populations of labeled/unlabeled RNA are present in your data.
#' A `fraction_design` table as one column for each mutation type (e.g., TC) present
#' in your cB file, and one column named "present". Each entry is either `TRUE` or
#' `FALSE`. The rows include all possible combinations of `TRUE` and `FALSE` for all
#' mutation types columns. A value of `TRUE` in a mutation type column represents
#' a population of reads that have high amounts (on average) of that mutation type.
#' For example, if your `fraction_design` table has mutation type columns "TC" and
#' "GA", the row with TC == `TRUE` and GA == `FALSE` represents a population of reads
#' with high T-to-C mutation content and low G-to-A mutation content. In other words,
#' these are reads from RNA synthesized in the presence of s4U but not s6G. If such
#' a population exists in your data, the "present" column for that row should have a
#' value of `TRUE`.
#'
#' @param mutrate_populations Character vector of the set of mutational populations
#' present in your data. For example, s4U fed data with standard nucleotide recoding
#' chemistry (e.g., TimeLapse, SLAM, TUC, AMUC, etc.) would have a `mutrate_populations`
#' of c("TC"). Dual labeling experiments with s4U and s6G feeds would have a `mutrate_populations`
#' of c("TC", "GA").
#' @return A `fraction_design` table that assumes that every possible combination of
#' mutational populations listed in `mutrate_populations` are present in your data.
#' The `present` column can be modified if this assumption is incorrect. This default
#' is chosen as it will in theory work for all analyses, it may just be unnecessarily
#' inefficient and estimate the abundance of populations that don't exist.
#' @examples
#'
#' # Standard, single-label NR-seq
#' fd <- create_fraction_design(c("TC"))
#'
#' # Dual-label NR-seq
#' fd2 <- create_fraction_design(c("TC", "GA"))
#'
#' # Adjust dual-label output for TILAC
#' fd2$present <- ifelse(fd2$TC & fd2$GA, FALSE, fd2$present)
#'
#'
#' @importFrom magrittr %>%
#' @export
create_fraction_design <- function(mutrate_populations){

  ### Hack to deal with devtools::check() NOTEs
  present <- NULL

  ### TO-DO: add custom types like TILAC and dual-label to create more realistic
  ### fraction design matrices with > 1 population.
  fraction_design <- dplyr::tibble(present = rep(TRUE,
                                                 times = 2^length(mutrate_populations)))

  fraction_list <- list()
  for(p in mutrate_populations){

    fraction_list[[p]] <- c(TRUE, FALSE)

  }

  fraction_design <- fraction_design %>%
    dplyr::bind_cols(dplyr::as_tibble(expand.grid(fraction_list))) %>%
    dplyr::select(!!mutrate_populations, present)

  return(fraction_design)

}




#' Estimate fractions of each RNA population.
#'
#' The first step of any NR-seq analysis is to figure out the fraction of reads
#' from each mutational population in your data. For example, if you are performing
#' a standard SLAM-seq or TimeLapse-seq experiment, this means estimating the fraction
#' of reads with high T-to-C mutation content, and the fraction with low T-to-C mutation
#' content. This is what `EstimateFractions` is for.
#'
#' `EstimateFractions` uses mixture modeling to estimate the fraction of reads from
#' each mutational population in your data, and this is done for each feature in your
#' data (i.e., combination of columns that specify genomic features from which reads
#' were derived). The set of mutational populations in your
#' data can be specified by providing a `fraction_design` object, described in more
#' depth above (also see `?create_fraction_design`). There are several flavors
#' of mixture modeling that can be performed by `EstimateFractions`. These are
#' as follows:
#'
#' 1. The default: global mutation rate parameters (global = same for all reads
#' from all features) are estimated for each sample by fitting a single two-component
#' mixture model to all of the reads in that sample. These are used to estimate
#' the fraction of reads from each feature that are from each mutational population,
#' also using a two-component mixture model.
#'    - With `Poisson` set to `TRUE`, this is a nucleotide-content adjusted Poisson
#'  mixture model, which is a more efficient alternative to binomial mixture modeling.
#'    - With `Poisson` set to `FALSE`, this is a binomial mixture model.
#' 2. Low mutation rate from -label: if `pold_from_nolabel` is `TRUE`, then the background,
#' no label mutation rates are estimated from -label samples. By default, a single set of
#' background mutation rates are estimated for all samples, but you can change this behavior
#' by setting `grouping_factors` to specify columns in your `metadf` by which samples
#' should be stratified.
#'    - This is a great strategy to use if you have low label incorporation rates or
#'    if you used a fairly short label time.
#' 3. Hierarchical model: if `strategy == "hierarchical"`, which is currently
#' only compatible for single-mutation type modeling (e.g., standard T-to-C mutation
#' modeling), then high T-to-C content mutation rates are estimated for each feature.
#' The global sample-wide estimates are used as an informative prior to increase the
#' accuracy of this process by avoiding extreme estimates.
#'
#' @param obj `EZbakRData` or `EZbakRArrowData` object
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB file.
#' @param mutrate_populations Character vector of the set of mutational populations
#' that you want to infer the rates of mutations for. By default, all mutation rates
#' are estimated for all populations present in cB.
#' @param fraction_design "Design matrix" specifying which RNA populations exist
#' in your samples. By default, this will be created automatically and will assume
#' that all combinations of the `mutrate_populations` you have requested to analyze are
#' present in your data. If this is not the case for your data, then you will have
#' to create one manually.
#'
#' If you call the function \code{create_fraction_design(...)}, providing a vector
#' of mutational population names as input, it will create a `fraction_design` table for
#' you, with the assumption that every single possible combination of mutational populations
#' is present in your data. You can then edit the `present` column as necessary to
#' get an appropriate `fraction_design` for your use case. See below for details on
#' the required contents of `fraction_design` and its interpretation.
#'
#' `fraction_design` must have one column per element of `mutrate_populations`,
#' with these columns sharing the name of the `mutrate_populations`. It must also have
#' one additional column named `present`. All elements of fraction_design should be
#' booleans (`TRUE` or `FALSE`). It should include all possible combinations of `TRUE`
#' and `FALSE` for the `mutrate_populations` columns. A `TRUE` in one of these columns
#' represents a population of RNA that is expected to have above background mutation
#' rates of that type. `present` will denote whether or not that population of RNA
#' is expected to exist in your data.
#'
#' For example, assume you are doing a typical TimeLapse-seq/SLAM-seq/TUC-seq/etc. experiment
#' where you have fed cells with s^4U and recoded any incorporated s^4U to a nucleotide
#' that reverse transcriptase will read as a cytosine. That means that `mutrate_populations`
#' will be "TC", since you want to estimate the fraction of RNA that was s^4U labeled, i.e.,
#' the fraction with high T-to-C mutation content. `fraction_design` will thus have two columns:
#' `TC` and `present`. It will also have two rows. One of these rows must have a value of
#' `TRUE` for `TC`, and the other must have a value of `FALSE`. The row with a value of `TRUE`
#' for `TC` represents the population of reads with high T-to-C mutation content,
#' i.e., the reads from RNA that were synthesized while s^4U was present. The row
#' with a value of `FALSE` for `TC` reprsents the population of reads with low T-to-C mutation
#' content, i.e., the reads from RNA that existed prior to s^4U labeling. Both of these
#' populations exist in your data, so the value of the `present` column should be `TRUE` for
#' both of these. See the lazily loaded `standard_fraction_design` object for an example
#' of what this tibble could look like. ("lazily loaded `standard_fraction_design` object"
#' means that if you run \code{print(standard_fraction_design)} after loading `EZbakR` with \code{library(EZbakR)},
#' then you can see its contents. More specifically, lazily loaded means that this table is not loaded
#' into memory until you ask for it, via something like a \code{print()} call.)
#'
#' As another example, consider TILAC, a NR-seq extension developed by the Simon lab. TILAC was originally
#' described in Courvan et al., 2022 (\doi{10.1093/nar/gkac693}). In this
#' method, two populations of RNA, one from s^4U fed cells and one from s^6G fed cells, are pooled
#' and prepped for sequencing together. This allows for internally controlled comparisons of RNA
#' abundance without spike-ins. s^4U is recoded to a cytosine analog by TimeLapse chemistry
#' (or similar chemistry) and s^6G is recoded to an adenine analog. Thus, `fraction_design` includes
#' columns called `TC` and `GA`. A unique aspect of the TILAC `fraction_design` table is that
#' one of the possible populations, `TC` and `GA` both `TRUE`, is denoted as not present (`present` = `FALSE`).
#' This is because there is no RNA that was exposed to both s^4U and s^6G, thus a population of reads
#' with both high T-to-C and G-to-A mutational content should not exist. To see an example
#' of what a TILAC `fraction_design` table could look like, see the lazily loaded
#' `tilac_fraction_design` object.
#'
#' @param Poisson If `TRUE`, use U-content adjusted Poisson mixture modeling strategy. Often provides
#' significant speed gain without sacrificing accuracy.
#' @param strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item hierarchical: Estimate feature-specific new read mutation
#'  rate, regularizing the feature-specific estimate with a sample-wide prior. Currently
#'  only compatible with single mutation type mixture modeling.
#'  }
#' @param filter_cols Which feature columns should be used to filter out feature-less
#' reads. The default value of "all" checks all feature columns for whether or not a
#' read failed to get assigned to said feature.
#' @param filter_condition Only two possible values for this make sense: \code{`&`} and \code{`|`}.
#' If set to \code{`&`}, then all features in `filter_cols` must have a "null" value (i.e., a value included
#' in `remove_features`) for the row to get filtered out. If set to \code{`|`}, then only
#' a single feature in `filter_cols` needs to have one of these "null" values to get filtered out.
#' @param remove_features All of the feature names that could indicate failed assignment of a read
#' to a given feature. the fastq2EZbakR pipeline uses a value of '__no_feature'.
#' @param remove_multi_features Logical; if `TRUE`, filter out reads assigned
#' ambiguously to multiple features (identified by containing `multi_feature_sep`,
#' e.g., `"featureA+featureB"`). Defaults to `FALSE` for backwards compatibility.
#' By default (`multi_feature_cols = NULL`), all analyzed `features` columns are
#' checked; pass a character vector to `multi_feature_cols` to restrict multi-feature
#' filtering to specific feature columns (e.g., `multi_feature_cols = "XF"` when
#' analyzing `features = c("XF", "TEC")` so multi-gene reads are dropped while
#' multi-isoform equivalence classes are retained). Cannot be `TRUE` at the same
#' time as `split_multi_features`.
#' @param split_multi_features If a set of reads maps ambiguously to multiple features,
#' should data for such reads be copied for each feature in the ambiguous set? If this is
#' `TRUE`, then `multi_feature_cols` also must be set. Examples where this should be set
#' to `TRUE` includes when analyzing exonic bins (concept defined in original DEXSeq paper),
#' exon-exon junctions, etc. Cannot be `TRUE` at the same time as `remove_multi_features`.
#' @param multi_feature_cols Character vector of columns that have the potential to
#' include assignment to multiple features. When `split_multi_features = TRUE`, this
#' argument is required and specifies which columns have their features split. When
#' `remove_multi_features = TRUE`, this optionally restricts which analyzed feature
#' columns are checked for multi-feature assignments (`NULL` or `"all"` checks all
#' analyzed features).
#' @param multi_feature_sep String representing how ambiguous feature assignments are
#' distinguished in the feature names (used by both `remove_multi_features` and
#' `split_multi_features`). For example, the default value of `"+"` denotes
#' that if a read maps to multiple features (call them featureA and featureB, for example),
#' then the feature column will have a value of `"featureA+featureB"` for that read.
#' @param pnew_prior_mean Mean for logit(pnew) prior.
#' @param pnew_prior_sd Standard deviation for logit(pnew) prior.
#' @param pold_prior_mean Mean for logit(pold) prior.
#' @param pold_prior_sd Standard deviation for logit(pold) prior.
#' @param hier_readcutoff If `strategy` == `hierarchical`, only features with this many reads
#' are used to infer the distribution of feature-specific labeled read mutation rates.
#' @param init_pnew_prior_sd If `strategy` == `hierarchical`, this is the initial logit(pnew)
#' prior standard deviation to regularize feature-specific labeled read mutation rate estimates.
#' @param pnew_prior_sd_min The minimum logit(pnew) prior standard deviation when `strategy`
#' is set to "hierarchcial". EZbakR will try to estimate this empirically as the standard deviation
#' of initial feature-specific logit(pnew) estimates using high coverage features, minus the
#' average uncertainty in the logit(pnew) estimates. As this difference can sometimes be negative,
#' a value of `pnew_prior_sd_min` will be imputed in that case.
#' @param pnew_prior_sd_max Similar to `pnew_prior_sd_min`, but now representing the
#' maximum allowed logit(pnew) prior sd.
#' @param pold_est Background mutation rate estimates if you have them. Can either be a single
#' number applied to all samples or a named vector of values, where the names should be sample
#' names.
#' @param pold_from_nolabel Fix background mutation rate estimate to mutation rates seen in -label samples.
#' By default, a single background rate is used for all samples, inferred from the average mutation rate
#' across all -label samples. The `grouping_factors` argument can be specified to use certain -label samples
#' to infer background mutation rates for certain sets of +label samples.
#' @param grouping_factors If `pold_from_nolabel` is TRUE, then `grouping_factors` will specify the
#' sample-detail columns in the metadf that should be used to group -label samples by. Average mutation
#' rates in each group of -label samples will be used as the background mutation rate estimate in
#' +label samples with the same values for the relevant metadf columns.
#' @param EBayes Logical; if `TRUE` (the default), calibrate the sample-specific
#' `logit(fraction)` prior location and scale empirically from the pure-data
#' maximum likelihood estimates of interior features in each labeled sample
#' (with heavy-tailed Pseudo-Huber regularization and profile-likelihood SE
#' flooring). If `FALSE`, use the legacy fixed `Normal(0, 1.7)` prior.
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @return An `EZbakRFractions` object, which is just an `EZbakRData` object
#' with a "fractions" list element. This will include tables of fraction estimates
#' (e.g., fraction of reads from the high T-to-C mutation rate population in a
#' standard single-label s4U experiment; termed fraction_highTC in the table) for
#' all features in all samples.
#' @import data.table
#' @importFrom magrittr %>%
#' @examples
#'
#' # Simulate data to analyze
#' simdata <- SimulateOneRep(30)
#'
#' # Create EZbakR input
#' metadf <- data.frame(sample = "sampleA", tl = 2)
#' ezbdo <- EZbakRData(simdata$cB, metadf)
#'
#' # Estimate fractions
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' @export
EstimateFractions <- function(obj, features = "all",
                              mutrate_populations = "all",
                              fraction_design = NULL,
                              Poisson = TRUE,
                              strategy = c("standard", "hierarchical"),
                              filter_cols = "all",
                              filter_condition = `&`,
                              remove_features = c("NA", "__no_feature"),
                              remove_multi_features = FALSE,
                              split_multi_features = FALSE,
                              multi_feature_cols = NULL,
                              multi_feature_sep = "+",
                              pnew_prior_mean = -2.94,
                              pnew_prior_sd = 0.3,
                              pold_prior_mean = -6.5,
                              pold_prior_sd = 0.5,
                              hier_readcutoff = 300,
                              init_pnew_prior_sd = 0.3,
                              pnew_prior_sd_min = 0.01,
                              pnew_prior_sd_max = 0.15,
                              pold_est = NULL,
                              pold_from_nolabel = FALSE,
                              grouping_factors = NULL,
                              EBayes = TRUE,
                              character_limit = 20,
                              overwrite = TRUE){


  UseMethod("EstimateFractions")


}



#' Estimate mutation rates
#'
#' Two component mixture models are fit to all data to estimate global high and
#' low mutation rates for all samples. Estimation of these mutation rates
#' are regularized through the use of weakly informative priors whose parameters can
#' be altered using the arguments defined below.
#'
#' @param obj An `EZbakRData` or `EZbakRArrowData` object
#' @param populations Character vector of the set of mutational populations
#' that you want to infer the fractions of. For example, say your cB file contains
#' columns tracking T-to-C and G-to-A
#' @param pnew_prior_mean logit-Normal mean for logit(pnew) prior.
#' @param pnew_prior_sd logit-Normal sd for logit(pnew) prior.
#' @param pold_prior_mean logit-Normal mean for logit(pold) prior.
#' @param pold_prior_sd logit-Normal sd for logit(pold) prior.
#' @param pold_est Background mutation rate estimates if you have them. Can either be a single
#' number applied to all samples or a named vector of values, where the names should be sample
#' names.
#' @param pold_from_nolabel Fix background mutation rate estimate to mutation rates seen in -label samples.
#' By default, a single background rate is used for all samples, inferred from the average mutation rate
#' across all -label samples. The `grouping_factors` argument can be specified to use certain -label samples
#' to infer background mutation rates for certain sets of +label samples.
#' @param grouping_factors If `pold_from_nolabel` is TRUE, then `grouping_factors` will specify the
#' sample-detail columns in the metadf that should be used to group -label samples by. Average mutation
#' rates in each group of -label samples will be used as the background mutation rate estimate in
#' +label samples with the same values for the relevant metadf columns.
#' @return `EZbakRData` object with an added `mutrates` slot containing estimated
#' high and low mutation rates for each mutation type modeled.
#' @import data.table
#' @examples
#'
#' # Simulate data to analyze
#' simdata <- SimulateOneRep(30)
#'
#' # Create EZbakR input
#' metadf <- data.frame(sample = "sampleA", tl = 2)
#' ezbdo <- EZbakRData(simdata$cB, metadf)
#'
#' # Estimate mutation rates
#' mutrates <- EstimateMutRates(ezbdo)
#'
#' @export
EstimateMutRates <- function(obj,
                            populations = "all",
                            pnew_prior_mean = -2.94,
                            pnew_prior_sd = 0.3,
                            pold_prior_mean = -6.5,
                            pold_prior_sd = 0.5,
                            pold_est = NULL,
                            pold_from_nolabel = FALSE,
                            grouping_factors = NULL
){

  UseMethod("EstimateMutRates")

}





###################################################
# EZBAKRDATA STANDARD METHODS
###################################################



#' @describeIn EstimateFractions Method for class **EZbakRData**
#' Estimates fractions using an entirely in memory object.
#' @export
EstimateFractions.EZbakRData <- function(obj, features = "all",
                                         mutrate_populations = "all",
                                         fraction_design = NULL,
                                         Poisson = TRUE,
                                         strategy = c("standard", "hierarchical"),
                                         filter_cols = "all",
                                         filter_condition = `&`,
                                         remove_features = c("NA", "__no_feature"),
                                         remove_multi_features = FALSE,
                                         split_multi_features = FALSE,
                                         multi_feature_cols = NULL,
                                         multi_feature_sep = "+",
                                         pnew_prior_mean = -2.94,
                                         pnew_prior_sd = 0.3,
                                         pold_prior_mean = -6.5,
                                         pold_prior_sd = 0.5,
                                         hier_readcutoff = 300,
                                         init_pnew_prior_sd = 0.3,
                                         pnew_prior_sd_min = 0.01,
                                         pnew_prior_sd_max = 0.15,
                                         pold_est = NULL,
                                         pold_from_nolabel = FALSE,
                                         grouping_factors = NULL,
                                         EBayes = TRUE,
                                         character_limit = 20,
                                         overwrite = TRUE){

  ### Hack to deal with devtools::check() NOTEs
  n <- params <- coverage <- pold <- pnew <- lpnew_uncert <- pnew_prior <- ..colvect <- NULL
  fit <- mixture_fit <- present <- ci_lower_fn <- ci_upper_fn <- NULL

  rm(..colvect)

  pnew_prior_sd_min <- pnew_prior_sd_min
  pnew_prior_sd_max <- pnew_prior_sd_max

  `.` <- list

  ### Check that input is valid

  args <- c(as.list(environment()))


  check_EstimateFractions_input(args)


  strategy <- match.arg(strategy)


  ### Estimate mutation rates
  message("Estimating mutation rates")
  obj <- EstimateMutRates(obj,
                           populations = mutrate_populations,
                          pnew_prior_mean = pnew_prior_mean,
                          pnew_prior_sd = pnew_prior_sd,
                          pold_prior_mean = pold_prior_mean,
                          pold_prior_sd = pold_prior_sd,
                          pold_est = pold_est,
                          pold_from_nolabel = pold_from_nolabel,
                          grouping_factors = grouping_factors)

  mutation_rates <- obj$mutation_rates

  ### Figure out which features will be analyzed

  # cB columns
  cB <- obj$cB
  cB_cols <- colnames(cB)

  # Mutation columns in cB
  mutcounts_in_cB <- find_mutcounts(obj)

  # Base count columns in cB
  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # feature columns
  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB,
                                           basecounts_in_cB,
                                           "sample", "n"))]

  # Need to determine which columns of the cB to group reads by
  if(features[1] == "all" & length(features) == 1){

    features_to_analyze <- features_in_cB

  }else{

    if(!all(features %in% features_in_cB)){

      stop("features includes columns that do not exist in your cB!")

    }else{

      features_to_analyze <- features

    }

  }


  ### Figure out which mutational populations to analyze

  if(length(mutrate_populations) == 1 && mutrate_populations == "all"){

    pops_to_analyze <- mutcounts_in_cB

  }else{

    if(!all(mutrate_populations %in% mutcounts_in_cB)){

      stop("mutrate_populations includes columns that do not exist in your cB!")

    }else{

      pops_to_analyze <- mutrate_populations

    }
  }

  # Get nucleotide counts that are needed
  necessary_basecounts <- paste0("n", substr(pops_to_analyze, start = 1, stop = 1))

  ### Create fraction_design data frame if necessary

  if(is.null(fraction_design)){

    fraction_design <- create_fraction_design(pops_to_analyze)

  }


  ### Filter out label-less samples

  metadf <- obj$metadf

  # What columns represent label times of interest?
  if(length(pops_to_analyze) == 1){

    tl_cols_possible <- c("tl", "tpulse")

  }else{

    tl_cols_possible <- c(paste0("tl_", pops_to_analyze),
                          paste0("tpulse_", pops_to_analyze))


  }

  tl_cols <- tl_cols_possible[tl_cols_possible %in% colnames(metadf)]

  ### Which samples should be filtered out

  metadf <- obj$metadf

  samples_with_no_label <- metadf %>%
    dplyr::rowwise() %>%
    dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) == 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample) %>%
    unlist() %>%
    unname()


  ### Estimate fraction new for each feature in each sample
  message("Summarizing data for feature(s) of interest")

  cB <- data.table::as.data.table(cB)

  # Keep only feature of interest
  cols_to_group <- c(necessary_basecounts, pops_to_analyze, "sample", features_to_analyze)
  cB <- cB[sample %in% metadf$sample][,.(n = sum(n)), by = cols_to_group]

  ### Filter out columns not mapping to any feature (easier and faster in data.table)
  if(filter_cols[1] == "all" & length(filter_cols) == 1){

    filter_cols <- features_to_analyze

  }

  cB <- cB[cB[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ]

  ### Filter out multi-feature mappers if necessary
  multi_feature_cols_to_use <- multi_feature_cols
  if(remove_multi_features){

    if(is.null(multi_feature_cols) ||
       (length(multi_feature_cols) == 1L && identical(multi_feature_cols, "all"))){
      multi_feature_cols_to_use <- features_to_analyze
    }else{
      if(!all(multi_feature_cols %in% features_to_analyze)){
        stop("multi_feature_cols includes columns that are not in the features being analyzed!")
      }
      multi_feature_cols_to_use <- multi_feature_cols
    }

    cB <- remove_multi_feature_rows(cB, multi_feature_cols_to_use, sep = multi_feature_sep)

    if(nrow(cB) == 0L){
      stop("No reads remain after filtering with remove_multi_features = TRUE!")
    }

  }


  # Pair down design matrix
  fraction_design <- dplyr::as_tibble(fraction_design)
  mutrate_design <- fraction_design %>%
    dplyr::filter(present) %>%
    dplyr::select(-present)

  # Summarize data
  if(Poisson){

    message("Averaging out the nucleotide counts for improved efficiency")
    cols_to_group_pois <- c("sample", pops_to_analyze, features_to_analyze)
    cols_to_avg <- setdiff(names(cB), c(cols_to_group_pois, "n"))

    cB <- cB[, c(lapply(.SD, function(x) sum(x*n)/sum(n) ),
                     .(n = sum(n))),
                 by = cols_to_group_pois, .SDcols = cols_to_avg]

  }


  ### Split multi feature mappers if necessary

  if(split_multi_features){

    if(!all(multi_feature_cols %in% features_to_analyze)){
      stop("multi_feature_cols includes columns that are not in the features being analyzed!")
    }

    cB <- split_features(cB, multi_feature_cols, c(cols_to_group, necessary_basecounts), sep = multi_feature_sep)

  }

  hier_mutrate_name <- NULL
  hier_prior_df <- NULL

  # TO-DO: Add an estimated runtime here based on the number of rows in the cB,
  # the number of features, and what model is being run.
  message("Estimating fractions")
  if(ncol(mutrate_design) == 1){

    ### Can optimize for two-component mixture model

    if(unlist(mutrate_design[1,1])){

      highpop <- 1

    }else{

      highpop <- 2

    }


    ### Add mutation rate info

    # Extract mutation rates
    mutrates <- setDT(obj$mutation_rates[[1]] %>%
      dplyr::select(-params))

    # Key for fast join
    setkey(mutrates, sample)
    setkey(cB, sample)


    # Join
    cB <- mutrates[cB]



    ### Estimate fraction new

    col_name <- paste0("logit_fraction_high", pops_to_analyze)
    uncertainty_col <- paste0("se_logit_fraction_high", pops_to_analyze)
    natural_col_name <- paste0("fraction_high", pops_to_analyze)
    ci_lower_col <- paste0("ci_lower_fraction_high", pops_to_analyze)
    ci_upper_col <- paste0("ci_upper_fraction_high", pops_to_analyze)

    if(strategy == "hierarchical"){

      message("FITTING HIERARCHICAL TWO-COMPONENT MIXTURE MODEL:")

      ### Get coverages to filter by; only want high coverage for feature-specific
      ### pnew estimation

      coverages <- cB[!(sample %in% samples_with_no_label)][,.(coverage = sum(n)),
                      by = c("sample", features_to_analyze)]
      coverages <- coverages[coverage > hier_readcutoff][, c("coverage") := .(NULL)]


      ### Get feature-specific pnew estimates

      message("Estimating distribution of feature-specific pnews")

      keyvector <- c("sample", features_to_analyze)
      setkeyv(coverages, keyvector)
      setkeyv(cB, keyvector)

      cB_highcov <- cB[coverages, nomatch = NULL]
      cB_highcov[, pnew_prior_init := logit(pnew)]

      feature_specific <- solve_tcmm_hier_dt(
        dt = cB_highcov,
        features_to_analyze = features_to_analyze,
        pop_col = pops_to_analyze,
        base_col = necessary_basecounts,
        samples_with_no_label = samples_with_no_label,
        Poisson = Poisson,
        pnew_prior_mean_col = "pnew_prior_init",
        pnew_prior_sd_arg = init_pnew_prior_sd,
        fraction_prior_sd = 1
      )

      feature_specific[, c("pnew", "lpnew_uncert") := .(inv_logit(p2), p2_u)]


      ### Hierarchical model using priors inferred from feature-specific pnew distribution

      message("Estimating fractions with feature-specific pnews")

      ## Mean should really be the estimated sample wide value
      global_est <- estimate_hier_pnew_prior(
        feature_specific = feature_specific,
        mutrates = mutrates,
        pnew_prior_sd_min = pnew_prior_sd_min,
        pnew_prior_sd_max = pnew_prior_sd_max
      )

      setkey(global_est, sample)
      setkey(cB, sample)

      cB_joined <- global_est[, .(sample, pnew_prior, pnew_prior_sd)][cB]

      feature_specific <- solve_tcmm_hier_dt(
        dt = cB_joined,
        features_to_analyze = features_to_analyze,
        pop_col = pops_to_analyze,
        base_col = necessary_basecounts,
        samples_with_no_label = samples_with_no_label,
        Poisson = Poisson,
        pnew_prior_mean_col = "pnew_prior",
        pnew_prior_sd_arg = "pnew_prior_sd",
        fraction_prior_sd = 1.7,
        EBayes = EBayes,
        compute_ci = TRUE
      )

      feature_specific[, c(col_name, "pnew",
                           uncertainty_col) := .(p1,
                                                 inv_logit(p2),
                                                 p1_u)]


      fns <- dplyr::as_tibble(feature_specific) %>%
        dplyr::mutate(!!natural_col_name := inv_logit(!!dplyr::sym(col_name)),
                      !!ci_lower_col := ci_lower_fn,
                      !!ci_upper_col := ci_upper_fn) %>%
        dplyr::select(sample, !!features_to_analyze, !!natural_col_name,
                      !!col_name, !!uncertainty_col, !!ci_lower_col, !!ci_upper_col, n)


      ### Save feature-specific pnews

      mutrates[, pnew := NULL]

      colvect <- c("sample", features_to_analyze, "pnew")

      feature_mutrates <- feature_specific[,..colvect][
             mutrates, on = "sample", nomatch = NULL
        ]

      current_name <- names(obj$mutation_rates)[1]
      hier_mutrate_name <- paste0(paste(features_to_analyze, collapse = "_"), "_", current_name)
      hier_prior_df <- dplyr::as_tibble(global_est[, .(sample, pnew_prior, pnew_prior_sd)])
      obj$mutation_rates[[hier_mutrate_name]] <- feature_mutrates



    }else{

      fns <- solve_tcmm_standard_dt(
        dt = cB,
        features_to_analyze = features_to_analyze,
        pop_col = pops_to_analyze,
        base_col = necessary_basecounts,
        samples_with_no_label = samples_with_no_label,
        Poisson = Poisson,
        col_name = col_name,
        uncertainty_col = uncertainty_col,
        natural_col_name = natural_col_name,
        EBayes = EBayes
      )

    }





  }else{


    # pops_to_analyze needs to match order as it appears in mutrate_design
    # and necessary_basecounts order needs to match relevant nucleotide as well
    pops_to_analyze <- colnames(mutrate_design)
    necessary_basecounts <- paste0("n", substr(pops_to_analyze, start = 1, stop = 1))

    fns <- solve_general_mixture_dt(
      dt = cB,
      features_to_analyze = features_to_analyze,
      mutrate_design = mutrate_design,
      mutcols = pops_to_analyze,
      basecols = necessary_basecounts,
      Poisson = Poisson,
      samples_with_no_label = samples_with_no_label,
      mutation_rates = mutation_rates
    )

  }


  ##### PROCESS OUTPUT AND RETURN


  message("Processing output")



  # What should output be named?
  fraction_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")

  if(nchar(fraction_vect) > character_limit){

    num_fractions <- length(obj[['fractions']])
    fraction_vect <- paste0("fractions", num_fractions + 1)

  }


  ### CONCERN: Worried about what happens if `character_limit`
  ### changes from run to run.
  ### ALLAYED: Behavior of `decide_output()` is such that
  ### it checks if there are identical existing output, and changes
  ### its name to that output rather than using the proposed name

  # Are there any metadata or fractions objects at this point?
  if(length(obj[['metadata']][['fractions']]) > 0){

    fraction_vect <- decide_output(obj,
                                   proposed_name = fraction_vect,
                                   type = "fractions",
                                   features = features_to_analyze,
                                   populations = pops_to_analyze,
                                   fraction_design = fraction_design,
                                   overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'fractions',
                               features = features_to_analyze,
                               populations = pops_to_analyze,
                               fraction_design = fraction_design,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }


  # Save
  obj[['fractions']][[fraction_vect]] <- dplyr::as_tibble(fns)

  # Save metadata
  obj[['metadata']][['fractions']][[fraction_vect]] <- list(features = features_to_analyze,
                                                            populations = pops_to_analyze,
                                                            fraction_design = fraction_design,
                                                            strategy = if(ncol(mutrate_design) == 1) strategy else "general",
                                                            Poisson = Poisson,
                                                            remove_multi_features = remove_multi_features,
                                                            split_multi_features = split_multi_features,
                                                            multi_feature_cols = multi_feature_cols_to_use,
                                                            multi_feature_sep = multi_feature_sep,
                                                            filter_cols = filter_cols,
                                                            remove_features = remove_features,
                                                            hier_mutrate_table = hier_mutrate_name,
                                                            hier_prior = hier_prior_df,
                                                            repeatID = repeatID)

  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)

}




#' @describeIn EstimateMutRates Method for class **EZbakRData**
#' Estimates mutation rates using a fully in memory object.
#' @export
EstimateMutRates.EZbakRData <- function(obj,
                                        populations = "all",
                                        pnew_prior_mean = -2.94,
                                        pnew_prior_sd = 0.3,
                                        pold_prior_mean = -6.5,
                                        pold_prior_sd = 0.5,
                                        pold_est = NULL,
                                        pold_from_nolabel = FALSE,
                                        grouping_factors = NULL
                             ){


  ### Hack to deal with devtools::check() NOTEs
  n <- params <- p1 <- p2 <- ..cols_to_keep <- NULL

  rm(..cols_to_keep)

  `.` <- list


  cB <- obj$cB
  metadf <- obj$metadf
  cB <- data.table::as.data.table(cB)


  # Figure out which mutation counts are in the cB
  mutcounts_in_cB <- find_mutcounts(obj)

  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # Which populations to analyze?
  if(length(populations) == 1 && populations == "all"){

    muts_analyze <- mutcounts_in_cB

  }else{

    muts_analyze <- populations

    if(!all(muts_analyze %in% mutcounts_in_cB)){
      stop("You specified mutrate_populations not present in your cB!")
    }


  }

  nucs_analyze <- basecounts_in_cB[which(mutcounts_in_cB %in% muts_analyze)]


  # Infer proportion of each population
  mutest <- vector(mode = "list", length = length(muts_analyze))

  for(i in seq_along(muts_analyze)){

    ### Infer polds from -label data and/or parse pold estimates

    # What columns represent label times of interest?
    if(length(muts_analyze) == 1){

      tl_cols_possible <- c("tl", "tpulse")

    }else{

      tl_cols_possible <- c(paste0("tl_", muts_analyze),
                            paste0("tpulse_", muts_analyze))


    }

    tl_cols <- tl_cols_possible[tl_cols_possible %in% colnames(metadf)]

    # Which samples are -label?
    samples_with_no_label <- metadf %>%
      dplyr::rowwise() %>%
      dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) == 0)) %>%
      dplyr::ungroup() %>%
      dplyr::select(sample) %>%
      unlist() %>%
      unname()

    # Infer background mutation rate from -label data
    if(pold_from_nolabel){

      metadf <- data.table::setDT(data.table::copy(metadf))
      setkey(cB, sample)
      setkey(metadf, sample)

      if(is.null(grouping_factors)){

        background_rates <- cB[metadf, nomatch = NULL][sample %in% samples_with_no_label][,
                                                       .(pold_est = sum(get(muts_analyze[i])*n)/sum(get(nucs_analyze[i])*n))
                                                       ]

        pold_dt <- data.table(pold_est = as.numeric(background_rates$pold_est[1]),
                              sample = metadf$sample[!(metadf$sample %in% samples_with_no_label)])

      }else{

        background_rates <- cB[metadf, nomatch = NULL][sample %in% samples_with_no_label][,
                                                       .(pold_est = sum(get(muts_analyze[i])*n)/sum(get(nucs_analyze[i])*n)),
                                                       by = grouping_factors
        ]

        cols_to_keep <- c("sample", grouping_factors)
        metagroup <- metadf[!(sample %in% samples_with_no_label)][, ..cols_to_keep]
        pold_dt <- background_rates[metagroup, on = grouping_factors, nomatch = NULL][,c(grouping_factors) := rep(NULL, times = length(grouping_factors))]

      }



    }else{

      # Create lookup table if pold_est is provided
      if(length(pold_est) == 1){

        pold_dt <- data.table(pold_est = pold_est,
                              sample = metadf$sample[!(metadf$sample %in% samples_with_no_label)])

      }else if(length(pold_est) > 1){

        samps_provided <- names(pold_est)
        label_samps <- metadf$sample[!(metadf$sample %in% samples_with_no_label)]

        if(!all(samps_provided %in% label_samps)){

          stop("Names of pold_est vector elements include samples not found in the metadf!")

        }

        if(!all(label_samps %in% samps_provided)){

          stop("Not all samples in metadf are found as names of the pold_est vector elements!")

        }

        pold_dt <- data.table::data.table(pold_est = pold_est,
                                          sample = samps_provided)

      }else{

        # Hack to filter out -label samples
        pold_dt <- data.table::data.table(sample = metadf$sample[!(metadf$sample %in% samples_with_no_label)])

      }


    }



    ### Estimate mutation rates

    group_cols <- c("sample", muts_analyze[i], nucs_analyze[i])

    mutest_dt <- cB[,.(n = sum(n)),
                    by = group_cols]


    mutest_temp <- mutest_dt[pold_dt, on = "sample", nomatch = NULL][,.(params = list(fit_tcmm(muts = get(muts_analyze[i]),
                                                         nucs = get(nucs_analyze[i]),
                                                         n = n,
                                                         Poisson = FALSE,
                                                         pnew_prior_mean = pnew_prior_mean,
                                                         pnew_prior_sd = pnew_prior_sd,
                                                         pold_prior_mean = pold_prior_mean,
                                                         pold_prior_sd = pold_prior_sd,
                                                       pold = pold_est))), by = sample]

    if(is.null(pold_est) & !pold_from_nolabel){

      mutest_temp[, c("p1", "p2") := .(inv_logit(sapply(params, `[[`, 2)),
                                       inv_logit(sapply(params, `[[`, 3)))]

      mutest_temp[,c("pold",
                     "pnew") := .(min(c(p1, p2)),
                                  max(c(p1, p2))), by = 1:nrow(mutest_temp)]

      mutest_temp[,c("p1", "p2") := .(NULL, NULL)]

      mutest[[i]] <- mutest_temp

    }else{

      mutest_temp[, c("p1") := .(inv_logit(sapply(params, `[[`, 2)))]

      mutest_temp <- mutest_temp[pold_dt, on = "sample", nomatch = NULL]

      mutest_temp[,c("pold",
                     "pnew") := .(pold_est,
                                  p1), by = 1:nrow(mutest_temp)]

      mutest_temp[,c("p1", "pold_est") := .(NULL, NULL)]

      mutest[[i]] <- mutest_temp

    }



  }

  names(mutest) <- muts_analyze


  obj$mutation_rates <- mutest


  # Add new class information
  if(!("EZbakRMutrates" %in% class(obj))){

    class(obj) <- c("EZbakRMutrates", class(obj))

  }


  return(obj)

}





###################################################
# ARROW DATABASE BACKEND FUNCTIONS
###################################################



#' @describeIn EstimateFractions Mehthod for class **EZbakRArrowData**
#' This is an alternative to the fully in memory **EZbakRData** `EstimateFractions` method that can help
#' with analyses of larger than RAM datasets. The provided "cB" is expected to be
#' an on-disk Arrow Dataset. Furthermore, it is expected to be partitioned by the
#' sample name, which will allow this method to read only a single-sample worth
#' of data into memory at a time. This can significantly reduce RAM usage. Input
#' object should be created with `EZbakRArrowData()`.
#' @import data.table
#' @import arrow
#' @importFrom magrittr %>%
#' @export
EstimateFractions.EZbakRArrowData <- function(obj, features = "all",
                                              mutrate_populations = "all",
                                              fraction_design = NULL,
                                              Poisson = TRUE,
                                              strategy = c("standard", "hierarchical"),
                                              filter_cols = "all",
                                              filter_condition = `&`,
                                              remove_features = c("NA", "__no_feature"),
                                              remove_multi_features = FALSE,
                                              split_multi_features = FALSE,
                                              multi_feature_cols = NULL,
                                              multi_feature_sep = "+",
                                              pnew_prior_mean = -2.94,
                                              pnew_prior_sd = 0.3,
                                              pold_prior_mean = -6.5,
                                              pold_prior_sd = 0.5,
                                              hier_readcutoff = 300,
                                              init_pnew_prior_sd = 0.3,
                                              pnew_prior_sd_min = 0.01,
                                              pnew_prior_sd_max = 0.15,
                                              pold_est = NULL,
                                              pold_from_nolabel = FALSE,
                                              grouping_factors = NULL,
                                              EBayes = TRUE,
                                              character_limit = 20,
                                              overwrite = TRUE){


  #browser()
  ### Hack to deal with devtools::check() NOTEs
  present <- n <- reads <- params <- coverage <- pold <- pnew <- sd <- lpnew_uncert <- pnew_prior <- ..colvect <- NULL
  fit <- mixture_fit <- ci_lower_fn <- ci_upper_fn <- NULL

  rm(..colvect)


  `.` <- list

  ### Check that input is valid

  args <- c(as.list(environment()))

  check_EstimateFractions_input(args)


  strategy <- match.arg(strategy)


  ### Estimate mutation rates
  message("Estimating mutation rates")
  obj <- EstimateMutRates(obj,
                          populations = mutrate_populations,
                          pnew_prior_mean = pnew_prior_mean,
                          pnew_prior_sd = pnew_prior_sd,
                          pold_prior_mean = pold_prior_mean,
                          pold_prior_sd = pold_prior_sd,
                          pold_est = pold_est,
                          pold_from_nolabel = pold_from_nolabel,
                          grouping_factors = grouping_factors)

  mutation_rates <- obj$mutation_rates

  ### Figure out which features will be analyzed

  # cB columns
  cB <- get_arrow_ds(obj)
  obj$cBds <- cB
  cBschema <- arrow::schema(cB)
  cB_cols <- names(cBschema)

  ### Vectors of potential column names

  # Mutation counts possible
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]
  mutcounts[mutcounts == "NA"] <- "NtoA"

  # Which mutcounts are in the cB?
  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]


  ##### BEGINNING OF CODE IDENTICALL TO EZBAKRDATA METHOD (SHOULD REFACTOR
  ##### SO A COMMON FUNCTION IS CALLED HERE BY BOTH)

  # Base count columns in cB
  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # feature columns
  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB,
                                             basecounts_in_cB,
                                             "sample", "n"))]

  # Need to determine which columns of the cB to group reads by
  if(features[1] == "all" & length(features) == 1){

    features_to_analyze <- features_in_cB

  }else{

    if(!all(features %in% features_in_cB)){

      stop("features includes columns that do not exist in your cB!")

    }else{

      features_to_analyze <- features

    }

  }


  ### Figure out which mutational populations to analyze

  if(length(mutrate_populations) == 1 && mutrate_populations == "all"){

    pops_to_analyze <- mutcounts_in_cB

  }else{

    if(!all(mutrate_populations %in% mutcounts_in_cB)){

      stop("mutrate_populations includes columns that do not exist in your cB!")

    }else{

      pops_to_analyze <- mutrate_populations

    }
  }

  # Get nucleotide counts that are needed
  necessary_basecounts <- paste0("n", substr(pops_to_analyze, start = 1, stop = 1))

  ### Create fraction_design data frame if necessary

  if(is.null(fraction_design)){

    fraction_design <- create_fraction_design(pops_to_analyze)

  }


  ### Filter out label-less samples

  metadf <- obj$metadf

  # What columns represent label times of interest?
  if(length(pops_to_analyze) == 1){

    tl_cols_possible <- c("tl", "tpulse")

  }else{

    tl_cols_possible <- c(paste0("tl_", pops_to_analyze),
                          paste0("tpulse_", pops_to_analyze))


  }

  tl_cols <- tl_cols_possible[tl_cols_possible %in% colnames(metadf)]


  ##### END OF CODE IDENTICAL TO EZBAKRDATA

  ### Which samples should be filtered out

  metadf <- obj$metadf

  samples_with_label <- metadf %>%
    dplyr::rowwise() %>%
    dplyr::filter(any(dplyr::c_across(dplyr::all_of(tl_cols)) > 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample) %>%
    unlist() %>%
    unname()

  samples_with_no_label <- metadf %>%
    dplyr::rowwise() %>%
    dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) == 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample) %>%
    unlist() %>%
    unname()



  ### Estimate fraction new for each feature in each sample
  message("Estimating fractions for each sample:")


  # Keep only feature of interest
  if(Poisson){

    cols_to_group <- c(pops_to_analyze, "sample", features_to_analyze)

  }else{

    cols_to_group <- c(necessary_basecounts, pops_to_analyze, "sample", features_to_analyze)

  }


  # Pair down design matrix
  fraction_design <- dplyr::as_tibble(fraction_design)
  mutrate_design <- fraction_design %>%
    dplyr::filter(present) %>%
    dplyr::select(-present)



  # Get read counts for -s4U samples
  all_samples <- metadf$sample
  samples_without_label <- all_samples[!(all_samples %in% samples_with_label)]

  fns <- vector(mode = "list",
                length = length(all_samples))


  col_name <- paste0("logit_fraction_high", pops_to_analyze)
  uncertainty_col <- paste0("se_logit_fraction_high", pops_to_analyze)
  natural_col_name <- paste0("fraction_high", pops_to_analyze)
  ci_lower_col <- paste0("ci_lower_fraction_high", pops_to_analyze)
  ci_upper_col <- paste0("ci_upper_fraction_high", pops_to_analyze)


  if(filter_cols[1] == "all" & length(filter_cols) == 1){

    filter_cols <- features_to_analyze

  }

  if(!all(filter_cols %in% features_to_analyze)){

    stop("Some of your filter_cols are not included in the features you want to analyze!")

  }

  multi_feature_cols_to_use <- multi_feature_cols
  if(remove_multi_features){
    if(is.null(multi_feature_cols) ||
       (length(multi_feature_cols) == 1L && identical(multi_feature_cols, "all"))){
      multi_feature_cols_to_use <- features_to_analyze
    }else{
      if(!all(multi_feature_cols %in% features_to_analyze)){
        stop("multi_feature_cols includes columns that are not in the features being analyzed!")
      }
      multi_feature_cols_to_use <- multi_feature_cols
    }
  }else if(split_multi_features){
    if(!all(multi_feature_cols %in% features_to_analyze)){
      stop("multi_feature_cols includes columns that are not in the features being analyzed!")
    }
  }


  hier_mutrate_name <- NULL
  hier_prior_df <- NULL
  if(ncol(mutrate_design) == 1 && strategy == "hierarchical"){
    hier_mutrate_name <- paste0(paste(features_to_analyze, collapse = "_"), "_", names(obj$mutation_rates)[1])
    obj$mutation_rates[[hier_mutrate_name]] <- NULL
  }

  for(s in seq_along(all_samples)){


    message(paste0("ANALYZING ", all_samples[s], "..."))


    if(all_samples[s] %in% samples_without_label){

      message("Counting reads")

      ctl_cols_to_group <- c("sample", features_to_analyze)

      message()

      sample_fns <- cB %>%
        dplyr::filter(sample == all_samples[s]) %>%
        dplyr::select(dplyr::all_of(c(ctl_cols_to_group, "n"))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ctl_cols_to_group))) %>%
        dplyr::summarise(n = sum(n)) %>%
        dplyr::collect() %>%
        dplyr::ungroup() %>%
        dplyr::mutate(n = as.numeric(n))

      ### Filter out columns not mapping to any feature (easier and faster in data.table)
      sample_fns <- setDT(sample_fns)
      sample_fns <- sample_fns[sample_fns[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ]

      ### Filter out multi-feature mappers if necessary
      if(remove_multi_features){

        message("Filtering out multi-feature mapping reads")

        sample_fns <- remove_multi_feature_rows(sample_fns, multi_feature_cols_to_use, sep = multi_feature_sep)

      }

      sample_fns <- dplyr::as_tibble(sample_fns)

      ### Split multi feature mappers if necessary
      if(split_multi_features){

        message("Splitting multi-feature mapping reads")


        sample_fns <- split_features(sample_fns, multi_feature_cols, ctl_cols_to_group, sep = multi_feature_sep)

      }


      if(ncol(mutrate_design) == 1){
        sample_fns <- dplyr::as_tibble(sample_fns) %>%
          dplyr::mutate(!!natural_col_name := 0,
                        !!col_name := -Inf,
                        !!uncertainty_col := 0,
                        !!ci_lower_col := 0,
                        !!ci_upper_col := 0)
      }else{
        nolabel_mix <- fit_general_mixture(
          dataset = NULL,
          mutrate_design = mutrate_design,
          twocomp = FALSE,
          sample = all_samples[s],
          samples_with_no_label = samples_without_label
        )
        nolabel_tbl <- dplyr::as_tibble(nolabel_mix)[rep(1L, nrow(sample_fns)), , drop = FALSE]
        sample_fns <- dplyr::bind_cols(
          dplyr::as_tibble(sample_fns),
          nolabel_tbl
        )
      }



    }else{



      if(Poisson){

        message("Extracting data for sample of interest and summarize out nucleotide content")


        sample_cB <- cB %>%
          dplyr::filter(sample == all_samples[s]) %>%
          dplyr::select(dplyr::all_of(unique(c(cols_to_group, necessary_basecounts, "n")))) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
          dplyr::summarise(reads = sum(n),
                           dplyr::across(dplyr::all_of(necessary_basecounts),
                                         ~ sum(as.numeric(.x) * as.numeric(n)) / sum(as.numeric(n)))) %>%
          dplyr::collect() %>%
          dplyr::ungroup() %>%
          dplyr::rename(n = reads) %>%
          dplyr::mutate(n = as.numeric(n),
                        dplyr::across(dplyr::all_of(unique(c(pops_to_analyze, necessary_basecounts))), as.numeric))

      }else{

        message("Extracting data for sample of interest")

        sample_cB <- cB %>%
          dplyr::filter(sample == all_samples[s]) %>%
          dplyr::select(dplyr::all_of(unique(c(cols_to_group, "n")))) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
          dplyr::summarise(n = sum(n)) %>%
          dplyr::collect()  %>%
          dplyr::ungroup() %>%
          dplyr::mutate(n = as.numeric(n),
                        dplyr::across(dplyr::all_of(unique(c(pops_to_analyze, necessary_basecounts))), as.numeric))
      }


      ### Filter out columns not mapping to any feature (easier and faster in data.table)
      sample_cB <- setDT(sample_cB)
      sample_cB <- sample_cB[sample_cB[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ]

      ### Filter out multi-feature mappers if necessary
      if(remove_multi_features){

        message("Filtering out multi-feature mapping reads")

        sample_cB <- remove_multi_feature_rows(sample_cB, multi_feature_cols_to_use, sep = multi_feature_sep)

      }

      sample_cB <- dplyr::as_tibble(sample_cB)

      ### Split multi feature mappers if necessary
      if(split_multi_features){

        message("Splitting multi-feature mapping reads")

        sample_cB <- split_features(sample_cB, multi_feature_cols, c(cols_to_group, necessary_basecounts), sep = multi_feature_sep)


      }



      sample_cB <- setDT(sample_cB)
      setkey(sample_cB, sample)


      if(ncol(mutrate_design) == 1){

        ### Can optimize for two-component mixture model

        if(unlist(mutrate_design[1,1])){

          highpop <- 1

        }else{

          highpop <- 2

        }


        ### Add mutation rate info
        message("Adding mutation rate estimation information")

        # Extract mutation rates
        mutrates <- setDT(obj$mutation_rates[[1]] %>%
                            dplyr::select(-params))

        # Key for fast join
        setkey(mutrates, sample)
        setkey(sample_cB, sample)

        # Join
        sample_cB <- sample_cB[mutrates, nomatch = NULL]


        ### Estimate fraction new
        message("Estimating fractions")

        ##### MOST OF THE REST OF THIS FUNCTION IS NEARLY IDENTICAL TO DEFAULT METHOD

        if(strategy == "hierarchical"){

          message("FITTING HIERARCHICAL TWO-COMPONENT MIXTURE MODEL:")

          ### Get coverages to filter by; only want high coverage for feature-specific
          ### pnew estimation

          coverages <- sample_cB[,.(coverage = sum(n)),
                                 by = c("sample", features_to_analyze)]
          coverages <- coverages[coverage > hier_readcutoff][, c("coverage") := .(NULL)]


          ### Get feature-specific pnew estimates

          message("Estimating distribution of feature-specific pnews")

          keyvector <- c("sample", features_to_analyze)
          setkeyv(coverages, keyvector)
          setkeyv(sample_cB, keyvector)

          sample_cB_highcov <- sample_cB[coverages, nomatch = NULL]
          sample_cB_highcov[, pnew_prior_init := logit(pnew)]

          feature_specific <- solve_tcmm_hier_dt(
            dt = sample_cB_highcov,
            features_to_analyze = features_to_analyze,
            pop_col = pops_to_analyze,
            base_col = necessary_basecounts,
            samples_with_no_label = samples_without_label,
            Poisson = Poisson,
            pnew_prior_mean_col = "pnew_prior_init",
            pnew_prior_sd_arg = init_pnew_prior_sd,
            fraction_prior_sd = 1
          )

          feature_specific[, c("pnew", "lpnew_uncert") := .(inv_logit(p2), p2_u)]


          ### Hierarchical model

          message("Estimating fractions with feature-specific pnews")

          global_est <- estimate_hier_pnew_prior(
            feature_specific = feature_specific,
            mutrates = mutrates[sample == all_samples[s]],
            pnew_prior_sd_min = pnew_prior_sd_min,
            pnew_prior_sd_max = pnew_prior_sd_max
          )

          setkey(global_est, sample)
          setkey(sample_cB, sample)

          sample_cB_joined <- sample_cB[global_est[, .(sample, pnew_prior, pnew_prior_sd)], nomatch = NULL]

          feature_specific <- solve_tcmm_hier_dt(
            dt = sample_cB_joined,
            features_to_analyze = features_to_analyze,
            pop_col = pops_to_analyze,
            base_col = necessary_basecounts,
            samples_with_no_label = samples_without_label,
            Poisson = Poisson,
            pnew_prior_mean_col = "pnew_prior",
            pnew_prior_sd_arg = "pnew_prior_sd",
            fraction_prior_sd = 1.7,
            EBayes = EBayes,
            compute_ci = TRUE
          )

          feature_specific[, c(col_name, "pnew",
                               uncertainty_col) := .(p1,
                                                     inv_logit(p2),
                                                     p1_u)]



          sample_fns <- dplyr::as_tibble(feature_specific) %>%
            dplyr::mutate(!!natural_col_name := inv_logit(!!dplyr::sym(col_name)),
                          !!ci_lower_col := ci_lower_fn,
                          !!ci_upper_col := ci_upper_fn) %>%
            dplyr::select(sample, !!features_to_analyze, !!natural_col_name,
                          !!col_name, !!uncertainty_col, !!ci_lower_col, !!ci_upper_col, n)


          ### Save feature-specific pnews

          mutrates[, pnew := NULL]

          colvect <- c("sample", features_to_analyze, "pnew")

          feature_mutrates <- feature_specific[,..colvect][
            mutrates, on = "sample", nomatch = NULL
          ]


          ### NOTE: NOT DOING A GREAT JOB AT HANDLING
          ### MULTIPLE DIFFERENT FEATURE SET HIERARCHICAL
          ### ANALYSES RIGHT NOW
          current_name <- names(obj$mutation_rates)[1]
          hier_mutrate_name <- paste0(paste(features_to_analyze, collapse = "_"), "_", current_name)
          hier_prior_df <- dplyr::bind_rows(
            hier_prior_df,
            dplyr::as_tibble(global_est[, .(sample, pnew_prior, pnew_prior_sd)])
          )

          # Kinda cutesy, but works because bind_rows(dataframe, NULL) = dataframe
          # so don't need a special condition for s == 1.
          obj$mutation_rates[[hier_mutrate_name]] <- dplyr::bind_rows(
            obj$mutation_rates[[hier_mutrate_name]],
            feature_mutrates
          )


        }else{

          sample_fns <- solve_tcmm_standard_dt(
            dt = sample_cB,
            features_to_analyze = features_to_analyze,
            pop_col = pops_to_analyze,
            base_col = necessary_basecounts,
            samples_with_no_label = samples_without_label,
            Poisson = Poisson,
            col_name = col_name,
            uncertainty_col = uncertainty_col,
            natural_col_name = natural_col_name,
            EBayes = EBayes
          )


        }


      }else{

        pops_to_analyze <- colnames(mutrate_design)
        necessary_basecounts <- paste0("n", substr(pops_to_analyze, start = 1, stop = 1))

        sample_fns <- solve_general_mixture_dt(
          dt = sample_cB,
          features_to_analyze = features_to_analyze,
          mutrate_design = mutrate_design,
          mutcols = pops_to_analyze,
          basecols = necessary_basecounts,
          Poisson = Poisson,
          samples_with_no_label = samples_without_label,
          mutation_rates = mutation_rates
        )

      }



    }




    fns[[s]] <- sample_fns

  }


  fns <- dplyr::bind_rows(fns)

  if(remove_multi_features && nrow(fns) == 0L){
    stop("No reads remain after filtering with remove_multi_features = TRUE!")
  }


  message("Processing output")



  ##### REST OF THIS CODE IS IDENTICAL TO AS IN DEFAULT METHOD


  # What should output be named?
  fraction_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")

  if(nchar(fraction_vect) > character_limit){

    num_fractions <- length(obj[['fractions']])
    fraction_vect <- paste0("fractions", num_fractions + 1)

  }

  # Are there any metadata or fractions objects at this point?
  if(length(obj[['metadata']][['fractions']]) > 0){

    fraction_vect <- decide_output(obj,
                                   proposed_name = fraction_vect,
                                   type = "fractions",
                                   features = features_to_analyze,
                                   populations = pops_to_analyze,
                                   fraction_design = fraction_design,
                                   overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'fractions',
                               features = features_to_analyze,
                               populations = pops_to_analyze,
                               fraction_design = fraction_design,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }


  # Save
  obj[['fractions']][[fraction_vect]] <- dplyr::as_tibble(fns)

  # Save metadata
  obj[['metadata']][['fractions']][[fraction_vect]] <- list(features = features_to_analyze,
                                                            populations = pops_to_analyze,
                                                            fraction_design = fraction_design,
                                                            strategy = if(ncol(mutrate_design) == 1) strategy else "general",
                                                            Poisson = Poisson,
                                                            remove_multi_features = remove_multi_features,
                                                            split_multi_features = split_multi_features,
                                                            multi_feature_cols = multi_feature_cols_to_use,
                                                            multi_feature_sep = multi_feature_sep,
                                                            filter_cols = filter_cols,
                                                            remove_features = remove_features,
                                                            hier_mutrate_table = hier_mutrate_name,
                                                            hier_prior = hier_prior_df,
                                                            repeatID = repeatID)

  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)



}



#' @describeIn EstimateMutRates Method for class **EZbakRArrowData**
#' Estimate mutation rates using a partially on-disk object.
#' @import data.table
#' @export
EstimateMutRates.EZbakRArrowData <- function(obj,
                                             populations = "all",
                                             pnew_prior_mean = -2.94,
                                             pnew_prior_sd = 0.3,
                                             pold_prior_mean = -6.5,
                                             pold_prior_sd = 0.5,
                                             pold_est = NULL,
                                             pold_from_nolabel = FALSE,
                                             grouping_factors = NULL
){


  #browser()
  ### Hack to deal with devtools::check() NOTEs
  n <- params <- p1 <- p2 <- ..cols_to_keep <- muts <- nucs <- NULL

  rm(..cols_to_keep)

  `.` <- list

  ### Figure out which features will be analyzed

  # cB columns
  cB <- get_arrow_ds(obj)
  obj$cBds <- cB
  cBschema <- arrow::schema(cB)
  cB_cols <- names(cBschema)

  ### Vectors of potential column names

  # Mutation counts possible
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]
  mutcounts[mutcounts == "NA"] <- "NtoA"

  # Which mutcounts are in the cB?
  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]

  # Base count columns in cB
  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))


  # Which populations to analyze?
  if(length(populations) == 1 && populations == "all"){

    muts_analyze <- mutcounts_in_cB

  }else{

    muts_analyze <- populations

    if(!all(muts_analyze %in% mutcounts_in_cB)){
      stop("You specified mutrate_populations not present in your cB!")
    }


  }

  nucs_analyze <- basecounts_in_cB[which(mutcounts_in_cB %in% muts_analyze)]


  # Infer proportion of each population
  mutest <- vector(mode = "list", length = length(muts_analyze))

  ### Which samples should be filtered out?

  metadf <- obj$metadf

  # What columns represent label times of interest?
  if(length(muts_analyze) == 1){

    tl_cols_possible <- c("tl", "tpulse")

  }else{

    tl_cols_possible <- c(paste0("tl_", muts_analyze),
                          paste0("tpulse_", muts_analyze))


  }

  tl_cols <- tl_cols_possible[tl_cols_possible %in% colnames(metadf)]


  samples_with_label <- metadf %>%
    dplyr::rowwise() %>%
    dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) > 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample) %>%
    unlist() %>%
    unname()


  ### STEPS:
  # 1) Extract sample's data and summarize
  # 2) Collect and pass to model as in default method
  for(i in seq_along(muts_analyze)){


    # Infer proportion of each population
    sample_mutest <- dplyr::tibble()

    group_cols <- c("sample", muts_analyze[i], nucs_analyze[i])

    ### Infer polds from -label data and/or parse pold estimates

    # Which samples are -label?
    samples_with_no_label <- metadf$sample[!(metadf$sample %in% samples_with_label)]

    if(pold_from_nolabel){

      ctl_cB <- dplyr::tibble()
      ctl_group_cols <- c("")

      if(length(samples_with_no_label) == 0){

        stop("Cannot set pold_from_nolabel to TRUE if you have no -label controls!")

      }

      ### Compile all -label data
      ctl_cols_sel <- unique(c("sample", "n", muts_analyze[i], nucs_analyze[i]))
      for(c in seq_along(samples_with_no_label)){


        ctl_cB <- cB %>%
          dplyr::filter(sample == samples_with_no_label[c]) %>%
          dplyr::select(dplyr::all_of(ctl_cols_sel)) %>%
          dplyr::group_by(sample) %>%
          dplyr::summarise(muts = sum(as.numeric(n) * as.numeric(!!dplyr::sym(muts_analyze[i]))),
                           nucs = sum(as.numeric(n) * as.numeric(!!dplyr::sym(nucs_analyze[i])))) %>%
          dplyr::collect() %>%
          dplyr::mutate(muts = as.numeric(muts),
                        nucs = as.numeric(nucs)) %>%
          dplyr::bind_rows(ctl_cB)


      }


      ### Infer background rates
      metadf <- data.table::setDT(data.table::copy(metadf))
      ctl_cB <- data.table::setDT(ctl_cB)
      setkey(ctl_cB, sample)
      setkey(metadf, sample)


      if(is.null(grouping_factors)){

        background_rates <- ctl_cB[metadf, nomatch = NULL][sample %in% samples_with_no_label][,
                                                       .(pold_est = sum(muts)/sum(nucs))
        ]

        pold_dt <- data.table(pold_est = as.numeric(background_rates$pold_est[1]),
                              sample = samples_with_label)

      }else{

        background_rates <- ctl_cB[metadf, nomatch = NULL][sample %in% samples_with_no_label][,
                                                       .(pold_est = sum(muts)/sum(nucs)),
                                                       by = grouping_factors
        ]

        cols_to_keep <- c("sample", grouping_factors)
        metagroup <- metadf[sample %in% samples_with_label][,..cols_to_keep]
        pold_dt <- background_rates[metagroup, on = grouping_factors, nomatch = NULL][,c(grouping_factors) := rep(NULL, times = length(grouping_factors))]

      }


    }else{

      # Create lookup table if pold_est is provided
      if(length(pold_est) == 1){

        pold_dt <- data.table(pold_est = pold_est,
                              sample = samples_with_label)

      }else if(length(pold_est) > 1){

        samps_provided <- names(pold_est)

        if(!all(samps_provided %in% samples_with_label)){

          stop("Names of pold_est vector elements include samples not found in the metadf!")

        }

        if(!all(samples_with_label %in% samps_provided)){

          stop("Not all samples in metadf are found as names of the pold_est vector elements!")

        }

        pold_dt <- data.table::data.table(pold_est = pold_est,
                                          sample = samps_provided)

      }else{

        # Dummy table because sample_cB will get merged with this one way or another
        pold_dt <- data.table::data.table(sample = samples_with_label)

      }

    }


    ### Estimate mutation rates

    for(s in seq_along(samples_with_label)){

      sample_cB <- cB %>%
        dplyr::filter(sample == samples_with_label[s]) %>%
        dplyr::select(dplyr::all_of(unique(c(group_cols, "n")))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(n = sum(n)) %>%
        dplyr::collect() %>%
        dplyr::ungroup() %>%
        dplyr::mutate(n = as.numeric(n),
                      dplyr::across(dplyr::all_of(c(muts_analyze[i], nucs_analyze[i])), as.numeric))


      mutest_dt <- data.table::setDT(sample_cB)


      mutest_temp <- mutest_dt[pold_dt, on = "sample", nomatch = NULL][,.(params = list(fit_tcmm(muts = get(muts_analyze[i]),
                                                         nucs = get(nucs_analyze[i]),
                                                         n = n,
                                                         Poisson = FALSE,
                                                         pnew_prior_mean = pnew_prior_mean,
                                                         pnew_prior_sd = pnew_prior_sd,
                                                         pold_prior_mean = pold_prior_mean,
                                                         pold_prior_sd = pold_prior_sd,
                                                         pold = pold_est))), by = sample]


      if(is.null(pold_est) & !pold_from_nolabel){

        mutest_temp[, c("p1", "p2") := .(inv_logit(sapply(params, `[[`, 2)),
                                         inv_logit(sapply(params, `[[`, 3)))]

        mutest_temp[,c("pold",
                       "pnew") := .(min(c(p1, p2)),
                                    max(c(p1, p2))), by = 1:nrow(mutest_temp)]

        mutest_temp[,c("p1", "p2") := .(NULL, NULL)]


      }else{

        mutest_temp[, c("p1") := .(inv_logit(sapply(params, `[[`, 2)))]

        mutest_temp <- mutest_temp[pold_dt, on = "sample", nomatch = NULL]

        mutest_temp[,c("pold",
                       "pnew") := .(pold_est,
                                    p1), by = 1:nrow(mutest_temp)]

        mutest_temp[,c("p1", "pold_est") := .(NULL, NULL)]

      }

      sample_mutest <- dplyr::bind_rows(sample_mutest, mutest_temp)

    }

    mutest[[i]] <- sample_mutest

  }


  names(mutest) <- muts_analyze


  obj$mutation_rates <- mutest


  # Add new class information
  if(!("EZbakRMutrates" %in% class(obj))){

    class(obj) <- c("EZbakRMutrates", class(obj))

  }


  return(obj)

}





###################################################
# MISCELLANEOUS HELPER FUNCTIONS
###################################################


# Filter out rows containing multi-feature assignments (e.g., "featureA+featureB")
remove_multi_feature_rows <- function(dt, cols_check, sep = "+"){

  if(nrow(dt) == 0L || length(cols_check) == 0L){
    return(dt)
  }

  has_multi <- rep(FALSE, nrow(dt))

  for(col in cols_check){
    vals <- dt[[col]]
    if(is.factor(vals)){
      lvls <- levels(vals)
      bad_lvl <- grepl(sep, lvls, fixed = TRUE)
      if(any(bad_lvl)){
        col_multi <- bad_lvl[as.integer(vals)]
        col_multi[is.na(col_multi)] <- FALSE
        has_multi <- has_multi | col_multi
      }
    }else{
      if(!is.character(vals)){
        vals <- as.character(vals)
      }
      col_multi <- grepl(sep, vals, fixed = TRUE)
      col_multi[is.na(col_multi)] <- FALSE
      has_multi <- has_multi | col_multi
    }
  }

  if(any(has_multi)){
    dt <- dt[!has_multi, ]
  }

  return(dt)

}


# Split up multi-feature sets
split_features <- function(cB, multi_feature_cols, grouping_cols, sep = "+"){

  ### Hack to deal with devtools::check() NOTEs
  n <- NULL

  `.` <- list

  if(nrow(cB) == 0L){
    return(data.table::as.data.table(cB))
  }

  sep_regex <- gsub("([.\\\\|()\\[\\]{}^$*+?])", "\\\\\\1", sep, perl = TRUE)

  ### Copy junction data more efficiently
  unique_juncs <- cB %>%
    dplyr::select(!!multi_feature_cols) %>%
    dplyr::distinct()


  ### have to move in and out of the tidyverse because this next line is much, much easier in data.table
  unique_juncs <- setDT(unique_juncs)
  new_temp_cols <- paste0("unrolled_", 1:length(multi_feature_cols))
  unique_juncs <- unique_juncs[,(new_temp_cols) := mget(multi_feature_cols)]


  ### Back to the tidyverse, because can't easily do this in data.table (sigh...)
  unique_juncs <- unique_juncs %>%
    tidyr::separate_rows(dplyr::all_of(new_temp_cols), sep = sep_regex)

  cB <- cB %>%
    dplyr::inner_join(unique_juncs, by = multi_feature_cols,
               relationship = "many-to-many")


  ### Back to data.table again because its just easier and more efficient
  setDT(cB)

  cB <- cB[,.(n = sum(n)), by = c(grouping_cols, new_temp_cols)][
    ,(multi_feature_cols) := mget(new_temp_cols)
  ][,(new_temp_cols) := NULL]


  return(cB)

}




# Find the mutation count columns in the cB
find_mutcounts <- function(obj){

  `.` <- list

  # Mutation count column names
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]
  mutcounts[mutcounts == "NA"] <- "NtoA"



  # Which mutcounts are in the cB?
  cB_cols <- colnames(obj$cB)


  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]

  return(mutcounts_in_cB)



}

# Softmax function
softmax <- function(vect){

  return(exp(vect)/sum(exp(vect)))

}



#################################################
# TWO COMPONENT MIXTURE MODELING HELPER FUNCTIONS
#################################################


log_add2 <- function(a, b){
  ## recycle to common length (base R handles this automatically in pmax)
  m   <- pmax(a, b)

  ## core computation
  out <- m + log1p(exp(-abs(a - b)))   # safe for |a-b| >> 0

  ## both a and b are -Inf  ->  log(0)
  out[is.infinite(m) & m < 0] <- -Inf
  out

}


# Two-component mixture model likelihood
tcmml <- function(param, muts, nucs, n,
                  pnew = NULL,
                  pold = NULL,
                  pnew_prior_mean = -2.94,
                  pnew_prior_sd = 0.3,
                  pold_prior_mean = -6.5,
                  pold_prior_sd = 0.5,
                  fraction_prior_mean = 0,
                  fraction_prior_sd = 1.7,
                  Poisson = TRUE){


  ### Parse relevant parameters to make later code cleaner

  fn <- inv_logit(param[1])

  count <- 2
  estimate_pnew <- FALSE
  estimate_pold <- FALSE
  if(is.null(pnew)){

    pnew <- inv_logit(param[count])
    count <- count + 1
    estimate_pnew <- TRUE

  }


  if(is.null(pold)){

    pold <- inv_logit(param[count])
    estimate_pold <- TRUE

  }


  ### Likelihoods

  if(Poisson){

    ll <- log_add2(log(fn) + stats::dpois(muts, nucs*pnew, log = TRUE),
                   log(1 - fn) + stats::dpois(muts, nucs*pold, log = TRUE))

  }else{

    ll <- log_add2(log(fn) + stats::dbinom(muts, nucs, pnew, log = TRUE),
                   log(1 - fn) + stats::dbinom(muts, nucs, pold, log = TRUE))

  }


  ### Add prior weights if necessary

  if(estimate_pnew & estimate_pold){

    if(pold > pnew){

      prior <- stats::dnorm(param[3], pnew_prior_mean, pnew_prior_sd,
                        log = TRUE) +
        stats::dnorm(param[2], pold_prior_mean, pold_prior_sd,
                        log = TRUE) +
        stats::dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
                        log = TRUE)

    }else{

      prior <- stats::dnorm(param[3], pold_prior_mean, pold_prior_sd,
                        log = TRUE) +
        stats::dnorm(param[2], pnew_prior_mean, pnew_prior_sd,
                        log = TRUE) +
        stats::dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
              log = TRUE)


    }


  }else if(estimate_pnew){

    prior <- stats::dnorm(param[2], pnew_prior_mean, pnew_prior_sd,
                        log = TRUE) +
      stats::dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
            log = TRUE)


  }else if(estimate_pold){

    prior <- stats::dnorm(param[2], pold_prior_mean, pold_prior_sd,
                    log = TRUE) +
      stats::dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
            log = TRUE)

  }else{

    prior <- stats::dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
                   log = TRUE)

  }


  # if(!is.finite(sum(n*log(ll)) + prior)){
  #
  #   browser()
  #
  # }


  ll <- -(sum(n*ll) + prior)



  return(ll)


}






# Fit all two-componenet mixture model variations
fit_tcmm <- function(muts, nucs, n, pnew = NULL, pold = NULL,
                     pnew_prior_mean = -2.94,
                     pnew_prior_sd = 0.3,
                     pold_prior_mean = -6.5,
                     pold_prior_sd = 0.5,
                     Poisson = TRUE,
                     fraction_prior_mean = 0,
                     fraction_prior_sd = 1.7,
                     return_fit = FALSE){


  init <- 0
  upper_bound <- 9
  lower_bound <- -9

  estimate_pnew <- FALSE
  estimate_pold <- FALSE
  if(is.null(pnew)){

    init <- append(init, -2)
    upper_bound <- append(upper_bound, 0)
    lower_bound <- append(lower_bound, -9)
    estimate_pnew <- TRUE

  }

  if(is.null(pold)){

    init <- append(init, -6.5)
    upper_bound <- append(upper_bound, 0)
    lower_bound <- append(lower_bound, -9)
    estimate_pold <- TRUE

  }

  muts <- as.numeric(muts)
  nucs <- as.numeric(nucs)
  n <- as.numeric(n)

  if(Poisson){
    c_const <- ifelse(muts == 0 & nucs == 0, 0,
                      muts * log(pmax(nucs, 1e-300)) - lgamma(muts + 1))
  }else{
    c_const <- lchoose(nucs, muts)
  }

  if(estimate_pnew && !estimate_pold){

    ll_old <- if(Poisson){
      stats::dpois(muts, nucs * pold, log = TRUE)
    }else{
      stats::dbinom(muts, nucs, pold, log = TRUE)
    }
    N_tot <- sum(n)
    inv_s1 <- 1 / (fraction_prior_sd * fraction_prior_sd)
    inv_s2 <- 1 / (pnew_prior_sd * pnew_prior_sd)
    prior_const <- -log(2 * pi) - log(fraction_prior_sd) - log(pnew_prior_sd)

    fn_2d <- function(param){
      e1 <- param[1]; e2 <- param[2]
      pn <- stats::plogis(e2)
      if(Poisson){
        log_pn <- if(e2 > 0) -log1p(exp(-e2)) else e2 - log1p(exp(e2))
        ll_new <- c_const + muts * log_pn - nucs * pn
      }else{
        sp2 <- if(e2 > 0) e2 + log1p(exp(-e2)) else log1p(exp(e2))
        ll_new <- c_const + muts * e2 - nucs * sp2
      }
      lt1 <- if(e1 > 0) -log1p(exp(-e1)) else e1 - log1p(exp(e1))
      lt0 <- if(e1 > 0) -e1 - log1p(exp(-e1)) else -log1p(exp(e1))
      ll <- log_add2(lt1 + ll_new, lt0 + ll_old)
      prior <- prior_const - 0.5 * (e1 - fraction_prior_mean)^2 * inv_s1 -
        0.5 * (e2 - pnew_prior_mean)^2 * inv_s2
      -(sum(n * ll) + prior)
    }

    gr_2d <- function(param){
      e1 <- param[1]; e2 <- param[2]
      theta <- stats::plogis(e1)
      pn <- stats::plogis(e2)
      if(Poisson){
        log_pn <- if(e2 > 0) -log1p(exp(-e2)) else e2 - log1p(exp(e2))
        ll_new <- c_const + muts * log_pn - nucs * pn
        d1 <- (1 - pn) * (muts - nucs * pn)
      }else{
        sp2 <- if(e2 > 0) e2 + log1p(exp(-e2)) else log1p(exp(e2))
        ll_new <- c_const + muts * e2 - nucs * sp2
        d1 <- muts - nucs * pn
      }
      delta <- ll_new - ll_old
      delta[is.nan(delta)] <- 0
      w <- stats::plogis(e1 + delta)
      nw <- n * w
      g1 <- N_tot * theta - sum(nw) + (e1 - fraction_prior_mean) * inv_s1
      g2 <- -sum(nw * d1) + (e2 - pnew_prior_mean) * inv_s2
      c(g1, g2)
    }

    fit <- stats::optim(par = init, fn = fn_2d, gr = gr_2d,
                        hessian = FALSE,
                        method = "L-BFGS-B", lower = lower_bound, upper = upper_bound)

    e1 <- fit$par[1]; e2 <- fit$par[2]
    theta <- stats::plogis(e1); pn <- stats::plogis(e2)
    if(Poisson){
      log_pn <- if(e2 > 0) -log1p(exp(-e2)) else e2 - log1p(exp(e2))
      ll_new <- c_const + muts * log_pn - nucs * pn
      d1 <- (1 - pn) * (muts - nucs * pn)
      d2 <- -pn * (1 - pn) * (muts + nucs * (1 - 2 * pn))
    }else{
      sp2 <- if(e2 > 0) e2 + log1p(exp(-e2)) else log1p(exp(e2))
      ll_new <- c_const + muts * e2 - nucs * sp2
      d1 <- muts - nucs * pn
      d2 <- -nucs * pn * (1 - pn)
    }
    delta <- ll_new - ll_old
    delta[is.nan(delta)] <- 0
    w <- stats::plogis(e1 + delta)
    nw <- n * w
    nww <- nw * (1 - w)
    h11 <- N_tot * theta * (1 - theta) - sum(nww) + inv_s1
    h12 <- -sum(nww * d1)
    h22 <- -sum(nww * d1 * d1 + nw * d2) + inv_s2
    fit$hessian <- matrix(c(h11, h12, h12, h22), nrow = 2, ncol = 2)

  }else if(estimate_pnew && estimate_pold){

    fn_3d <- function(param){
      e1 <- param[1]; e2 <- param[2]; e3 <- param[3]
      pn <- stats::plogis(e2); po <- stats::plogis(e3)
      if(Poisson){
        log_pn <- if(e2 > 0) -log1p(exp(-e2)) else e2 - log1p(exp(e2))
        log_po <- if(e3 > 0) -log1p(exp(-e3)) else e3 - log1p(exp(e3))
        ll_new <- c_const + muts * log_pn - nucs * pn
        ll_old <- c_const + muts * log_po - nucs * po
      }else{
        sp2 <- if(e2 > 0) e2 + log1p(exp(-e2)) else log1p(exp(e2))
        sp3 <- if(e3 > 0) e3 + log1p(exp(-e3)) else log1p(exp(e3))
        ll_new <- c_const + muts * e2 - nucs * sp2
        ll_old <- c_const + muts * e3 - nucs * sp3
      }
      lt1 <- if(e1 > 0) -log1p(exp(-e1)) else e1 - log1p(exp(e1))
      lt0 <- if(e1 > 0) -e1 - log1p(exp(-e1)) else -log1p(exp(e1))
      ll <- log_add2(lt1 + ll_new, lt0 + ll_old)
      if(po > pn){
        prior <- stats::dnorm(e3, pnew_prior_mean, pnew_prior_sd, log = TRUE) +
          stats::dnorm(e2, pold_prior_mean, pold_prior_sd, log = TRUE) +
          stats::dnorm(e1, fraction_prior_mean, fraction_prior_sd, log = TRUE)
      }else{
        prior <- stats::dnorm(e3, pold_prior_mean, pold_prior_sd, log = TRUE) +
          stats::dnorm(e2, pnew_prior_mean, pnew_prior_sd, log = TRUE) +
          stats::dnorm(e1, fraction_prior_mean, fraction_prior_sd, log = TRUE)
      }
      -(sum(n * ll) + prior)
    }

    fit <- stats::optim(par = init, fn = fn_3d,
                        hessian = TRUE,
                        method = "L-BFGS-B", lower = lower_bound, upper = upper_bound)

  }else{

    fit <- stats::optim(par = init, fn = tcmml,
                        muts = muts,
                        nucs = nucs,
                        n = n,
                        pnew = pnew,
                        pold = pold,
                        pnew_prior_mean = pnew_prior_mean,
                        pnew_prior_sd = pnew_prior_sd,
                        pold_prior_mean = pold_prior_mean,
                        pold_prior_sd = pold_prior_sd,
                        fraction_prior_mean = fraction_prior_mean,
                        fraction_prior_sd = fraction_prior_sd,
                        Poisson = Poisson,
                        hessian = TRUE,
                        method = "L-BFGS-B", lower = lower_bound, upper = upper_bound)

  }


  if(return_fit){

    return(fit)

  }else{

    uncertainty <- sqrt(abs(diag(solve(fit$hessian))))

    if(estimate_pnew & estimate_pold){

      return(list(p1 = fit$par[1], p2 = fit$par[2], p3 = fit$par[3],
                  p1_u = uncertainty[1], p2_u = uncertainty[2], p3_u = uncertainty[3]))

    }else if(estimate_pnew | estimate_pold){

      return(list(p1 = fit$par[1], p2 = fit$par[2],
                  p1_u = uncertainty[1], p2_u = uncertainty[2]))

    }else{

      return(list(p1 = fit$par[1], p1_u = uncertainty[1]))

    }

  }


}





###################################################
# HELPER FUNCTIONS FOR MULTI LABEL MIXTURE MODELING
###################################################




# Two-component mixture model; optimized for simplest use-case
two_comp_likelihood <- function(param, muts, nucs, Poisson = TRUE,
                                pnew, pold, n){

  if(Poisson){

    ll <- log_add2(
      log(inv_logit(param[1])) + stats::dpois(muts, lambda = pnew*nucs, log = TRUE ),
      log(1 - inv_logit(param[1])) + stats::dpois(muts, lambda = pold*nucs, log = TRUE )
    )

  }else{

    ll <- log_add2(
      log(inv_logit(param[1])) + stats::dbinom(x = muts,size = nucs, prob = pnew, log = TRUE),
      log(1 - inv_logit(param[1])) + stats::dbinom(x = muts, size = nucs, prob = pold, log = TRUE)
    )


  }

  return(-sum(n*ll) - stats::dnorm(param[1], log = TRUE))


}

# Empirical Bayes estimation of sample-specific logit(pnew) prior mean and SD
estimate_hier_pnew_prior <- function(feature_specific,
                                     mutrates,
                                     pnew_prior_sd_min,
                                     pnew_prior_sd_max){

  pnew <- pold <- lpnew_uncert <- pnew_prior <- pnew_prior_sd <- NULL
  `.` <- list

  prior_base <- data.table::as.data.table(mutrates)[, .(
    sample = sample,
    pnew_prior = logit(pnew),
    pold = pold
  )]

  if(nrow(feature_specific) == 0){
    prior_base[, pnew_prior_sd := pnew_prior_sd_min]
    return(prior_base)
  }

  fs_dt <- data.table::as.data.table(feature_specific)
  data.table::setkey(fs_dt, sample)
  data.table::setkey(prior_base, sample)

  fs_joined <- prior_base[fs_dt, nomatch = NULL]
  fs_filt <- fs_joined[is.finite(logit(pnew)) & pnew < 0.3 & pnew > pold]

  if(nrow(fs_filt) > 0){
    sd_dt <- fs_filt[, .(
      pnew_prior_sd = {
        if(.N >= 2){
          emp_var <- stats::var(logit(pnew), na.rm = TRUE) - mean(lpnew_uncert^2, na.rm = TRUE)
          if(is.na(emp_var) || emp_var <= 0) pnew_prior_sd_min else sqrt(emp_var)
        }else{
          pnew_prior_sd_min
        }
      }
    ), by = sample]

    data.table::setkey(sd_dt, sample)
    global_est <- sd_dt[prior_base]
  }else{
    global_est <- data.table::copy(prior_base)
    global_est[, pnew_prior_sd := pnew_prior_sd_min]
  }

  global_est[is.na(pnew_prior_sd), pnew_prior_sd := pnew_prior_sd_min]
  global_est[, pnew_prior_sd := pmin(pmax(pnew_prior_sd, pnew_prior_sd_min), pnew_prior_sd_max)]

  return(global_est[, .(sample, pnew_prior, pnew_prior_sd, pold)])
}


# Abstract the concept of an NR-seq mixture model to oblivion
generalized_likelihood <- function(param, dataset, Poisson = TRUE,
                                   mutrate_design, pnew, pold, highpop,
                                   mutcols, basecols, twocomp = TRUE,
                                   ref_idx = NULL){

  K <- nrow(mutrate_design)
  n_obs <- length(dataset[[1]])

  if(twocomp){

    likelihoods <- rep(0, times = n_obs)

    for(i in 1:K){

      L <- 1
      for(j in 1:ncol(mutrate_design)){
        rate_ij <- as.numeric(mutrate_design[i, j])*pnew[j] +
          pold[j]*(1 - as.numeric(mutrate_design[i, j]))
        if(Poisson){
          L <- L*stats::dpois(dataset[[mutcols[j]]],
                              lambda = rate_ij*dataset[[basecols[j]]])
        }else{
          L <- L*stats::dbinom(dataset[[mutcols[j]]],
                               dataset[[basecols[j]]],
                               prob = rate_ij)
        }
      }

      if(i == highpop){
        likelihoods <- likelihoods + L*inv_logit(param[1])
      }else{
        likelihoods <- likelihoods + L*(1 - inv_logit(param[1]))
      }

    }

    return(-sum(dataset[["n"]]*log(likelihoods)) - stats::dnorm(param[1], log = TRUE))

  }else{

    if(length(param) == K - 1){
      if(is.null(ref_idx)){
        ref_idx <- K
      }
      full_param <- numeric(K)
      full_param[-ref_idx] <- param
      full_param[ref_idx] <- 0
    }else{
      full_param <- param
    }

    m_par <- max(full_param)
    log_denom <- m_par + log(sum(exp(full_param - m_par)))
    log_fractions <- full_param - log_denom

    # Matrix of log component likelihoods (n_obs x K) for log-sum-exp stability
    comp_logL <- matrix(0, nrow = n_obs, ncol = K)

    for(i in 1:K){
      logL_i <- rep(log_fractions[i], times = n_obs)
      for(j in 1:ncol(mutrate_design)){
        rate_ij <- as.numeric(mutrate_design[i, j])*pnew[j] +
          pold[j]*(1 - as.numeric(mutrate_design[i, j]))
        if(Poisson){
          logL_i <- logL_i + stats::dpois(dataset[[mutcols[j]]],
                                          lambda = rate_ij*dataset[[basecols[j]]],
                                          log = TRUE)
        }else{
          logL_i <- logL_i + stats::dbinom(dataset[[mutcols[j]]],
                                           dataset[[basecols[j]]],
                                           prob = rate_ij,
                                           log = TRUE)
        }
      }
      comp_logL[, i] <- logL_i
    }

    row_max <- apply(comp_logL, 1, max)
    log_lik_obs <- row_max + log(rowSums(exp(comp_logL - row_max)))

    # Symmetric Dirichlet(0.5, ..., 0.5) prior in ALR space prevents boundary runaway
    # while remaining invariant to choice of reference component
    prior_ll <- 0.5 * sum(log_fractions)

    return(-sum(dataset[["n"]]*log_lik_obs) - prior_ll)

  }

}

# Softmax function (log-sum-exp stabilized)
softmax <- function(vect){

  m <- max(vect)
  if(!is.finite(m)){
    return(rep(1 / length(vect), times = length(vect)))
  }
  ev <- exp(vect - m)
  return(ev / sum(ev))

}

# Fit and process the output of the generalized likelihood model
fit_general_mixture <- function(dataset, Poisson = TRUE, mutrate_design, twocomp = FALSE,
                                pnew, pold, mutcols, basecols, highpop = NULL,
                                samples_with_no_label = NULL,
                                sample = NULL){

  if(!is.null(sample) && (sample %in% samples_with_no_label)){

    # Figure out what to call each of the fractions
    # Probably "p""muttype""old/new"_"muttype""old/new"
    muttypes <- colnames(mutrate_design)


    # Get vector of population statuses (new and old) for naming the output list
    population_list <- vector(mode = "list", length = nrow(mutrate_design))
    fvector <- rep(0, times = nrow(mutrate_design))
    ismvector <- rep(0, times = nrow(mutrate_design))
    for(i in 1:nrow(mutrate_design)){

      population_vector <- rep("", times = ncol(mutrate_design))
      count <- 1
      for(j in 1:ncol(mutrate_design)){

        if(as.logical(mutrate_design[i, j])){

          population_vector[count] <- "high"

        }else{

          population_vector[count] <- "low"

        }

        count <- count + 1


      }

      if(all(population_vector == "low")){
        fvector[i] <- 1
        ismvector[i] <- Inf
      }else{
        fvector[i] <- 0
        ismvector[i] <- -Inf
      }

      population_list[[i]] <- population_vector

    }


    # Make names for the output list that are easily interpretable (I hope)
    listnames <- lapply(population_list, function(x) paste0(x, muttypes))
    listnames <- sapply(listnames, function(x) paste(x, collapse = "_"),
                        simplify = "vector")

    listnames <- paste0(rep(c("invsoftmax_",
                              "fraction_",
                              "se_invsoftmax_",
                              "ci_lower_fraction_",
                              "ci_upper_fraction_"),
                            each = nrow(mutrate_design)),
                        listnames)

    outlist <- c(as.list(ismvector),
                 as.list(fvector),
                 as.list(rep(0, times = nrow(mutrate_design))),
                 as.list(fvector),
                 as.list(fvector))



    names(outlist) <- listnames

  }else{

    pnew <- pnew[match(colnames(mutrate_design), names(pnew))]
    pold <- pold[match(colnames(mutrate_design), names(pold))]


    if(twocomp){

      tryCatch(
        {
          fit <- stats::optim(par = 0,
                              fn = two_comp_likelihood,
                              dataset = dataset,
                              Poisson = Poisson,
                              pnew = pnew,
                              pold = pold,
                              mutcols = mutcols,
                              basecols = basecols,
                              upper = 7,
                              lower = -7,
                              method = "L-BFGS-B",
                              hessian = TRUE)

        },
        error = function(e){
          message(paste0(c(paste0("Fitting of the generalized mixture model failed for following data: "), utils::capture.output(dataset)), collapse = "\n"))
          stop("Error message from optim: ", e)
        }
      )


      outlist <- list(fit$par, logit(1 - inv_logit(fit$par)))

      names(outlist) <- c(paste0(mutcols, "high"),
                          paste0(mutcols, "low"))


    }else{

      K <- nrow(mutrate_design)
      low_rows <- which(rowSums(as.matrix(mutrate_design)) == 0)
      ref_idx <- if(length(low_rows) > 0) low_rows[1] else K

      tryCatch(
        {
          fit <- stats::optim(par = rep(0, times = K - 1),
                              fn = generalized_likelihood,
                              dataset = dataset,
                              Poisson = Poisson,
                              mutrate_design = mutrate_design,
                              pnew = pnew,
                              pold = pold,
                              twocomp = twocomp,
                              mutcols = mutcols,
                              basecols = basecols,
                              ref_idx = ref_idx,
                              lower = rep(-9, times = K - 1),
                              upper = rep(9, times = K - 1),
                              method = "L-BFGS-B",
                              hessian = TRUE)

        },
        error = function(e){
          message(paste0(c(paste0("Fitting of the generalized mixture model failed for following data: "), utils::capture.output(dataset)), collapse = "\n"))
          stop("Error message from optim: ", e)
        }
      )

      # Figure out what to call each of the fractions
      # Probably "p""muttype""old/new"_"muttype""old/new"
      muttypes <- colnames(mutrate_design)


      # Get vector of population statuses (new and old) for naming the output list
      population_list <- vector(mode = "list", length = K)
      for(i in 1:K){

        population_vector <- rep("", times = ncol(mutrate_design))
        count <- 1
        for(j in 1:ncol(mutrate_design)){

          if(as.logical(mutrate_design[i, j])){

            population_vector[count] <- "high"

          }else{

            population_vector[count] <- "low"

          }

          count <- count + 1


        }

        population_list[[i]] <- population_vector

      }


      # Make names for the output list that are easily interpretable (I hope)
      listnames <- lapply(population_list, function(x) paste0(x, muttypes))
      listnames <- sapply(listnames, function(x) paste(x, collapse = "_"),
                          simplify = "vector")

      listnames <- paste0(rep(c("invsoftmax_",
                                "fraction_",
                                "se_invsoftmax_",
                                "ci_lower_fraction_",
                                "ci_upper_fraction_"),
                              each = K),
                          listnames)

      full_alr <- numeric(K)
      full_alr[-ref_idx] <- fit$par
      full_alr[ref_idx] <- 0

      # Transform to centered log-ratio (CLR) coordinates so sum(clr_par) == 0,
      # softmax(clr_par) == softmax(full_alr), and every population gets a
      # reference-invariant SE from the nonsingular (K-1)x(K-1) ALR Hessian.
      clr_par <- full_alr - mean(full_alr)
      frac_est <- softmax(clr_par)

      uncertainties <- tryCatch(
        {
          cov_alr <- solve(fit$hessian)
          cov_full <- matrix(0, nrow = K, ncol = K)
          cov_full[-ref_idx, -ref_idx] <- cov_alr
          C_mat <- diag(K) - matrix(1 / K, nrow = K, ncol = K)
          cov_clr <- C_mat %*% cov_full %*% t(C_mat)
          sqrt(pmax(0, diag(cov_clr)))
        },
        error = function(e) {
          rep(Inf, times = K)
        }
      )

      # Profile-likelihood 95% CIs along relative-complement simplex rays
      n_obs <- nrow(dataset)
      E_1g <- matrix(0, nrow = n_obs, ncol = K)
      for(j in seq_len(ncol(mutrate_design))){
        mj <- as.numeric(dataset[[mutcols[j]]])
        bj <- as.numeric(dataset[[basecols[j]]])
        if(Poisson){
          ll_hi <- stats::dpois(mj, lambda = pnew[j] * bj, log = TRUE)
          ll_lo <- stats::dpois(mj, lambda = pold[j] * bj, log = TRUE)
        }else{
          ll_hi <- stats::dbinom(mj, bj, prob = pnew[j], log = TRUE)
          ll_lo <- stats::dbinom(mj, bj, prob = pold[j], log = TRUE)
        }
        for(k_i in seq_len(K)){
          E_1g[, k_i] <- E_1g[, k_i] + if(as.logical(mutrate_design[k_i, j])) ll_hi else ll_lo
        }
      }
      n_1g <- as.numeric(dataset[["n"]])
      N_1g <- sum(n_1g)
      ci_lo_vec <- numeric(K)
      ci_hi_vec <- numeric(K)
      for(k_i in seq_len(K)){
        other_idx <- setdiff(seq_len(K), k_i)
        rem_m <- sum(frac_est[other_idx])
        q_v <- if(rem_m > 1e-12) frac_est[other_idx] / rem_m else rep(1 / length(other_idx), length(other_idx))
        W_oth <- sweep(E_1g[, other_idx, drop = FALSE], 2L, log(pmax(q_v, 1e-300)), "+")
        w_mx <- apply(W_oth, 1L, max)
        w_mx_s <- ifelse(is.finite(w_mx), w_mx, 0)
        E_mk <- w_mx_s + log(rowSums(exp(W_oth - w_mx_s)))
        E_mk[!is.finite(w_mx)] <- -Inf
        pr_k <- profile_2comp_delta_vec(
          grp_id = rep(1L, n_obs),
          grp_end = n_obs,
          n = n_1g,
          delta = E_1g[, k_i] - E_mk,
          N_grp = N_1g
        )
        ci_lo_vec[k_i] <- pr_k$ci_lower_fn[1L]
        ci_hi_vec[k_i] <- pr_k$ci_upper_fn[1L]
      }

      outlist <- c(as.list(clr_par),
                   as.list(frac_est),
                   as.list(uncertainties),
                   as.list(ci_lo_vec),
                   as.list(ci_hi_vec))


      names(outlist) <- listnames

    }

  }


  return(outlist)


}



check_EstimateFractions_input <- function(args){

  ### features

  if(!is.character(args$features)){

    stop("features is not a character vector!")

  }

  ### Poisson

  if(!is.logical(args$Poisson)){

    stop("Poisson must be logical (TRUE or FALSE)!")

  }

  ### EBayes

  if(!is.null(args$EBayes) &&
     (!is.logical(args$EBayes) || length(args$EBayes) != 1L || is.na(args$EBayes))){

    stop("EBayes must be a single logical (TRUE or FALSE)!")

  }

  ### remove_multi_features & split_multi_features

  if(!is.null(args$remove_multi_features) &&
     (!is.logical(args$remove_multi_features) || length(args$remove_multi_features) != 1L || is.na(args$remove_multi_features))){

    stop("remove_multi_features must be a single logical (TRUE or FALSE)!")

  }

  if(!is.null(args$split_multi_features) &&
     (!is.logical(args$split_multi_features) || length(args$split_multi_features) != 1L || is.na(args$split_multi_features))){

    stop("split_multi_features must be a single logical (TRUE or FALSE)!")

  }

  if(isTRUE(args$remove_multi_features) && isTRUE(args$split_multi_features)){

    stop("Cannot set both remove_multi_features = TRUE and split_multi_features = TRUE!")

  }

  if(isTRUE(args$split_multi_features) &&
     (is.null(args$multi_feature_cols) || !is.character(args$multi_feature_cols) || length(args$multi_feature_cols) == 0L)){

    stop("multi_feature_cols must be provided as a non-empty character vector when split_multi_features = TRUE!")

  }

  if(!is.null(args$multi_feature_cols) &&
     (!is.character(args$multi_feature_cols) || length(args$multi_feature_cols) == 0L || any(is.na(args$multi_feature_cols)))){

    stop("multi_feature_cols must be NULL or a non-empty character vector!")

  }

  if(!is.null(args$multi_feature_sep) &&
     (!is.character(args$multi_feature_sep) || length(args$multi_feature_sep) != 1L ||
      is.na(args$multi_feature_sep) || !nzchar(args$multi_feature_sep))){

    stop("multi_feature_sep must be a single non-empty string!")

  }


}


#########################################################
# VECTORIZED NEWTON-RAPHSON MIXTURE SOLVERS (ZERO DEPS)
#########################################################

# Contiguous group sum with periodic cumsum resets to avoid floating-point drift
grp_sum_contig <- function(x, grp_end, block_groups = 2048L){
  n_g <- length(grp_end)
  if(n_g == 0L) return(numeric(0))
  n_rows <- length(x)
  if(n_g <= block_groups || n_rows <= 16384L){
    cs <- cumsum(x)
    if(n_g == 1L) return(cs[grp_end])
    return(c(cs[grp_end[1L]], cs[grp_end[-1L]] - cs[grp_end[-n_g]]))
  }
  out <- numeric(n_g)
  n_blocks <- as.integer(ceiling(n_g / block_groups))
  for(b in seq_len(n_blocks)){
    g_start <- (b - 1L) * block_groups + 1L
    g_end <- min(b * block_groups, n_g)
    r_start <- if(g_start == 1L) 1L else (grp_end[g_start - 1L] + 1L)
    r_end <- grp_end[g_end]
    cs <- cumsum(x[r_start:r_end])
    loc_end <- grp_end[g_start:g_end] - (r_start - 1L)
    out[g_start:g_end] <- if(length(loc_end) == 1L){
      cs[loc_end]
    }else{
      c(cs[loc_end[1L]], cs[loc_end[-1L]] - cs[loc_end[-length(loc_end)]])
    }
  }
  out
}


# Numerically stable softplus log(1 + exp(z))
softplus_vec <- function(z){
  pmax(z, 0) + log1p(exp(-abs(z)))
}


# Exact contiguous groupwise log-sum-exp (zero overflow/underflow)
grp_logsumexp <- function(log_w, grp_id, grp_end){
  ng <- length(grp_end)
  if(ng == 0L) return(numeric(0))
  grp_start <- c(1L, grp_end[-ng] + 1L)
  m_g <- numeric(ng)
  sizes <- grp_end - grp_start + 1L
  max_sz <- max(sizes)
  if(max_sz <= 64L){
    m_g <- log_w[grp_start]
    if(max_sz > 1L){
      for(offset in 1L:(max_sz - 1L)){
        active <- which(sizes > offset)
        if(length(active) == 0L) break
        m_g[active] <- pmax(m_g[active], log_w[grp_start[active] + offset])
      }
    }
  }else{
    small <- which(sizes <= 64L)
    if(length(small) > 0L){
      m_s <- log_w[grp_start[small]]
      ms_sz <- max(sizes[small])
      if(ms_sz > 1L){
        for(offset in 1L:(ms_sz - 1L)){
          act <- small[sizes[small] > offset]
          if(length(act) == 0L) break
          idx_s <- match(act, small)
          m_s[idx_s] <- pmax(m_s[idx_s], log_w[grp_start[act] + offset])
        }
      }
      m_g[small] <- m_s
    }
    big <- which(sizes > 64L)
    for(g in big){
      m_g[g] <- max(log_w[grp_start[g]:grp_end[g]])
    }
  }
  m_safe <- ifelse(is.finite(m_g), m_g, 0)
  s_exp <- grp_sum_contig(exp(log_w - m_safe[grp_id]), grp_end)
  res <- m_safe + log(s_exp)
  res[m_g == -Inf] <- -Inf
  res
}


# Vectorized pure-data MLE and exact 95% Wilks profile-likelihood CI for 2-comp mixture
profile_2comp_delta_vec <- function(grp_id, grp_end, n, delta, N_grp = NULL,
                                    crit = 0.5 * stats::qchisq(0.95, df = 1),
                                    eta_min = -25, eta_max = 25){
  ngroups <- length(grp_end)
  if(ngroups == 0L){
    return(list(
      eta_mle = numeric(0), info_mle = numeric(0),
      at_0 = logical(0), at_1 = logical(0),
      ci_lower_fn = numeric(0), ci_upper_fn = numeric(0),
      ci_lower_logit = numeric(0), ci_upper_logit = numeric(0)
    ))
  }
  if(is.null(N_grp)) N_grp <- grp_sum_contig(n, grp_end)
  delta_clean <- delta
  delta_clean[is.nan(delta_clean)] <- 0
  delta_clean <- pmin(pmax(delta_clean, -300), 300)

  log_n <- log(pmax(n, 1e-300))
  log_N <- log(pmax(N_grp, 1e-300))

  # Exact overflow-free simplex boundary score checks at theta = 0 and theta = 1
  lse0 <- grp_logsumexp(log_n + delta_clean, grp_id, grp_end)
  lse1 <- grp_logsumexp(log_n - delta_clean, grp_id, grp_end)

  at_0 <- lse0 <= log_N + 1e-13
  at_1 <- !at_0 & (lse1 <= log_N + 1e-13)

  eta_mle <- numeric(ngroups)
  eta_mle[at_0] <- -Inf
  eta_mle[at_1] <- Inf
  info_mle <- numeric(ngroups)

  # Relative log-likelihood l_rel(eta) = l(eta) - l(theta=0)
  ll_at_1 <- grp_sum_contig(n * delta_clean, grp_end)

  eval_rel <- function(eta_v, sg, se, sn, sdel, sN, want_ll = TRUE){
    theta <- stats::plogis(eta_v)
    z <- eta_v[sg] + sdel
    w <- stats::plogis(z)
    nw <- sn * w
    sum_nw <- grp_sum_contig(nw, se)
    g1 <- sum_nw - sN * theta
    info <- sN * theta * (1 - theta) - grp_sum_contig(nw * (1 - w), se)
    if(!want_ll) return(list(g1 = g1, info = info))
    ll_r <- grp_sum_contig(sn * (softplus_vec(z) - softplus_vec(eta_v[sg])), se)
    list(ll = ll_r, g1 = g1, info = info)
  }

  interior <- which(!(at_0 | at_1))
  if(length(interior) > 0L){
    active <- interior
    act_map <- integer(ngroups)
    act_map[active] <- seq_along(active)
    rk0 <- act_map[grp_id] > 0L
    sg <- act_map[grp_id[rk0]]
    gsz <- diff(c(0L, grp_end))[active]
    se <- cumsum(gsz)
    sn <- n[rk0]
    sdel <- delta_clean[rk0]
    sN <- N_grp[active]

    lo_a <- rep(eta_min, length(active))
    hi_a <- rep(eta_max, length(active))
    eta_a <- rep(0, length(active))

    g_lo <- eval_rel(lo_a, sg, se, sn, sdel, sN, want_ll = FALSE)$g1
    g_hi <- eval_rel(hi_a, sg, se, sn, sdel, sN, want_ll = FALSE)$g1
    pin_lo <- g_lo <= 0
    pin_hi <- !pin_lo & (g_hi >= 0)
    if(any(pin_lo | pin_hi)){
      if(any(pin_lo)) eta_mle[active[pin_lo]] <- eta_min
      if(any(pin_hi)) eta_mle[active[pin_hi]] <- eta_max
      keep <- which(!(pin_lo | pin_hi))
      if(length(keep) > 0L){
        remap <- integer(length(active))
        remap[keep] <- seq_along(keep)
        rk <- remap[sg] > 0L
        sg <- remap[sg[rk]]
        gsz <- gsz[keep]
        se <- cumsum(gsz)
        sn <- sn[rk]; sdel <- sdel[rk]; sN <- sN[keep]
        active <- active[keep]
        lo_a <- lo_a[keep]; hi_a <- hi_a[keep]; eta_a <- eta_a[keep]
      }else{
        active <- integer(0)
      }
    }

    if(length(active) > 0L){
      for(iter in seq_len(30L)){
        ev <- eval_rel(eta_a, sg, se, sn, sdel, sN, want_ll = FALSE)
        pos <- ev$g1 > 0
        lo_a[pos] <- pmax(lo_a[pos], eta_a[pos])
        hi_a[!pos] <- pmin(hi_a[!pos], eta_a[!pos])

        step <- ifelse(ev$info > 1e-12, ev$g1 / ev$info, NA_real_)
        cand <- eta_a + step
        bad <- is.na(cand) | cand <= lo_a | cand >= hi_a
        cand[bad] <- 0.5 * (lo_a[bad] + hi_a[bad])

        d_eta <- abs(cand - eta_a)
        eta_a <- cand
        conv <- (d_eta < 1e-9 & abs(ev$g1) < 1e-7) | ((hi_a - lo_a) < 1e-9)
        if(all(conv)){
          eta_mle[active] <- eta_a
          break
        }
        if(iter %% 3L == 0L && sum(conv) > 0.25 * length(conv)){
          done <- which(conv)
          eta_mle[active[done]] <- eta_a[done]
          rem <- which(!conv)
          active <- active[rem]
          eta_a <- eta_a[rem]; lo_a <- lo_a[rem]; hi_a <- hi_a[rem]
          sN <- sN[rem]; gsz <- gsz[rem]
          remap <- integer(length(conv))
          remap[rem] <- seq_along(rem)
          rk <- remap[sg] > 0L
          sg <- remap[sg[rk]]; se <- cumsum(gsz)
          sn <- sn[rk]; sdel <- sdel[rk]
        }
        if(iter == 30L) eta_mle[active] <- eta_a
      }
    }
  }

  ll_max <- numeric(ngroups)
  ll_max[at_0] <- 0
  ll_max[at_1] <- ll_at_1[at_1]
  if(length(interior) > 0L){
    act_map <- integer(ngroups)
    act_map[interior] <- seq_along(interior)
    rk0 <- act_map[grp_id] > 0L
    int_sz <- diff(c(0L, grp_end))[interior]
    ev_m <- eval_rel(eta_mle[interior], act_map[grp_id[rk0]], cumsum(int_sz),
                     n[rk0], delta_clean[rk0], N_grp[interior], want_ll = TRUE)
    ll_max[interior] <- pmax(ev_m$ll, 0, ll_at_1[interior])
    info_mle[interior] <- pmax(ev_m$info, 0)
  }

  ll_target <- ll_max - crit

  solve_prof_side <- function(idx_set, lo_init, hi_init, start_eta, increasing = TRUE){
    if(length(idx_set) == 0L) return(numeric(0))
    out_eta <- numeric(length(idx_set))
    active_loc <- seq_along(idx_set)
    act_map <- integer(ngroups)
    act_map[idx_set] <- active_loc
    rk0 <- act_map[grp_id] > 0L
    sg <- act_map[grp_id[rk0]]
    gsz <- diff(c(0L, grp_end))[idx_set]
    se <- cumsum(gsz)
    sn <- n[rk0]; sdel <- delta_clean[rk0]; sN <- N_grp[idx_set]
    targ <- ll_target[idx_set]

    lo_a <- lo_init
    hi_a <- hi_init
    eta_a <- pmin(pmax(start_eta, lo_a + 1e-6 * (hi_a - lo_a)), hi_a - 1e-6 * (hi_a - lo_a))

    for(iter in seq_len(45L)){
      ev <- eval_rel(eta_a, sg, se, sn, sdel, sN, want_ll = TRUE)
      f_val <- ev$ll - targ
      if(increasing){
        neg <- f_val < 0
        lo_a[neg] <- pmax(lo_a[neg], eta_a[neg])
        hi_a[!neg] <- pmin(hi_a[!neg], eta_a[!neg])
      }else{
        pos <- f_val > 0
        lo_a[pos] <- pmax(lo_a[pos], eta_a[pos])
        hi_a[!pos] <- pmin(hi_a[!pos], eta_a[!pos])
      }

      step <- ifelse(abs(ev$g1) > 1e-12, f_val / ev$g1, NA_real_)
      cand <- eta_a - step
      bad <- is.na(cand) | cand <= lo_a | cand >= hi_a
      cand[bad] <- 0.5 * (lo_a[bad] + hi_a[bad])

      d_eta <- abs(cand - eta_a)
      eta_a <- cand
      conv <- (abs(f_val) < 1e-9 & d_eta < 1e-8) | ((hi_a - lo_a) < 1e-9)
      if(all(conv)){
        out_eta[active_loc] <- eta_a
        break
      }
      if(iter %% 3L == 0L && sum(conv) > 0.25 * length(conv)){
        done <- which(conv)
        out_eta[active_loc[done]] <- eta_a[done]
        rem <- which(!conv)
        active_loc <- active_loc[rem]
        eta_a <- eta_a[rem]; lo_a <- lo_a[rem]; hi_a <- hi_a[rem]
        sN <- sN[rem]; targ <- targ[rem]; gsz <- gsz[rem]
        remap <- integer(length(conv))
        remap[rem] <- seq_along(rem)
        rk <- remap[sg] > 0L
        sg <- remap[sg[rk]]; se <- cumsum(gsz)
        sn <- sn[rk]; sdel <- sdel[rk]
      }
      if(iter == 45L) out_eta[active_loc] <- eta_a
    }
    out_eta
  }

  ci_lo_logit <- rep(-Inf, ngroups)
  ci_lo_fn <- numeric(ngroups)
  need_lo <- which(!at_0 & (ll_target > 0))
  if(length(need_lo) > 0L){
    hi_b <- ifelse(at_1[need_lo], 28, eta_mle[need_lo])
    lo_b <- rep(-28, length(need_lo))
    se_g <- ifelse(info_mle[need_lo] > 1e-6, 1 / sqrt(info_mle[need_lo]), 2.0)
    init_lo <- ifelse(at_1[need_lo], 3.0, eta_mle[need_lo] - 1.96 * se_g)
    eta_lo_res <- solve_prof_side(need_lo, lo_b, hi_b, init_lo, increasing = TRUE)
    ci_lo_logit[need_lo] <- eta_lo_res
    ci_lo_fn[need_lo] <- stats::plogis(eta_lo_res)
  }

  ci_hi_logit <- rep(Inf, ngroups)
  ci_hi_fn <- rep(1, ngroups)
  need_hi <- which(!at_1 & (ll_at_1 < ll_target))
  if(length(need_hi) > 0L){
    lo_b <- ifelse(at_0[need_hi], -28, eta_mle[need_hi])
    hi_b <- rep(28, length(need_hi))
    se_g <- ifelse(info_mle[need_hi] > 1e-6, 1 / sqrt(info_mle[need_hi]), 2.0)
    init_hi <- ifelse(at_0[need_hi], -3.0, eta_mle[need_hi] + 1.96 * se_g)
    eta_hi_res <- solve_prof_side(need_hi, lo_b, hi_b, init_hi, increasing = FALSE)
    ci_hi_logit[need_hi] <- eta_hi_res
    ci_hi_fn[need_hi] <- stats::plogis(eta_hi_res)
  }

  list(
    eta_mle = eta_mle,
    info_mle = info_mle,
    at_0 = at_0,
    at_1 = at_1,
    ci_lower_fn = ci_lo_fn,
    ci_upper_fn = ci_hi_fn,
    ci_lower_logit = ci_lo_logit,
    ci_upper_logit = ci_hi_logit
  )
}


# Per-sample Empirical Bayes calibration of logit(fraction) prior location and scale
calibrate_ebayes_logit_prior <- function(sample_grp, eta_mle, info_mle, at_0, at_1,
                                         default_mean = 0, default_sd = 1.7,
                                         sd_min = 0.8, sd_max = 2.5){
  ng <- length(eta_mle)
  if(ng == 0L){
    return(list(prior_mean = numeric(0), prior_sd = numeric(0)))
  }
  prior_mean <- rep(default_mean, ng)
  prior_sd <- rep(default_sd, ng)
  u_samps <- unique(sample_grp)

  for(s in u_samps){
    idx_s <- which(sample_grp == s)
    ok_base <- !at_0[idx_s] & !at_1[idx_s] & is.finite(eta_mle[idx_s])
    ok1 <- idx_s[ok_base & abs(eta_mle[idx_s]) < 7 & info_mle[idx_s] > 1]
    ok <- if(length(ok1) >= 10L){
      ok1
    }else{
      ok2 <- idx_s[ok_base & abs(eta_mle[idx_s]) < 8 & info_mle[idx_s] > 0.1]
      if(length(ok2) >= 5L) ok2 else integer(0)
    }
    if(length(ok) > 0L){
      eta_ok <- eta_mle[ok]
      info_ok <- info_mle[ok]
      mu_s <- pmin(pmax(stats::median(eta_ok), -5), 5)
      raw_var <- stats::var(eta_ok) - mean(1 / info_ok)
      sd_s <- if(is.na(raw_var) || raw_var <= 0){
        sd_min
      }else{
        pmin(pmax(sqrt(raw_var), sd_min), sd_max)
      }
      prior_mean[idx_s] <- mu_s
      prior_sd[idx_s] <- sd_s
    }
  }
  list(prior_mean = prior_mean, prior_sd = prior_sd)
}


# Globally convergent safeguarded Newton-bisection for 1D two-component mixture
fit_tcmm_known_vec <- function(grp_id, grp_end, muts, nucs, n, pnew, pold,
                               Poisson = TRUE,
                               prior_mean = 0, prior_sd = 1.7,
                               EBayes = FALSE, sample_grp = NULL,
                               lower = if(EBayes) -16 else -9,
                               upper = if(EBayes) 16 else 9,
                               tol = 1e-10, max_iter = 50L){
  ngroups <- length(grp_end)
  if(ngroups == 0L){
    return(list(
      logit_fn = numeric(0), se = numeric(0), n = numeric(0),
      ci_lower_fn = numeric(0), ci_upper_fn = numeric(0)
    ))
  }

  muts <- as.numeric(muts)
  nucs <- as.numeric(nucs)
  n <- as.numeric(n)

  if(Poisson){
    ll_new <- stats::dpois(muts, nucs * pnew, log = TRUE)
    ll_old <- stats::dpois(muts, nucs * pold, log = TRUE)
  }else{
    ll_new <- stats::dbinom(muts, nucs, pnew, log = TRUE)
    ll_old <- stats::dbinom(muts, nucs, pold, log = TRUE)
  }
  delta <- ll_new - ll_old
  delta[is.nan(delta)] <- 0
  delta <- pmin(pmax(delta, -300), 300)

  N_grp <- grp_sum_contig(n, grp_end)

  # Exact pure-data MLEs and 95% Wilks profile-likelihood CIs
  prof <- profile_2comp_delta_vec(grp_id, grp_end, n, delta, N_grp = N_grp)

  if(EBayes){
    if(is.null(sample_grp)) sample_grp <- rep(1L, ngroups)
    eb <- calibrate_ebayes_logit_prior(
      sample_grp = sample_grp,
      eta_mle = prof$eta_mle,
      info_mle = prof$info_mle,
      at_0 = prof$at_0,
      at_1 = prof$at_1,
      default_mean = prior_mean[1L],
      default_sd = prior_sd[1L]
    )
    pm_all <- eb$prior_mean
    is2_all <- 1 / (eb$prior_sd * eb$prior_sd)
  }else{
    pm_all <- rep_len(as.numeric(prior_mean), ngroups)
    ps_all <- rep_len(as.numeric(prior_sd), ngroups)
    is2_all <- 1 / (ps_all * ps_all)
  }

  M_ph <- 1.25

  eval_gh <- function(eta_vec, sub_grp, sub_end, sub_n, sub_delta, sub_N, sub_pm, sub_is2){
    theta <- stats::plogis(eta_vec)
    w <- stats::plogis(eta_vec[sub_grp] + sub_delta)
    nw <- sub_n * w
    nww <- nw * (1 - w)
    sum_nw <- grp_sum_contig(nw, sub_end)
    sum_nww <- grp_sum_contig(nww, sub_end)
    diff_e <- eta_vec - sub_pm
    if(EBayes){
      z_ph <- (diff_e * sub_is2) / M_ph
      rz_ph <- sqrt(1 + z_ph * z_ph)
      g_pr <- (diff_e * sub_is2) / rz_ph
      h_pr <- sub_is2 / (rz_ph * (1 + z_ph * z_ph))
    }else{
      g_pr <- diff_e * sub_is2
      h_pr <- sub_is2
    }
    g <- sub_N * theta - sum_nw + g_pr
    h <- sub_N * theta * (1 - theta) - sum_nww + h_pr
    list(g = g, h = h)
  }

  lo <- rep(lower, ngroups)
  hi <- rep(upper, ngroups)
  eta <- ifelse(prof$at_0, lower + 1,
                ifelse(prof$at_1, upper - 1,
                       pmin(pmax(prof$eta_mle, lower + 1e-4), upper - 1e-4)))

  gh_lo <- eval_gh(lo, grp_id, grp_end, n, delta, N_grp, pm_all, is2_all)
  gh_hi <- eval_gh(hi, grp_id, grp_end, n, delta, N_grp, pm_all, is2_all)

  at_lo <- gh_lo$g >= 0
  at_hi <- !at_lo & (gh_hi$g <= 0)
  eta[at_lo] <- lower
  eta[at_hi] <- upper
  h_final <- numeric(ngroups)
  h_final[at_lo] <- gh_lo$h[at_lo]
  h_final[at_hi] <- gh_hi$h[at_hi]

  interior <- which(!(at_lo | at_hi))
  if(length(interior) > 0L){
    active <- interior
    act_map <- integer(ngroups)
    act_map[active] <- seq_along(active)
    row_keep <- act_map[grp_id] > 0L
    sub_grp <- act_map[grp_id[row_keep]]
    grp_sizes <- diff(c(0L, grp_end))[active]
    sub_end <- cumsum(grp_sizes)
    sub_n <- n[row_keep]
    sub_delta <- delta[row_keep]
    sub_N <- N_grp[active]
    sub_pm <- pm_all[active]
    sub_is2 <- is2_all[active]
    eta_a <- eta[active]
    lo_a <- lo[active]
    hi_a <- hi[active]

    for(iter in seq_len(max_iter)){
      gh <- eval_gh(eta_a, sub_grp, sub_end, sub_n, sub_delta, sub_N, sub_pm, sub_is2)
      g <- gh$g
      h <- gh$h

      neg <- g < 0
      lo_a[neg] <- pmax(lo_a[neg], eta_a[neg])
      hi_a[!neg] <- pmin(hi_a[!neg], eta_a[!neg])

      step <- ifelse(h > 1e-12, g / h, NA_real_)
      cand <- eta_a - step
      bad <- is.na(cand) | cand <= lo_a | cand >= hi_a
      cand[bad] <- 0.5 * (lo_a[bad] + hi_a[bad])

      diff_eta <- abs(cand - eta_a)
      eta_a <- cand

      conv <- (diff_eta < tol & abs(g) < 1e-7) | ((hi_a - lo_a) < tol)
      if(all(conv)){
        eta[active] <- eta_a
        break
      }
      if(iter %% 3L == 0L && sum(conv) > 0.25 * length(conv)){
        done_local <- which(conv)
        eta[active[done_local]] <- eta_a[done_local]

        rem_local <- which(!conv)
        active <- active[rem_local]
        eta_a <- eta_a[rem_local]
        lo_a <- lo_a[rem_local]
        hi_a <- hi_a[rem_local]
        sub_N <- sub_N[rem_local]
        sub_pm <- sub_pm[rem_local]
        sub_is2 <- sub_is2[rem_local]
        grp_sizes <- grp_sizes[rem_local]

        remap <- integer(length(conv))
        remap[rem_local] <- seq_along(rem_local)
        rk <- remap[sub_grp] > 0L
        sub_grp <- remap[sub_grp[rk]]
        sub_end <- cumsum(grp_sizes)
        sub_n <- sub_n[rk]
        sub_delta <- sub_delta[rk]
      }
      if(iter == max_iter){
        eta[active] <- eta_a
      }
    }

    if(!EBayes){
      act_map <- integer(ngroups)
      act_map[interior] <- seq_along(interior)
      rk <- act_map[grp_id] > 0L
      int_sizes <- diff(c(0L, grp_end))[interior]
      gh_fin <- eval_gh(eta[interior], act_map[grp_id[rk]], cumsum(int_sizes),
                        n[rk], delta[rk], N_grp[interior],
                        pm_all[interior], is2_all[interior])
      h_final[interior] <- gh_fin$h
    }
  }

  if(EBayes){
    se_out <- posterior_se_from_ci(
      est = eta,
      ci_lo = prof$ci_lower_logit,
      ci_hi = prof$ci_upper_logit,
      prior_sd = eb$prior_sd,
      se_min = 1 / sqrt(pmax(N_grp, 1))
    )
  }else{
    theta_fin <- stats::plogis(eta)
    w_fin <- stats::plogis(eta[grp_id] + delta)
    h_fisher_final <- grp_sum_contig(n * (w_fin - theta_fin[grp_id])^2, grp_end)
    se_out <- sqrt(1 / pmax(h_final, h_fisher_final + is2_all, 1e-12))
  }

  list(
    logit_fn = eta,
    se = se_out,
    n = N_grp,
    ci_lower_fn = prof$ci_lower_fn,
    ci_upper_fn = prof$ci_upper_fn
  )
}


# Table-level wrapper for standard 2-component fraction estimation
solve_tcmm_standard_dt <- function(dt, features_to_analyze, pop_col, base_col,
                                   samples_with_no_label, Poisson,
                                   col_name, uncertainty_col, natural_col_name,
                                   EBayes = TRUE){
  sample <- n <- NULL
  `.` <- list

  ci_lower_col <- paste0("ci_lower_", natural_col_name)
  ci_upper_col <- paste0("ci_upper_", natural_col_name)

  grp_cols <- c("sample", features_to_analyze)
  dt_work <- data.table::as.data.table(dt)

  if(nrow(dt_work) == 0L){
    empty_dt <- dt_work[, ..grp_cols]
    empty_dt[[natural_col_name]] <- numeric(0)
    empty_dt[[col_name]] <- numeric(0)
    empty_dt[[uncertainty_col]] <- numeric(0)
    empty_dt[[ci_lower_col]] <- numeric(0)
    empty_dt[[ci_upper_col]] <- numeric(0)
    empty_dt[["n"]] <- numeric(0)
    return(dplyr::as_tibble(empty_dt))
  }

  is_nolabel <- dt_work$sample %in% samples_with_no_label
  out_parts <- list()

  if(any(is_nolabel)){
    nl_dt <- dt_work[is_nolabel, .(n = as.numeric(sum(n))), by = grp_cols]
    nl_dt[[natural_col_name]] <- 0
    nl_dt[[col_name]] <- -Inf
    nl_dt[[uncertainty_col]] <- 0
    nl_dt[[ci_lower_col]] <- 0
    nl_dt[[ci_upper_col]] <- 0
    out_cols <- c("sample", features_to_analyze, natural_col_name, col_name,
                  uncertainty_col, ci_lower_col, ci_upper_col, "n")
    out_parts[[length(out_parts) + 1L]] <- nl_dt[, ..out_cols]
  }

  if(any(!is_nolabel)){
    lb_dt <- dt_work[!is_nolabel]
    data.table::setorderv(lb_dt, grp_cols)
    grp_id <- data.table::rleidv(lb_dt, cols = grp_cols)
    n_rows <- nrow(lb_dt)
    grp_end <- if(n_rows > 1L){
      c(which(grp_id[-1L] != grp_id[-n_rows]), n_rows)
    }else{
      1L
    }

    fit_res <- fit_tcmm_known_vec(
      grp_id = grp_id,
      grp_end = grp_end,
      muts = lb_dt[[pop_col]],
      nucs = lb_dt[[base_col]],
      n = lb_dt[["n"]],
      pnew = lb_dt[["pnew"]],
      pold = lb_dt[["pold"]],
      Poisson = Poisson,
      EBayes = EBayes,
      sample_grp = lb_dt[["sample"]][grp_end]
    )

    lb_out <- lb_dt[grp_end, ..grp_cols]
    lb_out[[natural_col_name]] <- inv_logit(fit_res$logit_fn)
    lb_out[[col_name]] <- fit_res$logit_fn
    lb_out[[uncertainty_col]] <- fit_res$se
    lb_out[[ci_lower_col]] <- fit_res$ci_lower_fn
    lb_out[[ci_upper_col]] <- fit_res$ci_upper_fn
    lb_out[["n"]] <- fit_res$n
    out_parts[[length(out_parts) + 1L]] <- lb_out
  }

  res_dt <- data.table::rbindlist(out_parts, use.names = TRUE)
  data.table::setorderv(res_dt, grp_cols)
  dplyr::as_tibble(res_dt)
}


# Vectorized 2D damped projected Newton solver for hierarchical pnew model
solve_tcmm_hier_dt <- function(dt, features_to_analyze, pop_col, base_col,
                               samples_with_no_label, Poisson,
                               pnew_prior_mean_col, pnew_prior_sd_arg,
                               fraction_prior_sd = 1.7,
                               EBayes = FALSE,
                               compute_ci = FALSE){
  sample <- n <- pold <- NULL
  `.` <- list

  grp_cols <- c("sample", features_to_analyze)
  dt_work <- data.table::as.data.table(dt)

  if(nrow(dt_work) == 0L){
    empty_dt <- dt_work[, ..grp_cols]
    empty_dt[, c("p1", "p2", "p1_u", "p2_u", "pold", "n") :=
               .(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0))]
    if(compute_ci){
      empty_dt[, c("ci_lower_fn", "ci_upper_fn") := .(numeric(0), numeric(0))]
    }
    return(empty_dt)
  }

  is_nolabel <- dt_work$sample %in% samples_with_no_label
  out_parts <- list()

  if(any(is_nolabel)){
    nl_dt <- dt_work[is_nolabel, .(
      p1 = -Inf,
      p2 = -Inf,
      p1_u = 0,
      p2_u = 0,
      pold = as.numeric(pold[1L]),
      n = as.numeric(sum(n))
    ), by = grp_cols]
    if(compute_ci){
      nl_dt[, c("ci_lower_fn", "ci_upper_fn") := .(0, 0)]
    }
    out_parts[[length(out_parts) + 1L]] <- nl_dt
  }

  if(any(!is_nolabel)){
    lb_dt <- dt_work[!is_nolabel]
    data.table::setorderv(lb_dt, grp_cols)
    grp_id <- data.table::rleidv(lb_dt, cols = grp_cols)
    n_rows <- nrow(lb_dt)
    grp_end <- if(n_rows > 1L){
      c(which(grp_id[-1L] != grp_id[-n_rows]), n_rows)
    }else{
      1L
    }
    ngroups <- length(grp_end)

    muts <- as.numeric(lb_dt[[pop_col]])
    nucs <- as.numeric(lb_dt[[base_col]])
    n_vec <- as.numeric(lb_dt[["n"]])
    pold_vec <- as.numeric(lb_dt[["pold"]])

    N_grp <- grp_sum_contig(n_vec, grp_end)
    m2_all <- as.numeric(lb_dt[[pnew_prior_mean_col]][grp_end])
    s2_all <- if(is.character(pnew_prior_sd_arg)){
      as.numeric(lb_dt[[pnew_prior_sd_arg]][grp_end])
    }else{
      rep(as.numeric(pnew_prior_sd_arg), ngroups)
    }

    inv_s2_all <- 1 / (s2_all * s2_all)

    if(Poisson){
      ll_old <- stats::dpois(muts, nucs * pold_vec, log = TRUE)
      c_const <- ifelse(muts == 0 & nucs == 0, 0,
                        muts * log(pmax(nucs, 1e-300)) - lgamma(muts + 1))
    }else{
      ll_old <- stats::dbinom(muts, nucs, pold_vec, log = TRUE)
      c_const <- lchoose(nucs, muts)
    }

    if(EBayes){
      pnew_init_r <- stats::plogis(m2_all[grp_id])
      ll_new_init <- if(Poisson){
        stats::dpois(muts, nucs * pnew_init_r, log = TRUE)
      }else{
        stats::dbinom(muts, nucs, pnew_init_r, log = TRUE)
      }
      prof_eb <- profile_2comp_delta_vec(grp_id, grp_end, n_vec, ll_new_init - ll_old, N_grp = N_grp)
      eb_cal <- calibrate_ebayes_logit_prior(
        sample_grp = lb_dt[["sample"]][grp_end],
        eta_mle = prof_eb$eta_mle,
        info_mle = prof_eb$info_mle,
        at_0 = prof_eb$at_0,
        at_1 = prof_eb$at_1,
        default_mean = 0,
        default_sd = fraction_prior_sd
      )
      m1_all <- eb_cal$prior_mean
      s1_all <- eb_cal$prior_sd
      inv_s1_all <- 1 / (s1_all * s1_all)
      e1_lo_bound <- -16
      e1_hi_bound <- 16
    }else{
      m1_all <- rep(0, ngroups)
      s1_all <- rep(fraction_prior_sd, ngroups)
      inv_s1_all <- rep(1 / (fraction_prior_sd * fraction_prior_sd), ngroups)
      e1_lo_bound <- -9
      e1_hi_bound <- 9
    }

    M_ph <- 1.25

    eval_fgh_2d <- function(e1, e2, sg, se, sn, smuts, snucs, sc, sllo, sN,
                            sm1, is1, sm2, is2, use_ph = EBayes){
      theta <- stats::plogis(e1)
      pnew_g <- stats::plogis(e2)
      e1_r <- e1[sg]
      e2_r <- e2[sg]
      pn_r <- pnew_g[sg]

      if(Poisson){
        log_pn <- ifelse(e2_r > 0, -log1p(exp(-e2_r)), e2_r - log1p(exp(e2_r)))
        ll_new <- sc + smuts * log_pn - snucs * pn_r
        d1 <- (1 - pn_r) * (smuts - snucs * pn_r)
        d2 <- -pn_r * (1 - pn_r) * (smuts + snucs * (1 - 2 * pn_r))
      }else{
        sp2 <- ifelse(e2_r > 0, e2_r + log1p(exp(-e2_r)), log1p(exp(e2_r)))
        ll_new <- sc + smuts * e2_r - snucs * sp2
        d1 <- smuts - snucs * pn_r
        d2 <- -snucs * pn_r * (1 - pn_r)
      }

      lt1 <- ifelse(e1_r > 0, -log1p(exp(-e1_r)), e1_r - log1p(exp(e1_r)))
      lt0 <- ifelse(e1_r > 0, -e1_r - log1p(exp(-e1_r)), -log1p(exp(e1_r)))
      ll_obs <- log_add2(lt1 + ll_new, lt0 + sllo)

      delta <- ll_new - sllo
      delta[is.nan(delta)] <- 0
      w <- stats::plogis(e1_r + delta)
      ww <- w * (1 - w)

      diff1 <- e1 - sm1
      if(use_ph){
        z1 <- (diff1 * is1) / M_ph
        rz1 <- sqrt(1 + z1 * z1)
        pen1 <- ifelse(is1 > 0, (M_ph * M_ph / is1) * (rz1 - 1), 0)
        g1_pr <- (diff1 * is1) / rz1
        h11_pr <- is1 / (rz1 * (1 + z1 * z1))
      }else{
        pen1 <- 0.5 * diff1 * diff1 * is1
        g1_pr <- diff1 * is1
        h11_pr <- is1
      }

      f_val <- -grp_sum_contig(sn * ll_obs, se) +
        pen1 + 0.5 * (e2 - sm2)^2 * is2

      nw <- sn * w
      nww <- sn * ww
      g1 <- sN * theta - grp_sum_contig(nw, se) + g1_pr
      g2 <- -grp_sum_contig(nw * d1, se) + (e2 - sm2) * is2

      h11 <- sN * theta * (1 - theta) - grp_sum_contig(nww, se) + h11_pr
      h12 <- -grp_sum_contig(nww * d1, se)
      h22 <- -grp_sum_contig(nww * d1 * d1 + nw * d2, se) + is2

      list(f = f_val, g1 = g1, g2 = g2, h11 = h11, h12 = h12, h22 = h22)
    }

    e1_all <- if(EBayes) pmin(pmax(m1_all, -6), 6) else rep(0, ngroups)
    e2_all <- pmin(pmax(m2_all, -8.5), -0.1)

    active <- seq_len(ngroups)
    sub_grp <- grp_id
    grp_sizes <- diff(c(0L, grp_end))
    sub_end <- grp_end
    sub_n <- n_vec; sub_muts <- muts; sub_nucs <- nucs
    sub_c <- c_const; sub_llo <- ll_old
    sub_N <- N_grp
    sub_m1 <- m1_all; sub_is1 <- inv_s1_all
    sub_m2 <- m2_all; sub_is2 <- inv_s2_all
    e1_a <- e1_all; e2_a <- e2_all

    for(iter in seq_len(35L)){
      res <- eval_fgh_2d(e1_a, e2_a, sub_grp, sub_end, sub_n, sub_muts, sub_nucs,
                         sub_c, sub_llo, sub_N, sub_m1, sub_is1, sub_m2, sub_is2)
      pin1 <- (e1_a <= e1_lo_bound & res$g1 > 0) | (e1_a >= e1_hi_bound & res$g1 < 0)
      pin2 <- (e2_a <= -9 & res$g2 > 0) | (e2_a >= 0 & res$g2 < 0)
      pg1 <- ifelse(pin1, 0, res$g1)
      pg2 <- ifelse(pin2, 0, res$g2)
      conv <- pmax(abs(pg1), abs(pg2)) < 1e-7

      if(all(conv)){
        e1_all[active] <- e1_a
        e2_all[active] <- e2_a
        break
      }

      if(iter %% 3L == 0L && sum(conv) > 0.25 * length(conv)){
        done_loc <- which(conv)
        e1_all[active[done_loc]] <- e1_a[done_loc]
        e2_all[active[done_loc]] <- e2_a[done_loc]

        rem_loc <- which(!conv)
        active <- active[rem_loc]
        e1_a <- e1_a[rem_loc]; e2_a <- e2_a[rem_loc]
        sub_N <- sub_N[rem_loc]
        sub_m1 <- sub_m1[rem_loc]; sub_is1 <- sub_is1[rem_loc]
        sub_m2 <- sub_m2[rem_loc]; sub_is2 <- sub_is2[rem_loc]
        grp_sizes <- grp_sizes[rem_loc]

        remap <- integer(length(conv))
        remap[rem_loc] <- seq_along(rem_loc)
        rk <- remap[sub_grp] > 0L
        sub_grp <- remap[sub_grp[rk]]
        sub_end <- cumsum(grp_sizes)
        sub_n <- sub_n[rk]; sub_muts <- sub_muts[rk]; sub_nucs <- sub_nucs[rk]
        sub_c <- sub_c[rk]; sub_llo <- sub_llo[rk]

        res <- lapply(res, `[`, rem_loc)
        pin1 <- pin1[rem_loc]; pin2 <- pin2[rem_loc]
        pg1 <- pg1[rem_loc]; pg2 <- pg2[rem_loc]
      }

      h11 <- res$h11; h12 <- res$h12; h22 <- res$h22
      h12[pin1 | pin2] <- 0
      h11[pin1] <- 1; h22[pin2] <- 1

      tr <- h11 + h22
      det_h <- h11 * h22 - h12 * h12
      min_eig <- 0.5 * (tr - sqrt(pmax(0, tr * tr - 4 * det_h)))
      shift <- ifelse(min_eig < 1e-6, (1e-4 - min_eig), 0)
      h11 <- h11 + shift; h22 <- h22 + shift
      det_h <- h11 * h22 - h12 * h12

      d1_s <- -(h22 * pg1 - h12 * pg2) / det_h
      d2_s <- -(-h12 * pg1 + h11 * pg2) / det_h
      max_s <- pmax(abs(d1_s), abs(d2_s))
      sc_s <- ifelse(max_s > 2.5, 2.5 / max_s, 1)
      d1_s <- d1_s * sc_s; d2_s <- d2_s * sc_s

      c1 <- pmin(pmax(e1_a + d1_s, e1_lo_bound), e1_hi_bound)
      c2 <- pmin(pmax(e2_a + d2_s, -9), 0)
      alpha <- rep(1, length(active))

      for(bt in seq_len(8L)){
        f_new <- eval_fgh_2d(c1, c2, sub_grp, sub_end, sub_n, sub_muts, sub_nucs,
                             sub_c, sub_llo, sub_N, sub_m1, sub_is1, sub_m2, sub_is2)$f
        bad <- is.na(f_new) | (f_new > res$f + 1e-10)
        if(!any(bad)) break
        alpha[bad] <- alpha[bad] * 0.5
        if(bt >= 5L){
          gn <- pmax(sqrt(pg1[bad]^2 + pg2[bad]^2), 1e-8)
          c1[bad] <- pmin(pmax(e1_a[bad] - alpha[bad] * pg1[bad] / gn, e1_lo_bound), e1_hi_bound)
          c2[bad] <- pmin(pmax(e2_a[bad] - alpha[bad] * pg2[bad] / gn, -9), 0)
        }else{
          c1[bad] <- pmin(pmax(e1_a[bad] + alpha[bad] * d1_s[bad], e1_lo_bound), e1_hi_bound)
          c2[bad] <- pmin(pmax(e2_a[bad] + alpha[bad] * d2_s[bad], -9), 0)
        }
      }
      e1_a <- c1; e2_a <- c2
      if(iter == 35L){
        e1_all[active] <- e1_a
        e2_all[active] <- e2_a
      }
    }

    fin <- eval_fgh_2d(e1_all, e2_all, grp_id, grp_end, n_vec, muts, nucs,
                       c_const, ll_old, N_grp, m1_all, inv_s1_all, m2_all, inv_s2_all)
    pin1_f <- (e1_all <= e1_lo_bound & fin$g1 > 0) | (e1_all >= e1_hi_bound & fin$g1 < 0)
    pin2_f <- (e2_all <= -9 & fin$g2 > 0) | (e2_all >= 0 & fin$g2 < 0)
    pg_max <- pmax(abs(ifelse(pin1_f, 0, fin$g1)), abs(ifelse(pin2_f, 0, fin$g2)))
    det_fin <- fin$h11 * fin$h22 - fin$h12 * fin$h12
    p1_u <- sqrt(abs(fin$h22 / pmax(abs(det_fin), 1e-30)))
    p2_u <- sqrt(abs(fin$h11 / pmax(abs(det_fin), 1e-30)))

    # Fallback to fit_tcmm for any unconverged or ill-conditioned group
    need_fb <- which(is.na(pg_max) | pg_max > 1e-4 | det_fin <= 1e-14 |
                       !is.finite(p1_u) | !is.finite(p2_u))
    if(length(need_fb) > 0L){
      grp_start <- c(1L, grp_end[-ngroups] + 1L)
      for(g in need_fb){
        idx <- grp_start[g]:grp_end[g]
        fb <- fit_tcmm(
          muts = muts[idx], nucs = nucs[idx], n = n_vec[idx],
          Poisson = Poisson, pold = pold_vec[idx[1L]],
          pnew_prior_mean = m2_all[g], pnew_prior_sd = s2_all[g],
          fraction_prior_mean = m1_all[g], fraction_prior_sd = s1_all[g]
        )
        e1_all[g] <- fb$p1; e2_all[g] <- fb$p2
        p1_u[g] <- fb$p1_u; p2_u[g] <- fb$p2_u
      }
    }

    if(compute_ci){
      # 2D profile-likelihood 95% CIs profiling over feature-specific logit(pnew)
      pn_hat_r <- stats::plogis(e2_all[grp_id])
      ll_new_hat <- if(Poisson){
        stats::dpois(muts, nucs * pn_hat_r, log = TRUE)
      }else{
        stats::dbinom(muts, nucs, pn_hat_r, log = TRUE)
      }
      prof_init <- profile_2comp_delta_vec(grp_id, grp_end, n_vec, ll_new_hat - ll_old, N_grp = N_grp)

      eval_prof_e1 <- function(e1_v, e2_start, sg, se, sn, smuts, snucs, sc, sllo, sN, sm2, is2){
        z0 <- numeric(length(e1_v))
        e2_cur <- e2_start
        for(st_i in seq_len(3L)){
          ev_p <- eval_fgh_2d(e1_v, e2_cur, sg, se, sn, smuts, snucs, sc, sllo, sN,
                              z0, z0, sm2, is2, use_ph = FALSE)
          step2 <- ev_p$g2 / pmax(ev_p$h22, 1e-4)
          step2 <- pmin(pmax(step2, -1.5), 1.5)
          e2_cur <- pmin(pmax(e2_cur - step2, -9), 0)
        }
        ev_p <- eval_fgh_2d(e1_v, e2_cur, sg, se, sn, smuts, snucs, sc, sllo, sN,
                            z0, z0, sm2, is2, use_ph = FALSE)
        list(ll = -ev_p$f, d_ll = -ev_p$g1, e2 = e2_cur)
      }

      crit_h <- 0.5 * stats::qchisq(0.95, df = 1)
      ll_prof_0 <- grp_sum_contig(n_vec * ll_old, grp_end)
      pr_1 <- eval_prof_e1(rep(28, ngroups), e2_all, grp_id, grp_end, n_vec, muts, nucs,
                           c_const, ll_old, N_grp, m2_all, inv_s2_all)
      ll_prof_1 <- pr_1$ll

      e1_mle_init <- ifelse(prof_init$at_0, -25,
                            ifelse(prof_init$at_1, 25,
                                   pmin(pmax(prof_init$eta_mle, -25), 25)))
      pr_m <- eval_prof_e1(e1_mle_init, e2_all, grp_id, grp_end, n_vec, muts, nucs,
                           c_const, ll_old, N_grp, m2_all, inv_s2_all)
      pr_hat <- eval_prof_e1(e1_all, e2_all, grp_id, grp_end, n_vec, muts, nucs,
                             c_const, ll_old, N_grp, m2_all, inv_s2_all)
      use_hat <- pr_hat$ll > pr_m$ll
      e1_best <- ifelse(use_hat, e1_all, e1_mle_init)
      ll_int_best <- pmax(pr_m$ll, pr_hat$ll)
      best_is_0 <- (ll_prof_0 >= ll_int_best - 1e-10) & (ll_prof_0 >= ll_prof_1)
      best_is_1 <- !best_is_0 & (ll_prof_1 >= ll_int_best - 1e-10)
      e1_best[best_is_0] <- -28
      e1_best[best_is_1] <- 28
      ll_max_h <- pmax(ll_int_best, ll_prof_0, ll_prof_1)
      ll_targ_h <- ll_max_h - crit_h

      solve_hier_side <- function(idx_set, lo_b, hi_b, init_e1, increasing = TRUE){
        if(length(idx_set) == 0L) return(numeric(0))
        act_map <- integer(ngroups)
        act_map[idx_set] <- seq_along(idx_set)
        rk <- act_map[grp_id] > 0L
        sg <- act_map[grp_id[rk]]
        se <- cumsum(diff(c(0L, grp_end))[idx_set])
        sn <- n_vec[rk]; smuts <- muts[rk]; snucs <- nucs[rk]
        sc <- c_const[rk]; sllo <- ll_old[rk]
        sN <- N_grp[idx_set]; sm2 <- m2_all[idx_set]; is2 <- inv_s2_all[idx_set]
        targ <- ll_targ_h[idx_set]
        e2_st <- e2_all[idx_set]
        lo_a <- lo_b; hi_a <- hi_b
        eta_a <- pmin(pmax(init_e1, lo_a + 1e-5 * (hi_a - lo_a)), hi_a - 1e-5 * (hi_a - lo_a))

        for(it_h in seq_len(14L)){
          ev_s <- eval_prof_e1(eta_a, e2_st, sg, se, sn, smuts, snucs, sc, sllo, sN, sm2, is2)
          e2_st <- ev_s$e2
          fv <- ev_s$ll - targ
          if(increasing){
            neg <- fv < 0
            lo_a[neg] <- pmax(lo_a[neg], eta_a[neg])
            hi_a[!neg] <- pmin(hi_a[!neg], eta_a[!neg])
          }else{
            pos <- fv > 0
            lo_a[pos] <- pmax(lo_a[pos], eta_a[pos])
            hi_a[!pos] <- pmin(hi_a[!pos], eta_a[!pos])
          }
          stp <- ifelse(abs(ev_s$d_ll) > 1e-10, fv / ev_s$d_ll, NA_real_)
          cand <- eta_a - stp
          bad <- is.na(cand) | cand <= lo_a | cand >= hi_a
          cand[bad] <- 0.5 * (lo_a[bad] + hi_a[bad])
          d_e <- abs(cand - eta_a)
          eta_a <- cand
          if(all((abs(fv) < 1e-5 & d_e < 1e-5) | ((hi_a - lo_a) < 1e-6))) break
        }
        eta_a
      }

      ci_lo_logit_h <- rep(-Inf, ngroups)
      ci_lo_fn_h <- numeric(ngroups)
      need_lo_h <- which(!best_is_0 & (ll_prof_0 < ll_targ_h))
      if(length(need_lo_h) > 0L){
        init_lo_h <- ifelse(is.finite(prof_init$ci_lower_logit[need_lo_h]),
                            prof_init$ci_lower_logit[need_lo_h],
                            e1_best[need_lo_h] - 2.0)
        res_lo_h <- solve_hier_side(need_lo_h, rep(-28, length(need_lo_h)),
                                    e1_best[need_lo_h], init_lo_h, increasing = TRUE)
        ci_lo_logit_h[need_lo_h] <- res_lo_h
        ci_lo_fn_h[need_lo_h] <- stats::plogis(res_lo_h)
      }

      ci_hi_logit_h <- rep(Inf, ngroups)
      ci_hi_fn_h <- rep(1, ngroups)
      need_hi_h <- which(!best_is_1 & (ll_prof_1 < ll_targ_h))
      if(length(need_hi_h) > 0L){
        init_hi_h <- ifelse(is.finite(prof_init$ci_upper_logit[need_hi_h]),
                            prof_init$ci_upper_logit[need_hi_h],
                            e1_best[need_hi_h] + 2.0)
        res_hi_h <- solve_hier_side(need_hi_h, e1_best[need_hi_h],
                                    rep(28, length(need_hi_h)), init_hi_h, increasing = FALSE)
        ci_hi_logit_h[need_hi_h] <- res_hi_h
        ci_hi_fn_h[need_hi_h] <- stats::plogis(res_hi_h)
      }

      if(EBayes){
        p1_u <- posterior_se_from_ci(
          est = e1_all,
          ci_lo = ci_lo_logit_h,
          ci_hi = ci_hi_logit_h,
          prior_sd = s1_all,
          se_min = 1 / sqrt(pmax(N_grp, 1))
        )
      }
    }

    lb_out <- lb_dt[grp_end, ..grp_cols]
    lb_out[, c("p1", "p2", "p1_u", "p2_u", "pold", "n") :=
             .(e1_all, e2_all, p1_u, p2_u, pold_vec[grp_end], N_grp)]
    if(compute_ci){
      lb_out[, c("ci_lower_fn", "ci_upper_fn") := .(ci_lo_fn_h, ci_hi_fn_h)]
    }
    out_parts[[length(out_parts) + 1L]] <- lb_out
  }

  res_dt <- data.table::rbindlist(out_parts, use.names = TRUE)
  data.table::setorderv(res_dt, grp_cols)
  res_dt
}


# Vectorized K-component multi-label mixture solver (ALR Newton + CLR output)
solve_general_mixture_dt <- function(dt, features_to_analyze, mutrate_design,
                                     mutcols, basecols, Poisson,
                                     samples_with_no_label, mutation_rates){
  sample <- n <- NULL
  `.` <- list

  grp_cols <- c("sample", features_to_analyze)
  dt_work <- data.table::as.data.table(dt)

  K <- nrow(mutrate_design)
  J <- ncol(mutrate_design)
  muttypes <- colnames(mutrate_design)
  md_mat <- as.matrix(mutrate_design)

  pop_labels <- vapply(seq_len(K), function(i){
    pv <- ifelse(as.logical(md_mat[i, ]), "high", "low")
    paste(paste0(pv, muttypes), collapse = "_")
  }, character(1))

  ism_cols <- paste0("invsoftmax_", pop_labels)
  frac_cols <- paste0("fraction_", pop_labels)
  se_cols <- paste0("se_invsoftmax_", pop_labels)
  ci_lo_cols <- paste0("ci_lower_fraction_", pop_labels)
  ci_hi_cols <- paste0("ci_upper_fraction_", pop_labels)
  listnames <- c(ism_cols, frac_cols, se_cols, ci_lo_cols, ci_hi_cols)

  if(nrow(dt_work) == 0L){
    empty_dt <- dt_work[, ..grp_cols]
    for(cn in listnames) empty_dt[[cn]] <- numeric(0)
    empty_dt[["n"]] <- numeric(0)
    return(dplyr::as_tibble(empty_dt))
  }

  low_rows <- which(rowSums(md_mat) == 0)
  ref_idx <- if(length(low_rows) > 0L) low_rows[1L] else K
  free_idx <- setdiff(seq_len(K), ref_idx)
  D <- length(free_idx)
  C_clr <- diag(K) - matrix(1 / K, nrow = K, ncol = K)

  all_samps <- unique(dt_work$sample)
  out_parts <- vector("list", length(all_samps))

  for(s_idx in seq_along(all_samps)){
    samp <- all_samps[s_idx]
    s_dt <- dt_work[sample == samp]

    if(samp %in% samples_with_no_label){
      nl_dt <- s_dt[, .(n = as.numeric(sum(n))), by = grp_cols]
      n_nl <- nrow(nl_dt)
      is_all_low <- rowSums(md_mat) == 0
      for(k in seq_len(K)){
        nl_dt[[ism_cols[k]]] <- rep(if(is_all_low[k]) Inf else -Inf, n_nl)
      }
      for(k in seq_len(K)){
        nl_dt[[frac_cols[k]]] <- rep(if(is_all_low[k]) 1 else 0, n_nl)
      }
      for(k in seq_len(K)){
        nl_dt[[se_cols[k]]] <- rep(0, n_nl)
      }
      for(k in seq_len(K)){
        nl_dt[[ci_lo_cols[k]]] <- rep(if(is_all_low[k]) 1 else 0, n_nl)
      }
      for(k in seq_len(K)){
        nl_dt[[ci_hi_cols[k]]] <- rep(if(is_all_low[k]) 1 else 0, n_nl)
      }
      col_ord <- c("sample", features_to_analyze, listnames, "n")
      out_parts[[s_idx]] <- nl_dt[, ..col_ord]
      next
    }

    data.table::setorderv(s_dt, grp_cols)
    grp_id <- data.table::rleidv(s_dt, cols = grp_cols)
    n_rows <- nrow(s_dt)
    grp_end <- if(n_rows > 1L){
      c(which(grp_id[-1L] != grp_id[-n_rows]), n_rows)
    }else{
      1L
    }
    ngroups <- length(grp_end)
    n_vec <- as.numeric(s_dt[["n"]])
    N_grp <- grp_sum_contig(n_vec, grp_end)

    pnew_s <- vapply(muttypes, function(nm){
      mr <- mutation_rates[[nm]]
      mr$pnew[mr$sample == samp][1L]
    }, numeric(1))
    pold_s <- vapply(muttypes, function(nm){
      mr <- mutation_rates[[nm]]
      mr$pold[mr$sample == samp][1L]
    }, numeric(1))

    # Precompute n_rows x K emission log-likelihood matrix once per sample
    E_mat <- matrix(0, nrow = n_rows, ncol = K)
    for(j in seq_len(J)){
      mj <- as.numeric(s_dt[[mutcols[j]]])
      bj <- as.numeric(s_dt[[basecols[j]]])
      if(Poisson){
        ll_hi <- stats::dpois(mj, lambda = pnew_s[j] * bj, log = TRUE)
        ll_lo <- stats::dpois(mj, lambda = pold_s[j] * bj, log = TRUE)
      }else{
        ll_hi <- stats::dbinom(mj, bj, prob = pnew_s[j], log = TRUE)
        ll_lo <- stats::dbinom(mj, bj, prob = pold_s[j], log = TRUE)
      }
      for(k in seq_len(K)){
        E_mat[, k] <- E_mat[, k] + if(as.logical(md_mat[k, j])) ll_hi else ll_lo
      }
    }

    eval_gen <- function(alpha_mat, sg, se, sn, E_sub, sN, want_hess = TRUE){
      n_g <- nrow(alpha_mat)
      full_a <- matrix(0, nrow = n_g, ncol = K)
      full_a[, free_idx] <- alpha_mat
      row_m <- do.call(pmax, unname(split(full_a, col(full_a))))
      shifted_a <- full_a - row_m
      exp_a <- exp(shifted_a)
      sum_ea <- rowSums(exp_a)
      theta <- exp_a / sum_ea
      log_theta <- shifted_a - log(sum_ea)

      Z <- log_theta[sg, , drop = FALSE] + E_sub
      z_max <- do.call(pmax, unname(split(Z, col(Z))))
      bad_z <- !is.finite(z_max)
      z_max_safe <- ifelse(bad_z, 0, z_max)
      exp_Z <- exp(Z - z_max_safe)
      sum_ez <- rowSums(exp_Z)
      ll_obs <- z_max_safe + log(sum_ez)
      ll_obs[bad_z] <- -Inf
      resp <- exp_Z / sum_ez
      if(any(bad_z)) resp[bad_z, ] <- 1 / K

      f_val <- -grp_sum_contig(sn * ll_obs, se) - 0.5 * rowSums(log_theta)

      R_free <- matrix(0, nrow = n_g, ncol = D)
      grad <- matrix(0, nrow = n_g, ncol = D)
      wt_tot <- sN + 0.5 * K
      for(d in seq_len(D)){
        k_d <- free_idx[d]
        R_free[, d] <- grp_sum_contig(sn * resp[, k_d], se)
        grad[, d] <- wt_tot * theta[, k_d] - R_free[, d] - 0.5
      }

      if(!want_hess){
        return(list(f = f_val, grad = grad))
      }

      H_arr <- array(0, dim = c(n_g, D, D))
      for(a in seq_len(D)){
        ka <- free_idx[a]
        r_a <- resp[, ka]
        sum_raa <- grp_sum_contig(sn * r_a * r_a, se)
        H_arr[, a, a] <- wt_tot * theta[, ka] * (1 - theta[, ka]) - (R_free[, a] - sum_raa)
        if(a < D){
          for(b in (a + 1L):D){
            kb <- free_idx[b]
            sum_rab <- grp_sum_contig(sn * r_a * resp[, kb], se)
            hab <- -wt_tot * theta[, ka] * theta[, kb] + sum_rab
            H_arr[, a, b] <- hab
            H_arr[, b, a] <- hab
          }
        }
      }
      list(f = f_val, grad = grad, H = H_arr)
    }

    alpha_all <- matrix(0, nrow = ngroups, ncol = D)
    active <- seq_len(ngroups)
    sub_grp <- grp_id
    grp_sizes <- diff(c(0L, grp_end))
    sub_end <- grp_end
    sub_n <- n_vec
    E_sub <- E_mat
    sub_N <- N_grp
    alpha_a <- alpha_all

    for(iter in seq_len(35L)){
      ev <- eval_gen(alpha_a, sub_grp, sub_end, sub_n, E_sub, sub_N, want_hess = TRUE)
      pinned <- (alpha_a <= -9 & ev$grad > 0) | (alpha_a >= 9 & ev$grad < 0)
      pgrad <- ifelse(pinned, 0, ev$grad)
      max_pg <- do.call(pmax, unname(split(abs(pgrad), col(pgrad))))
      conv <- max_pg < 1e-7

      if(all(conv)){
        alpha_all[active, ] <- alpha_a
        break
      }

      if(iter %% 3L == 0L && sum(conv) > 0.25 * length(conv)){
        done_loc <- which(conv)
        alpha_all[active[done_loc], ] <- alpha_a[done_loc, , drop = FALSE]

        rem_loc <- which(!conv)
        active <- active[rem_loc]
        alpha_a <- alpha_a[rem_loc, , drop = FALSE]
        sub_N <- sub_N[rem_loc]
        grp_sizes <- grp_sizes[rem_loc]

        remap <- integer(length(conv))
        remap[rem_loc] <- seq_along(rem_loc)
        rk <- remap[sub_grp] > 0L
        sub_grp <- remap[sub_grp[rk]]
        sub_end <- cumsum(grp_sizes)
        sub_n <- sub_n[rk]
        E_sub <- E_sub[rk, , drop = FALSE]

        ev$f <- ev$f[rem_loc]
        ev$grad <- ev$grad[rem_loc, , drop = FALSE]
        ev$H <- ev$H[rem_loc, , , drop = FALSE]
        pinned <- pinned[rem_loc, , drop = FALSE]
        pgrad <- pgrad[rem_loc, , drop = FALSE]
      }

      n_act <- length(active)
      step_mat <- matrix(0, nrow = n_act, ncol = D)
      for(g in seq_len(n_act)){
        Hg <- matrix(ev$H[g, , ], nrow = D, ncol = D)
        pg_g <- pgrad[g, ]
        pin_g <- pinned[g, ]
        if(any(pin_g)){
          Hg[pin_g, ] <- 0
          Hg[, pin_g] <- 0
          diag(Hg)[pin_g] <- 1
        }
        eig_vals <- eigen(Hg, symmetric = TRUE, only.values = TRUE)$values
        min_ev <- min(eig_vals)
        if(!is.finite(min_ev) || min_ev < 1e-6){
          diag(Hg) <- diag(Hg) + (1e-4 - min(0, min_ev))
        }
        st <- tryCatch(-solve(Hg, pg_g), error = function(e) -pg_g / max(sqrt(sum(pg_g^2)), 1))
        ms <- max(abs(st))
        if(ms > 2.5) st <- st * (2.5 / ms)
        step_mat[g, ] <- st
      }

      cand <- pmin(pmax(alpha_a + step_mat, -9), 9)
      lam <- rep(1, n_act)
      for(bt in seq_len(8L)){
        f_new <- eval_gen(cand, sub_grp, sub_end, sub_n, E_sub, sub_N, want_hess = FALSE)$f
        bad <- is.na(f_new) | (f_new > ev$f + 1e-10)
        if(!any(bad)) break
        lam[bad] <- lam[bad] * 0.5
        if(bt >= 5L){
          gn <- pmax(sqrt(rowSums(pgrad[bad, , drop = FALSE]^2)), 1e-8)
          cand[bad, ] <- pmin(pmax(alpha_a[bad, , drop = FALSE] -
                                     (lam[bad] / gn) * pgrad[bad, , drop = FALSE], -9), 9)
        }else{
          cand[bad, ] <- pmin(pmax(alpha_a[bad, , drop = FALSE] +
                                     lam[bad] * step_mat[bad, , drop = FALSE], -9), 9)
        }
      }
      alpha_a <- cand
      if(iter == 35L){
        alpha_all[active, ] <- alpha_a
      }
    }

    fin <- eval_gen(alpha_all, grp_id, grp_end, n_vec, E_mat, N_grp, want_hess = TRUE)
    pin_f <- (alpha_all <= -9 & fin$grad > 0) | (alpha_all >= 9 & fin$grad < 0)
    pg_f <- ifelse(pin_f, 0, fin$grad)
    max_pg_f <- do.call(pmax, unname(split(abs(pg_f), col(pg_f))))

    # Fallback to L-BFGS-B with analytical gradient on precomputed E_mat if any group unconverged
    need_fb <- which(is.na(max_pg_f) | max_pg_f > 1e-4)
    if(length(need_fb) > 0L){
      grp_start <- c(1L, grp_end[-ngroups] + 1L)
      for(g in need_fb){
        idx <- grp_start[g]:grp_end[g]
        Eg <- E_mat[idx, , drop = FALSE]
        ng <- n_vec[idx]
        Ng <- N_grp[g]
        fn_1g <- function(par){
          eval_gen(matrix(par, nrow = 1L), rep(1L, length(idx)), length(idx),
                   ng, Eg, Ng, want_hess = FALSE)$f[1L]
        }
        gr_1g <- function(par){
          as.vector(eval_gen(matrix(par, nrow = 1L), rep(1L, length(idx)), length(idx),
                             ng, Eg, Ng, want_hess = FALSE)$grad)
        }
        opt_g <- stats::optim(par = rep(0, D), fn = fn_1g, gr = gr_1g,
                              method = "L-BFGS-B", lower = rep(-9, D), upper = rep(9, D))
        alpha_all[g, ] <- opt_g$par
      }
      fin <- eval_gen(alpha_all, grp_id, grp_end, n_vec, E_mat, N_grp, want_hess = TRUE)
    }

    full_alr <- matrix(0, nrow = ngroups, ncol = K)
    full_alr[, free_idx] <- alpha_all
    clr_mat <- full_alr - rowMeans(full_alr)

    row_m_c <- do.call(pmax, unname(split(clr_mat, col(clr_mat))))
    exp_c <- exp(clr_mat - row_m_c)
    frac_mat <- exp_c / rowSums(exp_c)

    se_mat <- matrix(0, nrow = ngroups, ncol = K)
    for(g in seq_len(ngroups)){
      Hg <- matrix(fin$H[g, , ], nrow = D, ncol = D)
      se_g <- tryCatch({
        cov_alr <- solve(Hg)
        cov_full <- matrix(0, nrow = K, ncol = K)
        cov_full[free_idx, free_idx] <- cov_alr
        cov_clr <- C_clr %*% cov_full %*% t(C_clr)
        sqrt(pmax(0, diag(cov_clr)))
      }, error = function(e) rep(Inf, K))
      se_mat[g, ] <- se_g
    }

    # Vectorized 95% Wilks profile-likelihood CIs for each population along complement simplex rays
    ci_lo_mat <- matrix(0, nrow = ngroups, ncol = K)
    ci_hi_mat <- matrix(1, nrow = ngroups, ncol = K)
    for(k in seq_len(K)){
      other_idx <- setdiff(seq_len(K), k)
      rem_mass <- rowSums(frac_mat[, other_idx, drop = FALSE])
      q_mat <- frac_mat[, other_idx, drop = FALSE] / ifelse(rem_mass > 1e-12, rem_mass, 1)
      if(any(rem_mass <= 1e-12)){
        q_mat[rem_mass <= 1e-12, ] <- 1 / length(other_idx)
      }
      log_q <- log(pmax(q_mat, 1e-300))
      W_oth <- log_q[grp_id, , drop = FALSE] + E_mat[, other_idx, drop = FALSE]
      w_mx <- do.call(pmax, unname(split(W_oth, col(W_oth))))
      w_mx_s <- ifelse(is.finite(w_mx), w_mx, 0)
      E_mk <- w_mx_s + log(rowSums(exp(W_oth - w_mx_s)))
      E_mk[!is.finite(w_mx)] <- -Inf
      pr_k <- profile_2comp_delta_vec(grp_id, grp_end, n_vec, E_mat[, k] - E_mk, N_grp = N_grp)
      ci_lo_mat[, k] <- pr_k$ci_lower_fn
      ci_hi_mat[, k] <- pr_k$ci_upper_fn
    }

    s_out <- s_dt[grp_end, ..grp_cols]
    for(k in seq_len(K)) s_out[[ism_cols[k]]] <- clr_mat[, k]
    for(k in seq_len(K)) s_out[[frac_cols[k]]] <- frac_mat[, k]
    for(k in seq_len(K)) s_out[[se_cols[k]]] <- se_mat[, k]
    for(k in seq_len(K)) s_out[[ci_lo_cols[k]]] <- ci_lo_mat[, k]
    for(k in seq_len(K)) s_out[[ci_hi_cols[k]]] <- ci_hi_mat[, k]
    s_out[["n"]] <- N_grp
    out_parts[[s_idx]] <- s_out
  }

  res_dt <- data.table::rbindlist(out_parts, use.names = TRUE)
  data.table::setorderv(res_dt, grp_cols)
  dplyr::as_tibble(res_dt)
}



