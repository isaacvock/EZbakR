#' Generate a `fraction_design` table for `EstimateFractions`
#'
#' @param mutrate_populations Character vector of the set of mutational populations
#' present in your data. For example, s4U fed data with standard nucleotide recoding
#' chemistry (e.g., TimeLapse, SLAM, TUC, AMUC, etc.) would have a `mutrate_populations`
#' of c("TC"). Dual labeling experiments with s4U and s6G feeds would have a `mutrate_populations`
#' of c("TC", "GA").
#' @return A `fraction_design` table that assumes that every possible combination of
#' mutational populations listed in `mutrate_populations` are present in your data.
#' The `present` column can be modified if this assumption is incorrect.
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




#' Estimate fractions of each RNA population
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
#' described in [Courvan et al., 2022](https://academic.oup.com/nar/article/50/19/e110/6677324). In this
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
#' @param split_multi_features If a set of reads maps ambiguously to multiple features,
#' should data for such reads be copied for each feature in the ambiguous set? If this is
#' `TRUE`, then `multi_feature_cols` also must be set. Examples where this should be set
#' to `TRUE` includes when analyzing exonic bins (concept defined in original DEXSeq paper),
#' exon-exon junctions, etc.
#' @param multi_feature_cols Character vector of columns that have the potential to
#' include assignment to multiple features. Only these columns will have their features split
#' if `split_multi_features` is `TRUE`.
#' @param multi_feature_sep String representing how ambiguous feature assignments are
#' distinguished in the feature names. For example, the default value of "+" denotes
#' that if a read maps to multiple features (call them featureA and featureB, for example),
#' then the feature column will have a value of "featureA+featureB" for that read.
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
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
EstimateFractions <- function(obj, features = "all",
                              mutrate_populations = "all",
                              fraction_design = NULL,
                              Poisson = TRUE,
                              strategy = c("standard", "hierarchical"),
                              filter_cols = "all",
                              filter_condition = `&`,
                              remove_features = c("NA", "__no_feature"),
                              split_multi_features = FALSE,
                              multi_feature_cols = NULL,
                              multi_feature_sep = "+",
                              pnew_prior_mean = -2.94,
                              pnew_prior_sd = 0.3,
                              pold_prior_mean = -6.5,
                              pold_prior_sd = 0.5,
                              hier_readcutoff = 300,
                              init_pnew_prior_sd = 0.8,
                              pnew_prior_sd_min = 0.01,
                              pnew_prior_sd_max = 0.15,
                              pold_est = NULL,
                              pold_from_nolabel = FALSE,
                              grouping_factors = NULL,
                              character_limit = 20,
                              overwrite = TRUE){


  UseMethod("EstimateFractions")


}



#' Estimate mutation rates
#'
#' Two component mixture models are fit to all data to estimate global high and
#' low mutation rates for all samples. Estimation of these mutation rates
#' are regularized through the use of informative priors whose parameters can
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
#' @import data.table
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



#' Estimate fractions of each RNA population with standard EZbakRData object
#'
#' This is the default fraction estimation strategy for EZbakR. cB is assumed to
#' be an in-memory table and `obj` is expected to be an `EZbakRData` object
#' created with `EZbakRData()`. For larger than RAM cB files, see the relevant
#' `EZbakRArrowData` method.
#'
#' @param obj `EZbakRData` object
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
#' described in [Courvan et al., 2022](https://academic.oup.com/nar/article/50/19/e110/6677324). In this
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
#' @param split_multi_features If a set of reads maps ambiguously to multiple features,
#' should data for such reads be copied for each feature in the ambiguous set? If this is
#' `TRUE`, then `multi_feature_cols` also must be set. Examples where this should be set
#' to `TRUE` includes when analyzing exonic bins (concept defined in original DEXSeq paper),
#' exon-exon junctions, etc.
#' @param multi_feature_cols Character vector of columns that have the potential to
#' include assignment to multiple features. Only these columns will have their features split
#' if `split_multi_features` is `TRUE`.
#' @param multi_feature_sep String representing how ambiguous feature assignments are
#' distinguished in the feature names. For example, the default value of "+" denotes
#' that if a read maps to multiple features (call them featureA and featureB, for example),
#' then the feature column will have a value of "featureA+featureB" for that read.
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
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
EstimateFractions.EZbakRData <- function(obj, features = "all",
                                         mutrate_populations = "all",
                                         fraction_design = NULL,
                                         Poisson = TRUE,
                                         strategy = c("standard", "hierarchical"),
                                         filter_cols = "all",
                                         filter_condition = `&`,
                                         remove_features = c("NA", "__no_feature"),
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
                                         character_limit = 20,
                                         overwrite = TRUE){

  ### Hack to deal with devtools::check() NOTEs
  n <- params <- coverage <- pold <- pnew <- lpnew_uncert <- pnew_prior <- ..colvect <- NULL
  fit <- mixture_fit <- present <- NULL

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

  if(mutrate_populations == "all"){

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

  cB <- setDT(cB)

  # Keep only feature of interest
  cols_to_group <- c(necessary_basecounts, pops_to_analyze, "sample", features_to_analyze)
  cB <- cB[sample %in% metadf$sample][,.(n = sum(n)), by = cols_to_group]

  ### Filter out columns not mapping to any feature (easier and faster in data.table)
  if(filter_cols[1] == "all" & length(filter_cols) == 1){

    filter_cols <- features_to_analyze

  }

  cB <- cB[cB[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ]


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

    cB <- split_features(cB, multi_feature_cols,  c(cols_to_group, necessary_basecounts))

  }

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

      feature_specific <- cB[
        coverages, nomatch = NULL
      ][,
        .(params = list(fit_tcmm(muts = get(pops_to_analyze),
                                      nucs = get(necessary_basecounts),
                                      n = n,
                                      Poisson = Poisson,
                                      pold = unique(pold),
                                      pnew_prior_mean = unique(pnew),
                                 fraction_prior_sd = 1,
                                      pnew_prior_sd = init_pnew_prior_sd)),
          pold = unique(pold)),
        by = c("sample", features_to_analyze)
      ]


      feature_specific[, c("pnew", "lpnew_uncert") := .( inv_logit(sapply(params, `[[`, 2)),
                                         sapply(params, `[[`, 4))]


      ### Hierarchical model using priors inferred from feature-specific pnew distribution

      message("Estimating fractions with feature-specific pnews")

      ## Mean should really be the estimated sample wide value
      pnew_prior_df <- mutrates %>%
        dplyr::rename(pnew_prior = pnew) %>%
        dplyr::select(sample, pnew_prior)
      global_est <- feature_specific[pnew_prior_df][pnew < 0.3 & pnew > max(cB$pold, na.rm = TRUE)][,.(pnew_prior = mean(logit(pnew_prior)),
                                    pnew_prior_sd = stats::sd(logit(pnew)) - mean(lpnew_uncert),
                                    pold = mean(pold)),
                                 by = sample]

      global_est[, pnew_prior_sd := ifelse(pnew_prior_sd < 0,
                                           pnew_prior_sd_min,
                                           ifelse(pnew_prior_sd > pnew_prior_sd_max,
                                                  pnew_prior_sd_max,
                                                  pnew_prior_sd))
      ]

      global_est[, pnew_prior_sd := min(pnew_prior_sd)]

      setkey(global_est, sample)
      setkey(cB, sample)

      feature_specific <- global_est[cB][,
         .(params = ifelse(!(unique(sample) %in% samples_with_no_label), list(fit_tcmm(muts = get(pops_to_analyze),
                                                                                   nucs = get(necessary_basecounts),
                                                                                   n = n,
                                                                                   Poisson = Poisson,
                                                                                   pold = unique(pold),
                                                                                   pnew_prior_mean = unique(pnew_prior),
                                                                                   pnew_prior_sd = unique(pnew_prior_sd))),
                                      list(list(p1 = -Inf,
                                                                                            p2 = -Inf,
                                                                                            p1_u = 0,
                                                                                            p2_u = 0))),
            n = sum(n)),
          by = c("sample", features_to_analyze)
      ]



      feature_specific[, c(col_name, "pnew",
                           uncertainty_col) := .(sapply(params, `[[`, 1),
                                               inv_logit(sapply(params, `[[`, 2)),
                                               sapply(params, `[[`, 3))]


      fns <- dplyr::as_tibble(feature_specific) %>%
        dplyr::mutate(!!natural_col_name := inv_logit(!!dplyr::sym(col_name))) %>%
        dplyr::select(sample, !!features_to_analyze, !!natural_col_name,
                      !!col_name, !!uncertainty_col, n)


      ### Save feature-specific pnews

      mutrates[, pnew := NULL]

      colvect <- c("sample", features_to_analyze, "pnew")

      feature_mutrates <- feature_specific[,..colvect][
             mutrates, on = "sample", nomatch = NULL
        ]

      current_name <- names(obj$mutation_rates)[1]
      obj$mutation_rates[[paste0(paste(features_to_analyze, collapse = "_"), "_", current_name)]] <- feature_mutrates



    }else{

      fns <- dplyr::as_tibble(cB) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", features_to_analyze)))) %>%
        dplyr::summarise(fit = ifelse(!(unique(sample) %in% samples_with_no_label),
                                      list(I(stats::optim(0,
                                                   fn = tcmml,
                                                   muts = !!dplyr::sym(pops_to_analyze),
                                                   nucs = !!dplyr::sym(necessary_basecounts),
                                                   pnew = pnew,
                                                   pold = pold,
                                                   pnew_prior_mean = pnew_prior_mean,
                                                   pnew_prior_sd = pnew_prior_sd,
                                                   pold_prior_mean = pold_prior_mean,
                                                   pold_prior_sd = pold_prior_sd,
                                                   Poisson = Poisson,
                                                   n = n,
                                                   lower = -9,
                                                   upper = 9,
                                                   method = "L-BFGS-B",
                                                   hessian = TRUE))),
                                      list(I(list(par = -Inf,
                                                  hessian = 1)))),
                         n = sum(n)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          !!col_name := purrr::map_dbl(fit, ~ .x$par[1]),
          !!uncertainty_col := purrr::map_dbl(fit, ~ sqrt(solve(.x$hessian)[1]))
        ) %>%
        dplyr::mutate(
          !!uncertainty_col := ifelse(
            !!dplyr::sym(col_name) == -Inf,
            0,
            !!dplyr::sym(uncertainty_col)
          )
        ) %>%
        dplyr::select(-fit) %>%
        dplyr::mutate(!!natural_col_name := inv_logit(!!dplyr::sym(col_name))) %>%
        dplyr::select(sample, !!features_to_analyze, !!natural_col_name,
                      !!col_name, !!uncertainty_col, n)

    }





  }else{


    # pops_to_analyze needs to match order as it appears in mutrate_design
    # and necessary_basecounts order needs to match relevant nucleotide as well
    pops_to_analyze <- colnames(mutrate_design)
    necessary_basecounts <- paste0("n", substr(pops_to_analyze, start = 1, stop = 1))

    fns <- cB[,.(mixture_fit = list(fit_general_mixture(dataset = .SD,
                                                        mutrate_design = mutrate_design,
                                                        mutcols = pops_to_analyze,
                                                        basecols = necessary_basecounts,
                                                        Poisson = Poisson,
                                                        twocomp = FALSE,
                                                        sample = unique(sample),
                                                        samples_with_no_label = samples_with_no_label,
                                                        pnew = sapply(pops_to_analyze,
                                                                      function(name) mutation_rates[[name]]$pnew[mutation_rates[[name]]$sample == sample ]),
                                                        pold = sapply(pops_to_analyze,
                                                                      function(name) mutation_rates[[name]]$pold[mutation_rates[[name]]$sample == sample ]) ) ),
                 n = sum(n)),
              by = c("sample", features_to_analyze)]


    # Unroll the fraction estimates
    fns <- fns %>%
      tidyr::unnest_wider(mixture_fit)

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
  if(length(obj[['metadata']]) > 0){

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
                                                            repeatID = repeatID)



  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)

}




#' Estimate mutation rates
#'
#' Two component mixture models are fit to all data to estimate global high and
#' low mutation rates for all samples. Estimation of these mutation rates
#' are regularized through the use of informative priors whose parameters can
#' be altered using the arguments defined below. This is the default method
#' which expects the input `obj` to include an in-memory cB table. For analyses
#' of larger than RAM datasets, see the `EZbakRArrowData` method.
#'
#' @param obj An `EZbakRData` object
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
#' @import data.table
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
  cB <- data.table::setDT(cB)




  # Figure out which mutation counts are in the cB
  mutcounts_in_cB <- find_mutcounts(obj)

  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # Which populations to analyze?
  if(populations == "all"){

    muts_analyze <- mutcounts_in_cB

  }else{

    muts_analyze <- populations

    if(!(mutcounts_in_cB %in% muts_analyze)){
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



#' Estimate fractions of each RNA population using Apache Arrow backend.
#'
#' This is an alternative to the default `EstimateFractions` method that can help
#' with analyses of larger than RAM datasets. The provided "cB" is expected to be
#' an on-disk Arrow Dataset. Furthermore, it is expected to be partitioned by the
#' sample name, which will allow this method to read only a single-sample worth
#' of data into memory at a time. This can significantly reduce RAM usage. Input
#' object should be created with `EZbakRArrowData()`.
#'
#' @param obj `EZbakRArrowData` object
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
#' described in [Courvan et al., 2022](https://academic.oup.com/nar/article/50/19/e110/6677324). In this
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
#' @param split_multi_features If a set of reads maps ambiguously to multiple features,
#' should data for such reads be copied for each feature in the ambiguous set? If this is
#' `TRUE`, then `multi_feature_cols` also must be set. Examples where this should be set
#' to `TRUE` includes when analyzing exonic bins (concept defined in original DEXSeq paper),
#' exon-exon junctions, etc.
#' @param multi_feature_cols Character vector of columns that have the potential to
#' include assignment to multiple features. Only these columns will have their features split
#' if `split_multi_features` is `TRUE`.
#' @param multi_feature_sep String representing how ambiguous feature assignments are
#' distinguished in the feature names. For example, the default value of "+" denotes
#' that if a read maps to multiple features (call them featureA and featureB, for example),
#' then the feature column will have a value of "featureA+featureB" for that read.
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
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
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
                                              split_multi_features = FALSE,
                                              multi_feature_cols = NULL,
                                              multi_feature_sep = "+",
                                              pnew_prior_mean = -2.94,
                                              pnew_prior_sd = 0.3,
                                              pold_prior_mean = -6.5,
                                              pold_prior_sd = 0.5,
                                              hier_readcutoff = 300,
                                              init_pnew_prior_sd = 0.8,
                                              pnew_prior_sd_min = 0.01,
                                              pnew_prior_sd_max = 0.15,
                                              pold_est = NULL,
                                              pold_from_nolabel = FALSE,
                                              grouping_factors = NULL,
                                              character_limit = 20,
                                              overwrite = TRUE){


  ### Hack to deal with devtools::check() NOTEs
  present <- n <- reads <- params <- coverage <- pold <- pnew <- sd <- lpnew_uncert <- pnew_prior <- ..colvect <- NULL
  fit <- mixture_fit <- NULL

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
  cB <- obj$cBds
  cBschema <- arrow::schema(cB)
  cB_cols <- names(cBschema)

  ### Vectors of potential column names

  # Mutation counts possible
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]

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

  if(mutrate_populations == "all"){

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


  if(filter_cols[1] == "all" & length(filter_cols) == 1){

    filter_cols <- features_to_analyze

  }

  if(!all(filter_cols %in% features_to_analyze)){

    stop("Some of your filter_cols are not included in the features you want to analyze!")

  }


  for(s in seq_along(all_samples)){


    message(paste0("ANALYZING ", all_samples[s], "..."))


    if(all_samples[s] %in% samples_without_label){

      message("Counting reads")

      ctl_cols_to_group <- c("sample", features_to_analyze)

      message()

      sample_fns <- cB %>%
        dplyr::filter(sample == all_samples[s]) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ctl_cols_to_group))) %>%
        dplyr::summarise(n = sum(n)) %>%
        dplyr::collect() %>%
        dplyr::ungroup()

      ### Filter out columns not mapping to any feature (easier and faster in data.table)
      sample_fns <- setDT(sample_fns)
      sample_fns <- dplyr::as_tibble(sample_fns[sample_fns[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ])

      ### Split multi feature mappers if necessary
      if(split_multi_features){

        message("Splitting multi-feature mapping reads")


        sample_fns <- split_features(sample_fns, multi_feature_cols, ctl_cols_to_group)

      }


      dplyr::as_tibble(sample_fns) %>%
        dplyr::mutate(!!natural_col_name := 0,
                      !!col_name := -Inf,
                      !!uncertainty_col := 0)



    }else{



      if(Poisson){

        message("Extracting data for sample of interest and summarize out nucleotide content")


        sample_cB <- cB %>%
          dplyr::filter(sample == all_samples[s]) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
          dplyr::summarise(reads = sum(n),
                           !!necessary_basecounts := sum(!!dplyr::sym(necessary_basecounts)*n)/sum(n)) %>%
          dplyr::collect() %>%
          dplyr::ungroup() %>%
          dplyr::rename(n = reads)

      }else{

        message("Extracting data for sample of interest")

        sample_cB <- cB %>%
          dplyr::filter(sample == all_samples[s]) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
          dplyr::summarise(n = sum(n)) %>%
          dplyr::collect()  %>%
          dplyr::ungroup()
      }


      ### Filter out columns not mapping to any feature (easier and faster in data.table)
      sample_cB <- setDT(sample_cB)
      sample_cB <- dplyr::as_tibble(sample_cB[sample_cB[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ])

      ### Split multi feature mappers if necessary
      if(split_multi_features){

        message("Splitting multi-feature mapping reads")

        sample_cB <- split_features(sample_cB, multi_feature_cols, c(cols_to_group, necessary_basecounts))


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

          feature_specific <- sample_cB[
            coverages, nomatch = NULL
          ][,
            .(params = list(fit_tcmm(muts = get(pops_to_analyze),
                                     nucs = get(necessary_basecounts),
                                     n = n,
                                     Poisson = Poisson,
                                     pold = unique(pold),
                                     pnew_prior_mean = unique(pnew),
                                     pnew_prior_sd = init_pnew_prior_sd)),
              pold = unique(pold)),
            by = c("sample", features_to_analyze)
          ]


          feature_specific[, c("pnew", "lpnew_uncert") := .( inv_logit(sapply(params, `[[`, 2)),
                                                             sapply(params, `[[`, 4))]


          ### Hierarchical model

          message("Estimating fractions with feature-specific pnews")

          global_est <- feature_specific[,.(pnew_prior = mean(logit(pnew)),
                                            pnew_prior_sd = stats::sd(logit(pnew)) - mean(lpnew_uncert),
                                            pold = mean(pold)),
                                         by = sample]

          global_est[, pnew_prior_sd := ifelse(pnew_prior_sd < 0,
                                               pnew_prior_sd_min,
                                               ifelse(pnew_prior_sd > pnew_prior_sd_max,
                                                      pnew_prior_sd_max,
                                                      pnew_prior_sd))
          ]

          global_est[, pnew_prior_sd := min(pnew_prior_sd)]

          setkey(global_est, sample)
          setkey(sample_cB, sample)

          feature_specific <- sample_cB[global_est, nomatch = NULL][,
                                                             .(params = dplyr::case_when(
                                                               !(unique(sample) %in% samples_with_no_label) ~ list(fit_tcmm(muts = get(pops_to_analyze),
                                                                                                                            nucs = get(necessary_basecounts),
                                                                                                                            n = n,
                                                                                                                            Poisson = Poisson,
                                                                                                                            pold = unique(pold),
                                                                                                                            pnew_prior_mean = unique(pnew_prior),
                                                                                                                            pnew_prior_sd = unique(pnew_prior_sd))),
                                                               .default = list(list(p1 = -Inf,
                                                                                    p2 = -Inf,
                                                                                    p1_u = 0,
                                                                                    p2_u = 0))),
                                                               n = sum(n)),
                                                             by = c("sample", features_to_analyze)
          ]


          feature_specific[, c(col_name, "pnew",
                               uncertainty_col) := .(sapply(params, `[[`, 1),
                                                     inv_logit(sapply(params, `[[`, 2)),
                                                     sapply(params, `[[`, 3))]



          sample_fns <- dplyr::as_tibble(feature_specific) %>%
            dplyr::mutate(!!natural_col_name := inv_logit(!!dplyr::sym(col_name))) %>%
            dplyr::select(sample, !!features_to_analyze, !!natural_col_name,
                          !!col_name, !!uncertainty_col, n)


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

          # Kinda cutesy, but works because bind_rows(dataframe, NULL) = dataframe
          # so don't need a special condition for s == 1.
          obj$mutation_rates[[paste0(paste(features_to_analyze, collapse = "_"), "_", current_name)]] <- dplyr::bind_rows(feature_mutrates,
                                                                                                                          obj$mutation_rates[[paste0(paste(features_to_analyze, collapse = "_"), "_", current_name)]])


        }else{

          ##### WHOLE CODE BLOCK NEARLY IDENTICAL TO AS IN DEFAULT METHOD
          # NOTE: Throws warning sometimes due to summarise returning multi-entries
          sample_fns <- dplyr::as_tibble(sample_cB) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", features_to_analyze)))) %>%
            dplyr::summarise(fit = ifelse(!(unique(sample) %in% samples_without_label),
                                          list(I(stats::optim(0,
                                                       fn = tcmml,
                                                       muts = !!dplyr::sym(pops_to_analyze),
                                                       nucs = !!dplyr::sym(necessary_basecounts),
                                                       pnew = pnew,
                                                       pold = pold,
                                                       pnew_prior_mean = pnew_prior_mean,
                                                       pnew_prior_sd = pnew_prior_sd,
                                                       pold_prior_mean = pold_prior_mean,
                                                       pold_prior_sd = pold_prior_sd,
                                                       Poisson = Poisson,
                                                       n = n,
                                                       lower = -9,
                                                       upper = 9,
                                                       method = "L-BFGS-B",
                                                       hessian = TRUE))),
                                          list(I(list(par = -Inf,
                                                      hessian = Inf)))),
                             n = sum(n)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
              !!col_name := purrr::map_dbl(fit, ~ .x$par[1]),
              !!uncertainty_col := purrr::map_dbl(fit, ~ sqrt(solve(.x$hessian)[1]))
            ) %>%
            dplyr::select(-fit) %>%
            dplyr::mutate(!!natural_col_name := inv_logit(!!dplyr::sym(col_name))) %>%
            dplyr::select(sample, !!features_to_analyze, !!natural_col_name,
                          !!col_name, !!uncertainty_col, n)


        }


      }else{

        Poisson <- FALSE


        sample_fns <- sample_cB[,.(mixture_fit = list(fit_general_mixture(dataset = .SD,
                                                                          mutrate_design = mutrate_design,
                                                                          mutcols = pops_to_analyze,
                                                                          basecols = necessary_basecounts,
                                                                          Poisson = Poisson,
                                                                          twocomp = FALSE,
                                                                          pnew = sapply(pops_to_analyze,
                                                                                        function(name) mutation_rates[[name]]$pnew[mutation_rates[[name]]$sample == sample ]),
                                                                          pold = sapply(pops_to_analyze,
                                                                                        function(name) mutation_rates[[name]]$pold[mutation_rates[[name]]$sample == sample ]) ) ),
                                   n = sum(n)),
                                by = c("sample", features_to_analyze)]


        # Unroll the fraction estimates
        sample_fns <- sample_fns %>%
          tidyr::unnest_wider(mixture_fit)

      }



    }




    fns[[s]] <- sample_fns

  }


  fns <- dplyr::bind_rows(fns)


  message("Processing output")



  ##### REST OF THIS CODE IS IDENTICAL TO AS IN DEFAULT METHOD


  # What should output be named?
  fraction_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")

  if(nchar(fraction_vect) > character_limit){

    num_fractions <- length(obj[['fractions']])
    fraction_vect <- paste0("fractions", num_fractions + 1)

  }

  # Are there any metadata or fractions objects at this point?
  if(length(obj[['metadata']]) > 0){

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
                                                            repeatID = repeatID)



  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)



}



#' Estimate mutation rates
#'
#' Two component mixture models are fit to all data to estimate global high and
#' low mutation rates for all samples. Estimation of these mutation rates
#' are regularized through the use of informative priors whose parameters can
#' be altered using the arguments defined
#'
#' This method expects the input object to contain an on-disk Arrow Dataset in place
#' of an in-memory cB file. Furtheromore, this dataset is expected to be partitioned by
#' the individual samples. This allows EZbakR to load only a single sample worth of data
#' in to memory at a time, which can significantly reduce RAM usage. The input object
#' should be created with `EZbakRArrowData()`.
#'
#' @param obj An `EZbakRData` object
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

  ### Hack to deal with devtools::check() NOTEs
  n <- params <- p1 <- p2 <- ..cols_to_keep <- muts <- nucs <- NULL

  rm(..cols_to_keep)

  `.` <- list

  ### Figure out which features will be analyzed

  # cB columns
  cB <- obj$cBds
  cBschema <- arrow::schema(cB)
  cB_cols <- names(cBschema)

  ### Vectors of potential column names

  # Mutation counts possible
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]

  # Which mutcounts are in the cB?
  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]

  # Base count columns in cB
  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))


  # Which populations to analyze?
  if(populations == "all"){

    muts_analyze <- mutcounts_in_cB

  }else{

    muts_analyze <- populations

    if(!(muts_analyze %in% mutcounts_in_cB)){
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

      ### Compile all -label data
      for(c in seq_along(samples_with_no_label)){


        ctl_cB <- cB %>%
          dplyr::filter(sample == samples_with_no_label[c]) %>%
          dplyr::group_by(sample) %>%
          dplyr::summarise(muts = sum(n*!!dplyr::sym(muts_analyze[i])),
                           nucs = sum(n*!!dplyr::sym(nucs_analyze[i]))) %>%
          dplyr::collect() %>%
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
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(n = sum(n)) %>%
        dplyr::collect()


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


# Split up multi-feature sets
split_features <- function(cB, multi_feature_cols, grouping_cols){

  ### Hack to deal with devtools::check() NOTEs
  n <- NULL

  `.` <- list

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
    tidyr::separate_rows(dplyr::all_of(new_temp_cols), sep = "\\+")

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

  cB <- obj$cB
  cB <- data.table::setDT(cB)

  # Mutation count column names
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]



  # Which mutcounts are in the cB?
  cB_cols <- colnames(cB)


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

    ll <- fn*stats::dpois(muts, nucs*pnew) + (1 - fn)*stats::dpois(muts, nucs*pold)

  }else{

    ll <- fn*stats::dbinom(muts, nucs, pnew) + (1 - fn)*stats::dbinom(muts, nucs, pold)

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




  ll <- -(sum(n*log(ll)) + prior)

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

    likelihoods <- inv_logit(param[1])*stats::dpois(muts,
                                             lambda = pnew*nucs ) +

      (1 - inv_logit(param[1]))*stats::dpois(muts,
                                      lambda = pold*nucs )

  }else{

    likelihoods <- inv_logit(param[1])*stats::dbinom(x = muts,
                                              size = nucs,
                                              prob = pnew) +

      (1 - inv_logit(param[1]))*stats::dbinom(x = muts,
                                       size = nucs,
                                       prob = pold)


  }

  return(-sum(n*log(likelihoods)) - stats::dnorm(param[1], log = TRUE))


}

# Abstract the concept of an NR-seq mixture model to oblivion
generalized_likelihood <- function(param, dataset, Poisson = TRUE,
                                   mutrate_design, pnew, pold, highpop,
                                   mutcols, basecols, twocomp = TRUE){


  # # # POTENTIAL VECTORIZATION
  # # # Assuming `rates` is KxN, `dataset` is a dataframe, `basecols` and `mutcols` are vectors of column names
  # # # Precompute rates for each population and mutation
  # rates_new <- sweep(mutrate_design, 2, pnew, `*`)
  # rates_old <- sweep(1 - mutrate_design, 2, pold, `*`)
  # rates <- rates_new + rates_old
  #
  # # Assuming `rates` is KxN, `dataset` is a dataframe, `basecols` and `mutcols` are vectors of column names
  #
  # # Step 1: Generate lambda matrix for a single data point (extendable to all points)
  # generate_lambda <- function(dataset_row) {
  #   sapply(1:length(basecols), function(j) rates[, j] * dataset_row[basecols[j]])
  # }
  #
  # # For all dataset points (assuming dataset is a dataframe and iterating over its rows)
  # lambda_matrix_list <- apply(dataset, 1, generate_lambda)
  #
  # # Step 2: Calculate Poisson likelihoods for each lambda
  # calculate_poisson_likelihoods <- function(lambda_matrix, dataset_row) {
  #   sapply(1:nrow(lambda_matrix), function(i) {
  #     dpois(dataset_row[mutcols], lambda_matrix[i, ])
  #   })
  # }
  #
  # # For all dataset points
  # poisson_likelihoods_list <- mapply(calculate_poisson_likelihoods, lambda_matrix_list,
  #                                    MoreArgs = list(dataset_row = dataset))
  #
  # # Step 3: Calculate total likelihood for each data point
  # calculate_total_likelihood <- function(poisson_likelihoods) {
  #   colSums(apply(poisson_likelihoods, 2, prod))
  # }
  #
  # total_likelihoods <- lapply(poisson_likelihoods_list, calculate_total_likelihood)
  #
  # # Assuming you need to sum these total likelihoods for a final value (depends on your model)
  # final_likelihood <- sum(unlist(total_likelihoods))



  likelihoods <- rep(0, times = length(dataset[[1]]))

  # Fractions of populations
  fractions <- softmax(param)


  # Get full likelihood for generalized mutational models
  for(i in 1:nrow(mutrate_design)){

    L <- 1

    for(j in 1:ncol(mutrate_design)){


      if(Poisson){

        # data.table syntax is an act of violence
        L <- L*stats::dpois(dataset[[mutcols[j]]],
                      lambda = (as.numeric(mutrate_design[i, j])*pnew[j] +
                        pold[j]*(1 - as.numeric(mutrate_design[i, j])))*dataset[[basecols[j]]] )


      }else{

        # data.table syntax is an act of violence
        L <- L*stats::dbinom(dataset[[mutcols[j]]],
                      dataset[[basecols[j]]],
                      prob = as.numeric(mutrate_design[i, j])*pnew[j] +
                        pold[j]*(1 - as.numeric(mutrate_design[i, j])))

      }



    }

    if(twocomp){

      if(i == highpop){

        likelihoods <- likelihoods + L*inv_logit(param[1])

      }else{

        likelihoods <- likelihoods + L*(1 - inv_logit(param[1]))

      }


    }else{

      likelihoods <- likelihoods + L*fractions[i]


    }


  }

  # Prior regularize
  if(twocomp){

    return(-sum(dataset[["n"]]*log(likelihoods)) - stats::dnorm(param[1], log = TRUE))

  }else{

    return(-sum(dataset[["n"]]*log(likelihoods)))

  }

}

# Softmax function
softmax <- function(vect){

  return(exp(vect)/sum(exp(vect)))

}

# Fit and process the output of the generalized likelihood model
fit_general_mixture <- function(dataset, Poisson = TRUE, mutrate_design, twocomp = FALSE,
                                pnew, pold, mutcols, basecols, highpop = NULL,
                                samples_with_no_label = NULL,
                                sample = NULL){

  # Hack to deal with devtools::check()
  fsvector <- NULL


  if(sample %in% samples_with_no_label){

    # Figure out what to call each of the fractions
    # Probably "p""muttype""old/new"_"muttype""old/new"
    muttypes <- colnames(mutrate_design)


    # Get vector of population statuses (new and old) for naming the output list
    population_list <- vector(mode = "list", length = nrow(mutrate_design))
    fvector <- rep(0, times = nrow(mutrate_design))
    ismvector <- fsvector
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

      if(all(population_vector) == "low"){
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
                              "se_invsoftmax_"),
                            each = nrow(mutrate_design)),
                        listnames)

    outlist <- c(as.list(ismvector),
                 as.list(fvector),
                 as.list(rep(0, times = nrow(mutrate_design))) )



    names(outlist) <- listnames

  }else{

    pnew <- pnew[match(colnames(mutrate_design), names(pnew))]
    pold <- pold[match(colnames(mutrate_design), names(pold))]


    if(twocomp){

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

      outlist <- list(fit$par, logit(1 - inv_logit(fit$par)))

      names(outlist) <- c(paste0(mutcols, "high"),
                          paste0(mutcols, "low"))


    }else{

      fit <- stats::optim(par = rep(0, times = nrow(mutrate_design)),
                          fn = generalized_likelihood,
                          dataset = dataset,
                          Poisson = Poisson,
                          mutrate_design = mutrate_design,
                          pnew = pnew,
                          pold = pold,
                          twocomp = twocomp,
                          mutcols = mutcols,
                          basecols = basecols,
                          method = "L-BFGS-B",
                          hessian = TRUE)

      # Figure out what to call each of the fractions
      # Probably "p""muttype""old/new"_"muttype""old/new"
      muttypes <- colnames(mutrate_design)


      # Get vector of population statuses (new and old) for naming the output list
      population_list <- vector(mode = "list", length = nrow(mutrate_design))
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

        population_list[[i]] <- population_vector

      }


      # Make names for the output list that are easily interpretable (I hope)
      listnames <- lapply(population_list, function(x) paste0(x, muttypes))
      listnames <- sapply(listnames, function(x) paste(x, collapse = "_"),
                          simplify = "vector")

      listnames <- paste0(rep(c("invsoftmax_",
                                "fraction_",
                                "se_invsoftmax_"),
                              each = length(fit$par)),
                          listnames)


      uncertainties <- tryCatch(
        {
          sqrt(abs(diag(solve(fit$hessian))))
        },
        error = function(e) {
          rep(Inf, times = length(fit$par))
        }
      )

      outlist <- c(as.list(fit$par),
                   as.list(softmax(fit$par)),
                   as.list(uncertainties) )


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


}


