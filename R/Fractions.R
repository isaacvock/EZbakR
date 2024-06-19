#' Generate a `fraction_design` table for `EstimateFractions`
#'
#' @param mutrate_populations Character vector of the set of mutational populations
#' present in your data.
#' @return A `fraction_design` table that assumes that every possible combination of
#' mutational populations listed in `mutrate_populations` are present in your data.
#' The `present` column can be modified if this assumption is incorrect
#' @export
create_fraction_design <- function(mutrate_populations){

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
#' @param obj EZbakRDataobject
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param mutrate_populations Character vector of the set of mutational populations
#' that you want to infer the rates of mutations for.
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
#' significant performance gain without sacrificing accuracy.
#' @param strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item hierarchical: Estimate feature-specific new read mutation
#'  rate, regularizing the feature-specific estimate with a sample-wide prior. Currently
#'  only compatible with single mutation type mixture modeling.
#'  }
#' @param split_multi_features If a set of reads maps ambiguously to multiple features,
#' should data for such reads be copied for each feature in the ambiguous set? If this is
#' `TRUE`, then `multi_feature_cols` also must be set.
#' @param multi_feature_cols Character vector of columns that have the potential to
#' include assignment to multiple features. Only these columns will have their feaures split
#' if `split_multi_features` is `TRUE`.
#' @param multi_feature_sep String representing how ambiguous feature assignments are
#' distinguished in the feature names. For example, the default value of "+" denotes
#' that if a read maps to multiple features (call them featureA and featureB, for example),
#' then the feature column will have a value of "featureA+featureB".
#' @param pnew_prior_mean logit-Normal mean for logit(pnew) prior.
#' @param pnew_prior_sd logit-Normal sd for logit(pnew) prior.
#' @param pold_prior_mean logit-Normal mean for logit(pold) prior.
#' @param pold_prior_sd logit-Normal sd for logit(pold) prior.
#' @param hier_readcutoff If `strategy` == `hierarchical`, only features with this many reads
#' are used to infer the distribution of feature-specific labeled read mutation rates.
#' @param init_pnew_prior_sd If `strategy` == `hierarchical`, this is the initial logit(pnew)
#' prior standard deviation to regularize feature-specific labeled read mutation rate estimates.
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param overwrite If TRUE and an fractions estimate output already exists that
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
                              character_limit = 20,
                              overwrite = TRUE){


  UseMethod("EstimateFractions")


}



#' Estimate mutation rates
#'
#' @param obj EZbakRDataobject
#' @param populations Character vector of the set of mutational populations
#' that you want to infer the fractions of. For example, say your cB file contains
#' columns tracking T-to-C and G-to-A
#' @param strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model fit to all reads in a given sample.
#'  \item hierarchical (NOT YET IMPLEMENTED): Estimate feature-specific mutation
#'  rate with standard, regularizing the feature-specific
#'  estimate with a sample-wide prior.
#'  \item smalec (NOT YET IMPLEMENTED): Estimate two old read mutation rates, as was done in
#'  Smalec et al., 2023. Idea is that alignment artifacts can give rise to a
#'  high mutation rate old read population that should be accounted for
#'  to avoid overestimating the fraction of new reads
#' }
#' @param pnew_prior_mean logit-Normal mean for logit(pnew) prior.
#' @param pnew_prior_sd logit-Normal sd for logit(pnew) prior.
#' @param pold_prior_mean logit-Normal mean for logit(pold) prior.
#' @param pold_prior_sd logit-Normal sd for logit(pold) prior.
#' @import data.table
#' @export
EstimateMutRates <- function(obj,
                                        populations = "all",
                                        strategy = "standard",
                                        pnew_prior_mean = -2.94,
                                        pnew_prior_sd = 0.3,
                                        pold_prior_mean = -6.5,
                                        pold_prior_sd = 0.5
){

  UseMethod("EstimateMutRates")

}





###################################################
# EZBAKRDATA STANDARD METHODS
###################################################



#' Estimate fractions of each RNA population with standard EZbakRData object
#'
#' @import data.table
#' @importFrom magrittr %>%
#' @export
EstimateFractions.EZbakRData <- function(obj, features = "all",
                              mutrate_populations = "all",
                              fraction_design = NULL,
                              Poisson = TRUE,
                              strategy = c("standard", "hierarchical"),
                              split_multi_features = FALSE,
                              multi_feature_cols = NULL,
                              multi_feature_sep = "+",
                              pnew_prior_mean = -2.94,
                              pnew_prior_sd = 0.3,
                              pold_prior_mean = -6.5,
                              pold_prior_sd = 0.5,
                              hier_readcutoff = 300,
                              init_pnew_prior_sd = 0.8,
                              character_limit = 20,
                              overwrite = TRUE){

  `.` <- list

  ### Check that input is valid

  args <- c(as.list(environment()))

  if(!is(obj, "EZbakRData")){

    stop("obj is not an EZbakRData object!")

  }

  check_EstimateFractions_input(args)


  strategy <- match.arg(strategy)


  ### Estimate mutation rates
  message("Estimating mutation rates")
  obj <- EstimateMutRates(obj,
                           populations = mutrate_populations,
                           strategy = strategy,
                          pnew_prior_mean = pnew_prior_mean,
                          pnew_prior_sd = pnew_prior_sd,
                          pold_prior_mean = pold_prior_mean,
                          pold_prior_sd = pold_prior_sd)

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
  if(features == "all"){

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
  cB <- cB[,.(n = sum(n)), by = cols_to_group]

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

    cB <- split_features(cB, multi_feature_cols)

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
    cB <- cB[mutrates, nomatch = NULL]



    ### Estimate fraction new

    col_name <- paste0("logit_fraction_high", pops_to_analyze)
    uncertainty_col <- paste0("se_logit_fraction_high", pops_to_analyze)
    natural_col_name <- paste0("fraction_high", pops_to_analyze)

    if(strategy == "hierarchical"){


      message("FITTING HIERARCHICAL TWO-COMPONENT MIXTURE MODEL:")

      ### Get coverages to filter by; only want high coverage for feature-specific
      ### pnew estimation

      coverages <- cB[,.(coverage = sum(n)),
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
                                      pnew_prior_sd = init_pnew_prior_sd)),
          pold = unique(pold)),
        by = c("sample", features_to_analyze)
      ]


      feature_specific[, c("pnew", "lpnew_uncert") := .( inv_logit(sapply(params, `[[`, 2)),
                                         sapply(params, `[[`, 4))]


      ### Hierarchical model

      message("Estimating fractions with feature-specific pnews")

      global_est <- feature_specific[,.(pnew_prior = mean(logit(pnew)),
                                    pnew_prior_sd = sd(logit(pnew)) - mean(lpnew_uncert),
                                    pold = mean(pold)),
                                 by = sample]

      global_est[, pnew_prior_sd := ifelse(pnew_prior_sd < 0,
                                           pnew_prior_sd_min,
                                           pnew_prior_sd)
      ]

      global_est[, pnew_prior_sd := min(pnew_prior_sd)]

      setkey(global_est, sample)
      setkey(cB, sample)

      feature_specific <- cB[global_est, nomatch = NULL][,
                                                     .(params = list(fit_tcmm(muts = get(pops_to_analyze),
                                                                                   nucs = get(necessary_basecounts),
                                                                                   n = n,
                                                                                   Poisson = Poisson,
                                                                                   pold = unique(pold),
                                                                                   pnew_prior_mean = unique(pnew_prior),
                                                                                   pnew_prior_sd = unique(pnew_prior_sd))),
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
      obj$mutation_rates[[2]] <- feature_mutrates
      names(obj$mutration_rates) <- c(current_name, paste0("feature_", current_name))



    }else{

      fns <- dplyr::as_tibble(cB) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", features_to_analyze)))) %>%
        dplyr::summarise(fit = ifelse(!(unique(sample) %in% samples_with_no_label),
                                      list(I(optim(0,
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


    fns <- cB[,.(mixture_fit = list(fit_general_mixture(dataset = .SD,
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

  }


  # Save
  obj[['fractions']][[fraction_vect]] <- dplyr::as_tibble(fns)

  # Save metadata
  obj[['metadata']][['fractions']][[fraction_vect]] <- list(features = features_to_analyze,
                                                            populations = pops_to_analyze,
                                                            fraction_design = fraction_design)



  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)

}




#' Estimate mutation rates
#'
#' @param obj EZbakRDataobject
#' @param populations Character vector of the set of mutational populations
#' that you want to infer the fractions of. For example, say your cB file contains
#' columns tracking T-to-C and G-to-A
#' @param strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model fit to all reads in a given sample.
#'  \item hierarchical (NOT YET IMPLEMENTED): Estimate feature-specific mutation
#'  rate with standard, regularizing the feature-specific
#'  estimate with a sample-wide prior.
#'  \item smalec (NOT YET IMPLEMENTED): Estimate two old read mutation rates, as was done in
#'  Smalec et al., 2023. Idea is that alignment artifacts can give rise to a
#'  high mutation rate old read population that should be accounted for
#'  to avoid overestimating the fraction of new reads
#' }
#' @param pnew_prior_mean logit-Normal mean for logit(pnew) prior.
#' @param pnew_prior_sd logit-Normal sd for logit(pnew) prior.
#' @param pold_prior_mean logit-Normal mean for logit(pold) prior.
#' @param pold_prior_sd logit-Normal sd for logit(pold) prior.
#' @import data.table
#' @export
EstimateMutRates.EZbakRData <- function(obj,
                             populations = "all",
                             strategy = "standard",
                             pnew_prior_mean = -2.94,
                             pnew_prior_sd = 0.3,
                             pold_prior_mean = -6.5,
                             pold_prior_sd = 0.5
                             ){

  `.` <- list


  cB <- obj$cB
  cB <- data.table::setDT(cB)

  # Figure out which mutation counts are in the cB
  mutcounts_in_cB <- find_mutcounts(obj)

  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # Which populations to analyze?
  if(populations == "all"){

    muts_analyze <- mutcounts_in_cB

  }else{

    muts_analyze <- populations

  }

  nucs_analyze <- basecounts_in_cB[which(mutcounts_in_cB %in% muts_analyze)]


  # Infer proportion of each population
  mutest <- vector(mode = "list", length = length(muts_analyze))

  for(i in seq_along(muts_analyze)){

    group_cols <- c("sample", muts_analyze[i], nucs_analyze[i])

    mutest_dt <- cB[,.(n = sum(n)),
                    by = group_cols]

    mutest_temp <- mutest_dt[,.(params = list(fit_tcmm(muts = get(muts_analyze[i]),
                                                         nucs = get(nucs_analyze[i]),
                                                         n = n,
                                                         Poisson = Poisson,
                                                         pnew_prior_mean = pnew_prior_mean,
                                                         pnew_prior_sd = pnew_prior_sd,
                                                         pold_prior_mean = pold_prior_mean,
                                                         pold_prior_sd = pold_prior_sd))), by = sample]

    mutest_temp[, c("p1", "p2") := .(inv_logit(sapply(params, `[[`, 2)),
                                     inv_logit(sapply(params, `[[`, 3)))]

    mutest_temp[,c("pold",
                   "pnew") := .(min(c(p1, p2)),
                                        max(c(p1, p2))), by = 1:nrow(mutest_temp)]

    mutest_temp[,c("p1", "p2") := .(NULL, NULL)]

    mutest[[i]] <- mutest_temp

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



#' Estimate fractions of each RNA population with standard EZbakRData object
#'
#' @import data.table
#' @import arrow
#' @importFrom magrittr %>%
#' @export
EstimateFractions.EZbakRArrowData <- function(obj, features = "all",
                                              mutrate_populations = "all",
                                              fraction_design = NULL,
                                              Poisson = TRUE,
                                              strategy = c("standard", "hierarchical"),
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
                                              character_limit = 20,
                                              overwrite = TRUE){

  `.` <- list

  ### Check that input is valid

  args <- c(as.list(environment()))

  if(!is(obj, "EZbakRData")){

    stop("obj is not an EZbakRData object!")

  }

  check_EstimateFractions_input(args)


  strategy <- match.arg(strategy)


  ### Estimate mutation rates
  message("Estimating mutation rates")
  obj <- EstimateMutRates(obj,
                          populations = mutrate_populations,
                          strategy = strategy,
                          pnew_prior_mean = pnew_prior_mean,
                          pnew_prior_sd = pnew_prior_sd,
                          pold_prior_mean = pold_prior_mean,
                          pold_prior_sd = pold_prior_sd)

  mutation_rates <- obj$mutation_rates

  ### Figure out which features will be analyzed

  # cB columns
  cB <- obj$cBds
  cBschema <- arrow::schema(cBds)
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
  if(features == "all"){

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
    dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) > 0)) %>%
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



  for(s in seq_along(all_samples)){



    if(s %in% samples_without_label){

      ctl_cols_to_group <- c("sample", features_to_analyze)


      sample_fns <- cB %>%
        dplyr::filter(sample == all_samples[s]) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ctl_cols_to_group))) %>%
        dplyr::summarise(n = sum(n)) %>%
        dplyr::collect() %>%
        dplyr::filter(!dplyr::if_all(dplyr::all_of(features_to_analyze), ~ .x %in% c("NA", "__no_feature")))


      ### Split multi feature mappers if necessary
      if(split_multi_features){

        sample_fns <- split_features(sample_fns, multi_feature_cols)

      }


      dplyr::as_tibble(sample_fns) %>%
        dplyr::mutate(!!natural_col_name := 0,
                      !!col_name := -Inf,
                      !!uncertainty_col := 0)



    }else{

      message(paste0("Analyzing ", all_samples[s], "..."))

      if(Poisson){

        sample_cB <- cB %>%
          dplyr::filter(sample == all_samples[s]) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
          dplyr::summarise(reads = sum(n),
                           !!necessary_basecounts := sum(!!dplyr::sym(necessary_basecounts)*n)/sum(n)) %>%
          dplyr::collect() %>%
          dplyr::rename(n = reads) %>%
          dplyr::filter(!dplyr::if_all(dplyr::all_of(features_to_analyze), ~ .x %in% c("NA", "__no_feature")))

      }else{

        sample_cB <- cB %>%
          dplyr::filter(sample == all_samples[s]) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
          dplyr::summarise(n = sum(n)) %>%
          dplyr::collect() %>%
          dplyr::filter(!dplyr::if_all(dplyr::all_of(features_to_analyze), ~ .x %in% c("NA", "__no_feature")))

      }

      sample_cB <- setDT(sample_cB)


      ### Split multi feature mappers if necessary
      if(split_multi_features){

        sample_cB <- split_features(sample_cB, multi_feature_cols)


      }



      setkey(sample_cB, sample)


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
        sample_cB <- sample_cB[mutrates, nomatch = NULL]


        ### Estimate fraction new


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
                                            pnew_prior_sd = sd(logit(pnew)) - mean(lpnew_uncert),
                                            pold = mean(pold)),
                                         by = sample]

          global_est[, pnew_prior_sd := ifelse(pnew_prior_sd < 0,
                                               pnew_prior_sd_min,
                                               pnew_prior_sd)
          ]

          global_est[, pnew_prior_sd := min(pnew_prior_sd)]

          setkey(global_est, sample)
          setkey(sample_cB, sample)

          feature_specific <- sample_cB[global_est, nomatch = NULL][,
                                                                    .(params = list(fit_tcmm(muts = get(pops_to_analyze),
                                                                                             nucs = get(necessary_basecounts),
                                                                                             n = n,
                                                                                             Poisson = Poisson,
                                                                                             pold = unique(pold),
                                                                                             pnew_prior_mean = unique(pnew_prior),
                                                                                             pnew_prior_sd = unique(pnew_prior_sd))),
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
          obj$mutation_rates[[2]] <- feature_mutrates
          names(obj$mutration_rates) <- c(current_name, paste0("feature_", current_name))


        }else{

          ##### WHOLE CODE BLOCK NEARLY IDENTICAL TO AS IN DEFAULT METHOD
          sample_fns <- dplyr::as_tibble(sample_cB) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", features_to_analyze)))) %>%
            dplyr::summarise(fit = ifelse(!(unique(sample) %in% samples_with_no_label),
                                          list(I(optim(0,
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

        stop("Arrow backend is not currently compatible with > 1 mutation type modeling")

        ### Split multi feature mappers if necessary
        if(split_multi_features){

          sample_cB <- split_features(sample_cB, multi_feature_cols)


        }


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

  }


  # Save
  obj[['fractions']][[fraction_vect]] <- dplyr::as_tibble(fns)

  # Save metadata
  obj[['metadata']][['fractions']][[fraction_vect]] <- list(features = features_to_analyze,
                                                            populations = pops_to_analyze,
                                                            fraction_design = fraction_design)



  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)



}



#' Estimate mutation rates
#'
#' @param obj EZbakRDataobject
#' @param populations Character vector of the set of mutational populations
#' that you want to infer the fractions of. For example, say your cB file contains
#' columns tracking T-to-C and G-to-A
#' @param strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model fit to all reads in a given sample.
#'  \item hierarchical (NOT YET IMPLEMENTED): Estimate feature-specific mutation
#'  rate with standard, regularizing the feature-specific
#'  estimate with a sample-wide prior.
#'  \item smalec (NOT YET IMPLEMENTED): Estimate two old read mutation rates, as was done in
#'  Smalec et al., 2023. Idea is that alignment artifacts can give rise to a
#'  high mutation rate old read population that should be accounted for
#'  to avoid overestimating the fraction of new reads
#' }
#' @param pnew_prior_mean logit-Normal mean for logit(pnew) prior.
#' @param pnew_prior_sd logit-Normal sd for logit(pnew) prior.
#' @param pold_prior_mean logit-Normal mean for logit(pold) prior.
#' @param pold_prior_sd logit-Normal sd for logit(pold) prior.
#' @import data.table
#' @export
EstimateMutRates.EZbakRArrowData <- function(obj,
                                        populations = "all",
                                        strategy = "standard",
                                        pnew_prior_mean = -2.94,
                                        pnew_prior_sd = 0.3,
                                        pold_prior_mean = -6.5,
                                        pold_prior_sd = 0.5
){

  `.` <- list

  ### Figure out which features will be analyzed

  # cB columns
  cB <- obj$cBds
  cBschema <- arrow::schema(cBds)
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

  }

  nucs_analyze <- basecounts_in_cB[which(mutcounts_in_cB %in% muts_analyze)]


  # Infer proportion of each population
  mutest <- vector(mode = "list", length = length(muts_analyze))

  ### Which samples should be filtered out?

  metadf <- obj$metadf

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
    sample_mutest <- tibble()

    group_cols <- c("sample", muts_analyze[i], nucs_analyze[i])

    for(s in seq_along(samples_with_label)){

      cB <- cB %>%
        dplyr::filter(sample == samples_with_label[s]) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(n = sum(n)) %>%
        dplyr::collect()

      mutest_dt <- data.table::setDT(cB)

      mutest_temp <- mutest_dt[,.(params = list(fit_tcmm(muts = get(muts_analyze[i]),
                                                         nucs = get(nucs_analyze[i]),
                                                         n = n,
                                                         Poisson = Poisson,
                                                         pnew_prior_mean = pnew_prior_mean,
                                                         pnew_prior_sd = pnew_prior_sd,
                                                         pold_prior_mean = pold_prior_mean,
                                                         pold_prior_sd = pold_prior_sd))), by = sample]

      mutest_temp[, c("p1", "p2") := .(inv_logit(sapply(params, `[[`, 2)),
                                       inv_logit(sapply(params, `[[`, 3)))]

      mutest_temp[,c("pold",
                     "pnew") := .(min(c(p1, p2)),
                                  max(c(p1, p2))), by = 1:nrow(mutest_temp)]

      mutest_temp[,c("p1", "p2") := .(NULL, NULL)]

      sample_mutest <- bind_rows(sample_mutest, mutest_temp)

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
split_features <- function(cB, multi_feature_cols){

  ### Copy junction data more efficiently
  unique_juncs <- cB %>%
    dplyr::select(!!multi_feature_cols) %>%
    dplyr::distinct()


  ### have to move in and out of the tidyverse because this next line is much, much easier in data.table
  unique_juncs <- setDT(unique_juncs)
  new_temp_cols <- paste0("unrolled_", 1:length(multi_feature_cols))
  unique_juncs[,new_temp_cols := mget(multi_feature_cols)]


  ### Back to the tidyverse, because can't easily do this in data.table (sigh...)
  unique_juncs %>%
    tidyr::separate_rows(dplyr::all_of(multi_feature_cols), sep = "\\+")

  sample_cB <- sample_cB %>%
    inner_join(unique_juncs, by = multi_feature_cols,
               relationship = "many-to-many")


  ### Back to data.table again because its just easier and more efficient
  setDT(sample_cB)

  sample_cB[
    ,multi_feature_cols := mget(new_temp_cols)
  ][,new_temp_cols := NULL]


  return(sample_cB)

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
                  fraction_prior_sd = 1,
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

    ll <- fn*dpois(muts, nucs*pnew) + (1 - fn)*dpois(muts, nucs*pold)

  }else{

    ll <- fn*dbinom(muts, nucs, pnew) + (1 - fn)*dbinom(muts, nucs, pold)

  }


  ### Add prior weights if necessary

  if(estimate_pnew & estimate_pold){

    if(pold > pnew){

      prior <- dnorm(param[2], pnew_prior_mean, pnew_prior_sd,
                        log = TRUE) +
        dnorm(param[3], pold_prior_mean, pold_prior_sd,
                        log = TRUE) +
        dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
                        log = TRUE)

    }else{

      prior <- dnorm(param[2], pold_prior_mean, pold_prior_sd,
                        log = TRUE) +
        dnorm(param[3], pnew_prior_mean, pnew_prior_sd,
                        log = TRUE) +
        dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
              log = TRUE)


    }


  }else if(estimate_pnew){

    prior <- dnorm(param[2], pnew_prior_mean, pnew_prior_sd,
                        log = TRUE) +
      dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
            log = TRUE)


  }else if(estimate_pold){

    prior <- dnorm(param[2], pold_prior_mean, pold_prior_sd,
                    log = TRUE) +
      dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
            log = TRUE)

  }else{

    prior <- dnorm(param[1], fraction_prior_mean, fraction_prior_sd,
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
                     fraction_prior_sd = 1,
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

    uncertainty <- sqrt(diag(solve(fit$hessian)))

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

    likelihoods <- inv_logit(param[1])*dpois(muts,
                                             lambda = pnew*nucs ) +

      (1 - inv_logit(param[1]))*dpois(muts,
                                      lambda = pold*nucs )

  }else{

    likelihoods <- inv_logit(param[1])*dbinom(x = muts,
                                              size = nucs,
                                              prob = pnew) +

      (1 - inv_logit(param[1]))*dbinom(x = muts,
                                       size = nucs,
                                       prob = pold)


  }

  return(-sum(n*log(likelihoods)) - dnorm(param[1], log = TRUE))


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
        L <- L*dpois(dataset[[mutcols[j]]],
                      lambda = (as.numeric(mutrate_design[i, j])*pnew[j] +
                        pold[j]*(1 - as.numeric(mutrate_design[i, j])))*dataset[[basecols[j]]] )


      }else{

        # data.table syntax is an act of violence
        L <- L*dbinom(dataset[[mutcols[j]]],
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

    return(-sum(dataset[["n"]]*log(likelihoods)) - dnorm(param[1], log = TRUE))

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
                                pnew, pold, mutcols, basecols, highpop = NULL){


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
                        method = "L-BFGS-B")

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
                        method = "L-BFGS-B")

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


    outlist <- as.list(fit$par)


    names(outlist) <- listnames

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


