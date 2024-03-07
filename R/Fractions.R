#' Generate a `fraction_design` table for `EstimateFractions`
#'
#' @param mutrate_populations Character vector of the set of mutational populations
#' present in your data.
#' @return A `fraction_design` table that assumes that every possible combination of
#' mutational populations listed in `mutrate_populations` are present in your data.
#' The `present` column can be modified if this assumption is incorrect
create_fraction_design <- function(mutrate_populations){

  fraction_design <- dplyr::tibble(present = rep(TRUE,
                                                 times = 2^length(mutrate_populations)))

  fraction_list <- list()
  for(p in mutrate_populations){

    fraction_list[[p]] <- c("TRUE", "FALSE")

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
#' will use all feature columns.
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
#' This is because there is no RNA was exposed to both s^4U and s^6G, thus a population of reads
#' with both high T-to-C and G-to-A mutational content should not exist. To see an example
#' of what a TILAC `fraction_design` table could look like, see the lazily loaded
#' `tilac_fraction_design` object.
#'
#'
#' @param Poisson Use U-content adjusted Poisson mixture modeling strategy. Can see
#' significant performance gain without sacrificing accuracy.
#' @param strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item hierarchical (NOT YET IMPLEMENTED): Estimate feature-specific mutation
#'  rate with standard, regularizing the feature-specific
#'  estimate with a sample-wide prior.
#'  \item smalec (NOT YET IMPLEMENTED): Estimate two old read mutation rates, as was done in
#'  Smalec et al., 2023. Idea is that alignment artifacts can give rise to a
#'  high mutation rate old read population that should be accounted for
#'  to avoid overestimating the fraction of new reads
#' }
#' @import data.table
#' @importFrom magrittr %>%
EstimateFractions <- function(obj, features = "all",
                              mutrate_populations = "all",
                              fraction_design = NULL,
                              Poisson = TRUE,
                              strategy = "standard"){

  `.` <- list

  ### Check that input is valid

  args <- c(as.list(environment()))

  # check_EstimateFractions_input(args)


  ### Estimate mutation rates
  message("Estimating mutation rates")
  obj <- EstimateMutRates(obj,
                           populations = mutrate_populations,
                           strategy = strategy)

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

  ### Create fraction_design data frame if necessary

  if(is.null(fraction_design)){

    fraction_design <- create_fraction_design(pops_to_analyze)

  }


  ### Estimate fraction new for each feature in each sample
  message("Summarizing data for feature(s) of interest")

  cB <- setDT(cB)

  # Keep only feature of interest
  cols_to_group <- c(basecounts_in_cB, mutcounts_in_cB, "sample", features_to_analyze)
  cB[,.(n = sum(n)), by = cols_to_group]

  # Pair down design matrix
  fraction_design <- as_tibble(fraction_design)
  mutrate_design <- fraction_design %>%
    dplyr::filter(present) %>%
    dplyr::select(-present)

  # Summarize data
  if(Poisson){

    message("Averaging out the nucleotide counts for improved efficiency")
    cols_to_group_pois <- c("sample", mutcounts_in_cB, features_to_analyze)
    cols_to_avg <- setdiff(names(cB), c(cols_to_group_pois, "n"))

    cB <- cB[, c(lapply(.SD, function(x) sum(x*n)/sum(n) ),
                     .(n = sum(n))),
                 by = cols_to_group_pois, .SDcols = cols_to_avg]

  }

  # TO-DO: Add an estimated runtime here based on the number of rows in the cB,
  # the number of features, and what model is being run.
  message("Estimating fractions")
  if(ncol(mutrate_design) == 1){

    if(unlist(mutrate_design[1,1])){

      highpop <- 1

    }else{

      highpop <- 2

    }

    # dplyr rewrite that may somehow be much faster...
    # fns <- as_tibble(cB) %>%
    #   dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", features_to_analyze)))) %>%
    #   dplyr::do(dplyr::tibble(mixture_fit = list(dataset = .,
    #                                       mutrate_design = mutrate_design,
    #                                       mutcols = mutcounts_in_cB,
    #                                       basecols = basecounts_in_cB,
    #                                       Poisson = Poisson,
    #                                       highpop = highpop,
    #                                       twocomp = TRUE,
    #                                       pnew = sapply(mutcounts_in_cB,
    #                                                     function(name) mutation_rates[[name]]$pnew),
    #                                       pold = sapply(mutcounts_in_cB,
    #                                                     function(name) mutation_rates[[name]]$pold))))

    fns <- cB[,.(mixture_fit = list(fit_general_mixture(dataset = .SD,
                                                        mutrate_design = mutrate_design,
                                                        mutcols = mutcounts_in_cB,
                                                        basecols = basecounts_in_cB,
                                                        Poisson = Poisson,
                                                        highpop = highpop,
                                                        twocomp = TRUE,
                                                        pnew = sapply(mutcounts_in_cB,
                                                                      function(name) mutation_rates[[name]]$pnew),
                                                        pold = sapply(mutcounts_in_cB,
                                                                      function(name) mutation_rates[[name]]$pold)) )),
              by = c("sample", features_to_analyze)]


  }else{

    fns <- cB[,.(mixture_fit = list(fit_general_mixture(dataset = .SD,
                                                        mutrate_design = mutrate_design,
                                                        mutcols = mutcounts_in_cB,
                                                        basecols = basecounts_in_cB,
                                                        Poisson = Poisson,
                                                        twocomp = FALSE,
                                                        pnew = sapply(mutcounts_in_cB,
                                                                      function(name) mutation_rates[[name]]$pnew),
                                                        pold = sapply(mutcounts_in_cB,
                                                                      function(name) mutation_rates[[name]]$pold)) )),
              by = c("sample", features_to_analyze)]

  }




  message("Processing output")
  # Unroll the fraction estimates
  fns <- fns %>%
    tidyr::unnest_wider(mixture_fit)


  # What should output be named?
  fraction_vect <- paste(c("fractions", features_to_analyze), collapse = "_")

  obj[[fraction_vect]] <- fns

  # Add new class information
  if(!("EZbakRFractions" %in% class(obj))){

    class(obj) <- c("EZbakRFractions", class(obj))

  }

  return(obj)

}

# check_EstimateFractions_input(args){
#
#
#
# }



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
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item hierarchical (NOT YET IMPLEMENTED): Estimate feature-specific mutation
#'  rate with standard, regularizing the feature-specific
#'  estimate with a sample-wide prior.
#'  \item smalec (NOT YET IMPLEMENTED): Estimate two old read mutation rates, as was done in
#'  Smalec et al., 2023. Idea is that alignment artifacts can give rise to a
#'  high mutation rate old read population that should be accounted for
#'  to avoid overestimating the fraction of new reads
#' }
#' @import data.table
EstimateMutRates <- function(obj,
                             populations = "all",
                             strategy = "standard"
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

  # Helper function to run optim and organize output
  estimate_ps <- function(muts, nucs, n){

    fit <- stats::optim(par = c(0, -2, -4), fn = standard_binom_mix, muts = muts, nucs = nucs, n = n,
                 method = "L-BFGS-B", lower = c(-7, -7, -7), upper = c(7, 7, 7))
    return(list(p1 = fit$par[1], p2 = fit$par[2], p3 = fit$par[3]))

  }


  # Infer proportion of each population
  mutest <- vector(mode = "list", length = length(muts_analyze))

  for(i in seq_along(muts_analyze)){

    group_cols <- c("sample", muts_analyze[i], nucs_analyze[i])

    mutest_dt <- cB[,.(n = sum(n)),
                    by = group_cols]

    mutest_temp <- mutest_dt[,.(params = list(estimate_ps(muts = get(muts_analyze[i]),
                                                         nucs = get(nucs_analyze[i]),
                                                         n = n))), by = sample]

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


# Logit and sigmoid functions that I use a lot
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))


# Simple binomial mixture likelihood
standard_binom_mix <- function(param, muts, nucs, n){

  logl <- sum(n*log( inv_logit(param[1])*dbinom(muts, nucs, inv_logit(param[2]))
                     + (1 - inv_logit(param[1]))*dbinom(muts, nucs, inv_logit(param[3])) ) )

  return(-logl)

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

# Abstract the concept of an NR-seq mixture model to oblivion
generalized_likelihood <- function(param, dataset, Poisson = TRUE,
                                   mutrate_design, pnew, pold, highpop,
                                   mutcols, basecols, twocomp = TRUE){


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
fit_general_mixture <- function(dataset, Poisson = TRUE, mutrate_design, twocomp,
                                pnew, pold, mutcols, basecols, highpop = NULL){


  if(twocomp){

    fit <- stats::optim(par = 0,
                        fn = generalized_likelihood,
                        dataset = dataset,
                        Poisson = Poisson,
                        mutrate_design = mutrate_design,
                        pnew = pnew,
                        pold = pold,
                        mutcols = mutcols,
                        basecols = basecols,
                        highpop = highpop,
                        twocomp = twocomp,
                        upper = 7,
                        lower = -7,
                        method = "L-BFGS-B")

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

  }




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
  listnames <- lapply(population_list, function(x) paste0(muttypes, x))
  listnames <- sapply(listnames, function(x) paste(x, collapse = "_"),
                      simplify = "vector")


  if(twocomp){

    outlist <- list(fit$par, logit(1 - inv_logit(fit$par)))


  }else{

    outlist <- as.list(fit$par)

  }

  names(outlist) <- listnames

  return(outlist)


}



