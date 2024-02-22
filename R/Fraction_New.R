#' Estimate fractions of each RNA population
#'
#' @param obj EZbakRDataobject
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns.
#' @param mutrate_populations Character vector of the set of mutational populations
#' that you want to infer the rates of mutations for.
#' @param fraction_design "Design matrix" specifying which RNA populations exist
#' in your samples
#' @param pnew_strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item hierarchical: Estimate feature-specific mutation rate with standard , regularizing the feature-specific
#'  estimate with a sample-wide prior.
#'  \item smalec: Estimate two old read mutation rates, as was done in
#'  Smalec et al., 2023. Idea is that alignment artifacts can give rise to a
#'  high mutation rate old read population that should be accounted for
#'  to avoid overestimating the fraction of new reads
#' }
#' @param pold_strategy String denoting which old read mutation rate estimation strategy
#' @import data.table
#' @importFrom magrittr %>%
EstimateFractions <- function(obj, features = "all",
                              mutrate_populations = "all",
                              fraction_design = NULL,
                              pnew_strategy = "standard",
                              pold_strategy = "standard"){

  `.` <- list

  ### Check that input is valid

  args <- c(as.list(environment()))

  # check_EstimateFractions_input(args)


  ### Estimate mutation rates
  obj <- EstimateMutRates(obj,
                                     populations = mutrate_populations,
                                     pnew_strategy = pnew_strategy,
                                     pold_strategy = pold_strategy)

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


  ### Estimate fraction new for each feature in each sample
  cB <- setDT(cB)

  # Pair down design matrix
  fraction_design <- as_tibble(fraction_design)
  mutrate_design <- fraction_design %>%
    dplyr::filter(present) %>%
    dplyr::select(-present)

  # Estimate fractions
  fns <- cB[,.(mixture_fit = list(fit_general_mixture(dataset = .SD,
                                                      mutrate_design = mutrate_design,
                                                      mutcols = mutcounts_in_cB,
                                                      basecols = basecounts_in_cB,
                                                      pnew = sapply(mutcounts_in_cB,
                                                                    function(name) mutation_rates[[name]]$pnew),
                                                      pold = sapply(mutcounts_in_cB,
                                                                    function(name) mutation_rates[[name]]$pold)) )),
            by = c("sample", features_to_analyze)]


  # Unroll the fraction estimates
  fns <- cbind(fns %>% dplyr::select(-mixture_fit),
               as.data.frame(fns$mixture_fit))


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
#' @param pnew_strategy String denoting which new read mutation rate estimation strategy to use.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item intronic: Estimate sample-wide pnew from intronic mutation rate. Have to
#'  specify which feature in the EZbakRDataobject cB represents the intronic reads
#'  \item hierarchical: Estimate feature-specific pnew with standard , regularizing the feature-specific
#'  estimate with a sample-wide prior.
#' }
#' @import data.table
#' @param pold_strategy String denoting which old read mutation rate estimation strategy
EstimateMutRates <- function(obj,
                             populations = "all",
                             pnew_strategy = "standard",
                             pold_strategy = "standard"
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

    fit <- optim(par = c(0, -2, -4), fn = mixture_lik, muts = muts, nucs = nucs, n = n,
                 method = "L-BFGS-B", lower = c(-7, -7, -7), upper = c(7, 7, 7))
    return(list(p1 = fit$par[1], p2 = fit$par[2], p3 = fit$par[3]))

  }


  # Infer proportion of each population
  mutest <- vector(mode = "list", length = length(muts_analyze))

  for(i in seq_along(muts_analyze)){

    mutest_temp <- cB[,.(params = list(estimate_ps(muts = get(muts_analyze[i]),
                                                         nucs = get(nucs_analyze[i]),
                                                         n = n)))]

    mutest_temp[, c("p1", "p2") := .(inv_logit(sapply(params, `[[`, 2)),
                                     inv_logit(sapply(params, `[[`, 3)))]

    mutest_temp[,c("pold",
                   "pnew") := .(min(c(p1, p2)),
                                        max(c(p1, p2)))]

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
generalized_likelihood <- function(param, dataset,
                                   mutrate_design, pnew, pold,
                                   mutcols, basecols){


  likelihoods <- rep(0, times = length(dataset[[1]]))

  # Fractions of populations
  fractions <- softmax(param)


  # Get full likelihood for generalized mutational models
  for(i in 1:nrow(mutrate_design)){

    L <- 1

    for(j in 1:ncol(mutrate_design)){

      # data.table syntax is an act of violence
      L <- L*dbinom(dataset[[mutcols[j]]],
                    dataset[[basecols[j]]],
                    prob = as.numeric(mutrate_design[i, j])*pnew[j] +
                      pold[j]*(1 - as.numeric(mutrate_design[i, j])))



    }

    likelihoods <- likelihoods + L*fractions[i]

  }

  return(-sum(dataset[["n"]]*log(likelihoods)))

}

# Fit and process the output of the generalized likelihood model
fit_general_mixture <- function(dataset, mutrate_design,
                                pnew, pold, mutcols, basecols){


  fit <- optim(par = rep(0, times = nrow(mutrate_design)),
                         fn = generalized_likelihood,
                         dataset = dataset,
                         mutrate_design = mutrate_design,
                         pnew = pnew,
                         pold = pold,
                         mutcols = mutcols,
                         basecols = basecols)

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



  outlist <- as.list(fit$par)
  names(outlist) <- listnames

  return(outlist)


}

