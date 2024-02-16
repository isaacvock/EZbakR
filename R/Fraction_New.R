#' Estimate fractions of each RNA population
#'
#' @param obj EZbakRDataobject
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns.
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
EstimateFractions <- function(obj, features = "all",
                              populations = "all",
                              pnew_strategy = "standard",
                              pold_strategy = "standard"){

  mutation_rates <- EstimateMutRates(obj,
                                     pnew_strategy = pnew_strategy,
                                     pold_strategy = pold_strategy)





}



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

    mutest_temp[,c(paste0("pold_", muts_analyze[i]),
                   paste0("pnew_", muts_analyze[i])) := .(min(c(p1, p2)),
                                        max(c(p1, p2)))]

    mutest_temp[,c("p1", "p2") := .(NULL, NULL)]

    mutest[[i]] <- mutest_temp

  }


  return(mutest)

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
