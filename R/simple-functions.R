# Generic for inferring features from EZbakR objects
get_features <- function(obj, objtype = "cB"){

  if(objtype == "cB"){

    stop("Not defined yet!")

  }else if(objtype == "fractions"){

    fractions <- obj

    ### What are the features?

    fraction_cols <- colnames(fractions)

    substrings <- fraction_cols[grepl("fraction_", fraction_cols)]

    features <- fraction_cols[!(fraction_cols %in% c(substrings, "n", "sample"))]

    if(!(length(features) > 0)){

      rlang::abort(
        "No features are defined in your fractions data frame!",
        class = "no_fractions_features"
      )

    }


    return(features)

  }else{

    stop("objtype is undefined!")

  }



}

#' Simple simulation function
#'
#' @param nreads Number of reads to simulate
#' @param fn Fraction of reads that are new in simulation. Whether a read will
#' be new will be determined by a draw from a Bernoulli(fn) distribution.
#' @param pnew T-to-C mutation rate in new reads
#' @param pold T-to-C mutation rate in old reads
#' @param rlen Length of simulated reads
#' @param Ucont Fraction of nucleotides in simulated reads that are Ts (U in RNA)
#' @importFrom magrittr %>%
#' @export
SimpleSim <- function(nreads = 1000, fn = 0.5,
                      pnew = 0.05, pold = 0.001,
                      rlen = 100, Ucont = 0.25){

  ### CHECK INPUT

  args <- c(as.list(environment()))


  check_SimpleSim_input(args)


  ### SIMULATE DATA

  # The dumb hack that gets check() to shut up
  TC <- nT <- NULL


  # Determine newness status of each read
    # 1 = new; 0 = old
  newness <- stats::rbinom(nreads, size = 1, prob = fn)

  # Simulate a number of Ts for each read
  nTs <- stats::rbinom(nreads, size = rlen, prob = Ucont)

  # Simulate mutational content in each read
    # Cute trick is used here to simulate a different mutation rate for
    # new reads (pnew) and old reads (pold).
  nTC <- stats::rbinom(nreads, size = nTs, prob = newness*pnew + (1-newness)*pold)

  # Aggregate data into a data frame
  sim_df <- dplyr::tibble(nT = nTs,
                   TC = nTC) %>%
    dplyr::group_by(nT, TC) %>%
    dplyr::count()

  return(sim_df)


}


check_SimpleSim_input <- function(args){


  ### nreads
  if(!is.numeric(args$nreads)){
    stop("nreads must be numeric!")
  }

  if(args$nreads <= 0){
    stop("nreads must be > 0!")
  }


  ### fn
  if(!is.numeric(args$fn)){
    stop("fn must be numeric!")
  }

  if(args$fn < 0){
    stop("fn must be >= 0!")
  }

  if(args$fn > 1){
    stop("fn must be <= 1!")
  }


  ### pnew
  if(!is.numeric(args$pnew)){
    stop("pnew must be numeric!")
  }

  if(args$pnew < 0){
    stop("pnew must be >= 0!")
  }

  if(args$pnew > 1){
    stop("pnew must be <= 1!")
  }


  ### pold
  if(!is.numeric(args$pnew)){
    stop("pnew must be numeric!")
  }

  if(args$pnew < 0){
    stop("pold must be >= 0!")
  }

  if(args$pold > 1){
    stop("pold must be <= 1!")
  }


  ### pold vs. pnew
  if(args$pold >= args$pnew){
    stop("pold must be strictly less than pnew!")
  }


}


SimpleMultiPopSim <- function(nreads = 1000, populations = c("TC"),
                              high_mutrates = c(0.05), low_mutrates = c(0.002),
                              fractions = NULL,
                              Nuc_cont = 0.25, rlen = 100, ngenes = 1){


  ### label vectors with population identities
  names(high_mutrates) <- populations
  names(low_mutrates) <- populations
  num_pops <- length(populations)


  ### Number of reads for each gene
  nreads_per_gene <- as.vector(rmultinom(1,
                               size = nreads,
                               prob = rep(1, times = ngenes)))


  ### Simulate or extract proportions for each population
  # ordered in such a way that elements 1:ngenes are fractions for population 1,
  # elements (ngenes+1):2*ngenes are factions for population 2, etc.
  if(is.null(fractions)){

    fractions <- rep(1/(2^num_pops + 1),
                     times = 2^num_pops)


  }


  ### Simulate data for each population

  # Truth tables for presence of each population
  mutlist <- vector(mode = "list", length = length(populations))

  for(p in seq_along(populations)){

    mutlist[[p]] <- 0:1

  }

  names(mutlist) <- populations

  # Create mutational design matrix
  mutdesign <- expand.grid(mutlist)


  # Nucleotide content
  nucleotides <- unique(substr(populations, start = 1, stop = 1))
  nucleotide_cnts <- paste0("n", nucleotides)

  # Mutation
  datalist <- vector(mode = "list", length = (num_pops*length(nucleotides) + 2))

  names(datalist) <- c(populations, nucleotide_cnts, "sample", "gene")

  datalist[["sample"]] <- "Simulation"
  datalist[["gene"]] <- rep(1:ngenes, times = nreads_per_gene)

  browser()
  # Populate nucleotide counts
  for(n in seq_along(nucleotide_cnts))

    datalist[[nucleotide_cnts[n]]] <- round(rlen*Nuc_cont)

  for(p in seq_along(populations)){



    if(length(fractions) > num_pops){

      index_low <- (p-1)*ngenes + 1
      index_high <- p*ngenes

      newness <- rbinom(nreads,
                        size = 1,
                        rep(fractions[index_low:index_high],
                            times = nreads_per_gene))

    }else{

      identity <- rmultinom(1, size = nreads,
                            probs = fractions)

      newnewss

    }

    datalist[[populations[p]]] <- rbinom(nreads,
                            size = round(rlen*Nuc_cont),
                            prob = newness*(high_mutrates[p]) +
                              (1 - newness)*(low_mutrates[p]))



  }


  data_df <- dplyr::as_tibble(as.data.frame(datalist))

  data_df <- data_df %>%
    dplyr::group_by_all() %>%
    dplyr::count()

  return(data_df)





}
