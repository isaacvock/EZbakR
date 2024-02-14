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
