#' Simple simulation function
#'
#' @param nreads Number of reads to simulate
#' @param fn Fraction of reasd that are new in simulation. Whether a read will
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
