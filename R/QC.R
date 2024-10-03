EZQC <- function(obj,
                 mutrate_populations = "all",
                 ...){

  UseMethod("EZQC")

}


EZQC.EZbakRFractions <- function(obj,
                                 mutrate_populations = "all",
                                 ...){


  ### Check raw mutation rates
  graw <- check_raw_mutation_rates(obj)


  ### Check labeled and unlabeled mutation rates
  glabel <- check_plabeled(obj)


  ### Check read count replicate correlation
  gread <- check_read_count_corr(obj)


  ### Check fraction labeled distribution
  gfn <- check_fl_dist(obj)


  ### Check fraction labeled replicate correlation
  gfn <- check_fl_corr(obj)



}


EZQC.EZbakRData <- function(obj,
                            mutrate_populations = "all",
                            ...,
                            features = "all"){

  ### Check raw mutation rates
  graw <- check_raw_mutation_rates(obj)


  ### Check labeled and unlabeled mutation rates
  mutrates <- EstimateMutRates(obj,
                               populations = mutrate_populations)


  glabel <- check_plabeled(obj, mutrate_populations
                           mutrates)


  ### Check read count replicate correlation
  gread <- check_read_count_corr(obj)


}


check_raw_mutation_rates <- function(){

}

check_plabeled <- function(obj, mutrate_populations,
                           mutrates = NULL){

  # Gotta fetch mutrates
  if(is.null(mutrates)){

    mutrates <- obj$mutation_rates

  }

  # Plot distribution of pold and pnew for each sample
  glist <- vector(mode = "list",
                  length = length(mutrates))

}

check_read_count_corr <- function(){


}


check_fl_dist <- function(){

}

check_fl_corr <- function(){


}
