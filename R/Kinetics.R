#' Generic function for estimating kinetic parameters.
#'
#' @param obj An `EZbakRFractions` object, which is an `EZbakRData` object on
#' which `EstimateFractions()` has been run.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param strategy Kinetic parameter estimation strategy.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item pulse-chase (NOT YET IMPLEMENTED): Estimate kdeg for a pulse-chase experiment. A kdeg will be estimated
#'  for each time point at which label was present. This includes any pulse-only samples,
#'  as well as all samples including a chase after the pulse.
#'  \item custom (NOT YET IMPLEMENTED): Provide a custom function that takes
#'  fraction estimates as input and produces as output kinetic parameter estimates.
#' }
#' @param quant_name Name of quantification tool appended to table name of interest.
#' This is only relevant if you are providing isoform-specific estimates, in which
#' case the isoform quantification tool's name may need to be provided in order
#' for EZbakR to uniquely identify the table of interest. Even in that case though,
#' this should only have to be non-null in the case where you have performed isoform-specific
#' fraction estimation with more than one quantification tool's output.
#' @return `EZbakRKinetics` object.
#' @import data.table
#' @export
EstimateKinetics <- function(obj,
                             strategy = c("standard", "tilac"),
                             features = NULL,
                             populations = NULL,
                             fraction_design = NULL,
                             character_limit = 20,
                             overwrite = TRUE){

  ### Check that input is valid

  # EZbakRData object on which EstimateFractions() has been run?
  if(!is(obj, "EZbakRFractions")){

    if(is(obj, "EZbakRData")){

      stop("obj is not an EZbakRFractions object! Run `obj <- EstimateFractions(obj, ...)`,
           where ... represents optional parameters, before running `EstimateKinetics`.")

    }else{

      stop("obj is not an EZbakRData object with fraction estimates!")

    }


  }

  # Analysis strategy
  strategy <- match.arg(strategy)

  # features
  if(!is.null(features)){

    if(!is.character(features)){

      stop("features is not a character vector!")

    }

  }




  ### "Method dispatch"

  if(strategy == "standard"){

    obj <- Standard_kinetic_estimation(obj,
                                features = features,
                                populations = populations,
                                fraction_design = fraction_design,
                                character_limit = character_limit,
                                overwrite = overwrite)


  }else if(strategy == "tilac"){

    obj <- tilac_ratio_estimation(obj,
                                  features = features,
                                  populations = populations,
                                  fraction_design = fraction_design,
                                  character_limit = character_limit,
                                  overwrite = overwrite)

  }


  return(obj)

}


# kdeg = -log(1 - fn)/tl
# ksyn = (normalized read count)*kdeg
Standard_kinetic_estimation <- function(obj, features = NULL,
                                        populations = NULL,
                                        fraction_design = NULL,
                                        character_limit = 20,
                                        overwrite = TRUE){



  `.` <- list


  ### Figure out which fraction new estimates to use

  # Function is in Helpers.R
  fractions_name <- EZget(obj,
                          features = features,
                          populations = populations,
                          fraction_design = fraction_design,
                          returnNameOnly = TRUE)


  # Get fractions
  kinetics <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[["metadata"]][["fractions"]][[fractions_name]][["features"]]



  ### Estimate kdegs

  # Determine which column to use for kinetic parameter estimation
  fraction_cols <- colnames(kinetics)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  setDT(data.table::copy(kinetics))

  # Add label time info
  metadf <- obj$metadf

  metadf <- setDT(data.table::copy(metadf))
  setkey(metadf, sample)
  setkey(kinetics, sample)

  kinetics <- kinetics[metadf[,c("sample", "tl")], nomatch = NULL]

  # Make sure no -s4U controls made it through
  kinetics <- kinetics[tl > 0]


  kinetics[, kdeg := -log(1 - get(fraction_of_interest))/tl]
  kinetics[, log_kdeg := log(kdeg)]


  ### Normalize read counts


  # TO-DO: ALLOW USERS TO JUST USE THE TPM FROM ISOFORM QUANTIFICATION
  reads_norm <- get_normalized_read_counts(obj = obj,
                                           features_to_analyze = features_to_analyze,
                                           fractions_name = fractions_name)



  ### Estimate uncertainty in log(kdeg)

  lkdeg_uncert <- function(fn, se_lfn){

    deriv <- (1/log(1 - fn)) * (1 / (1 - fn)) * (fn*(1 - fn))

    uncert <- abs(deriv)*se_lfn

    return(uncert)

  }

  se_of_interest <- paste0("se_logit_", fraction_of_interest)

  kinetics[, se_log_kdeg := lkdeg_uncert(fn = get(fraction_of_interest),
                                         se_lfn = get(se_of_interest))]


  ### Estimate ksyn

  # Merge with kinetics
  setkeyv(reads_norm, c("sample", features_to_analyze))
  setkeyv(kinetics, c("sample", features_to_analyze))

  cols_to_keep <- c("sample", features_to_analyze, "normalized_reads", "scale_factor")
  kinetics_cols_to_keep <- c(cols_to_keep, "kdeg", "log_kdeg", "se_log_kdeg", "n")

  kinetics <- kinetics[reads_norm[,..cols_to_keep], ..kinetics_cols_to_keep, nomatch = NULL]

  # Estimate ksyn
  kinetics[, ksyn := normalized_reads*kdeg]
  kinetics[, log_ksyn := log(ksyn)]

  # Estimate uncertainty (assuming normalized_reads ~ Poisson(normalized_reads)/scale_factor)
  kinetics[, se_log_ksyn := sqrt( (1/(normalized_reads*scale_factor)) + se_log_kdeg^2)]



  ##### PROCESS OUTPUT AND RETURN

  # What should output be named?
  kinetics_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")

  if(nchar(kinetics_vect) > character_limit){

    num_kinetics <- length(obj[['kinetics']])
    kinetics_vect <- paste0("kinetics", num_kinetics + 1)

  }


  # Are there any metadata objects
  if(length(obj[['metadata']][['kinetics']]) > 0){

    kinetics_vect <- decide_output(obj,
                                   proposed_name = kinetics_vect,
                                   type = "kinetics",
                                   features = features_to_analyze,
                                   kstrat = "standard",
                                   overwrite = overwrite)

  }


  readcount_vect <- kinetics_vect

  if(length(obj[['metadata']][['readcounts']]) > 0){

    readcount_vect <- decide_output(obj,
                                   proposed_name = readcount_vect,
                                   type = "readcounts",
                                   features = features_to_analyze,
                                   counttype = "TMM_normalized",
                                   overwrite = overwrite)

  }

  obj[["kinetics"]][[kinetics_vect]] <- dplyr::as_tibble(kinetics) %>%
    dplyr::select(sample, !!features_to_analyze, kdeg, log_kdeg, se_log_kdeg, ksyn, log_ksyn, se_log_ksyn, normalized_reads, n)

  # Eventually want to add count matrix output
  obj[["readcounts"]][[readcount_vect]] <- dplyr::as_tibble(reads_norm) %>%
    dplyr::select(sample, !!features_to_analyze, n, normalized_reads, geom_mean, scale_factor)


  obj[["metadata"]][["kinetics"]][[kinetics_vect]] <- list(features = features_to_analyze,
                                                           kstrat = "standard")
  obj[["metadata"]][["readcounts"]][[readcount_vect]] <- list(features = features_to_analyze,
                                                             counttype = "TMM_normalized")


  if(!is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)


}


# Estimate TILAC ratio
tilac_ratio_estimation <- function(obj,
                                   features = features,
                                   populations = populations,
                                   fraction_design = fraction_design,
                                   character_limit = character_limit){


  `.` <- list

  ### Figure out which fraction new estimates to use

  # Need to determine which columns of the cB to group reads by
  if(is.null(features)){

    fractions_name <- names(obj)[grepl("fractions_", names(obj))]

    if(length(fractions_name) > 1){

      stop("There is more than one fractions estimate data frame; therefore,
           you need to explicit specify a `features` vector to let EZbakR
           know which of these you would like to use!")

    }

    features_to_analyze <- unname(unlist(strsplit(fractions_name, "_")))
    lf <- length(features_to_analyze)
    features_to_analyze <- features_to_analyze[2:lf]

  }else{

    supposed_fractions_name <- paste0("fractions_", paste(features, collapse = "_"))

    if(!(supposed_fractions_name %in% names(obj))){

      stop("features do not have an associated fractions data frame!")

    }else{

      features_to_analyze <- features

    }

  }

  # Name of fractions table to use
  fractions_table_name <- paste(c("fractions", features_to_analyze), collapse = "_")

  # Get fractions
  kinetics <- obj[[fractions_table_name]]



  ### Estimate TILAC ratio

  # Determine which column to use for kinetic parameter estimation
  fraction_cols <- colnames(kinetics)
  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  # TODO: COULD ADD PARAMETER TO ALLOW USER TO DECIDE WHICH FRACTION COLUMNS TO USE
  if(length(fraction_of_interest) != 2){

    stop("There are not exactly two `high` mutational content fractions!
         Are you sure this is TILAC data?")

  }

  setDT(data.table::copy(kinetics))

  # Add label time info
  metadf <- obj$metadf
  metacols <- colnames(metadf)

  metadf <- setDT(data.table::copy(metadf))
  setkey(metadf, sample)
  setkey(kinetics, sample)

  tl_cols <- metacols[grepl("tl", metacols)]


  cols_keep <- c("sample", tl_cols)
  kinetics <- kinetics[metadf[,..cols_keep], nomatch = NULL]

  # Make sure no -s4U controls made it through
    # Not sure how to easily do this filtering in data.table
  kinetics <- setDT(dplyr::as_tibble(kinetics) %>%
    dplyr::filter(dplyr::if_all(dplyr::starts_with("tl"), ~ . > 0)))


  kinetics[, TILAC_ratio := get(fraction_of_interest[1])/get(fraction_of_interest[2])]
  kinetics[, log_TILAC_ratio := TILAC_ratio]



  ### Get normalized read counts

  reads_norm <- get_normalized_read_counts(obj,
                                           features_to_analyze = features_to_analyze)


  ### Estimate ksyn

  # Merge with kinetics
  setkeyv(reads_norm, c("sample", features_to_analyze))
  setkeyv(kinetics, c("sample", features_to_analyze))

  cols_to_keep <- c("sample", features_to_analyze, "reads", "normalized_reads")

  kinetics <- kinetics[reads_norm[,..cols_to_keep], nomatch = NULL]


  # Figure out what to name output
  kinetics_name <- paste(features_to_analyze, collapse = "_")
  reads_name <- paste(c("readcounts", features_to_analyze), collapse = "_")

  obj[['kinetics']][[kinetics_name]] <- dplyr::as_tibble(kinetics)

  # Eventually want to add count matrix output
  obj[['readcounts']][[reads_name]] <- list(reads_df = reads_norm)

  if(!is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)



}




