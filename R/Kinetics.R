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
#' @return `EZbakRKinetics` object.
#' @import data.table
#' @export
EstimateKinetics <- function(obj,
                             features = NULL,
                             strategy = c("standard", "tilac"),
                             quant_name = NULL){

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
                                quant_name = quant_name)


  }else if(strategy == "tilac"){

    obj <- tilac_ratio_estimation(obj,
                           features = features)

  }


  return(obj)

}


# kdeg = -log(1 - fn)/tl
# ksyn = (normalized read count)*kdeg
Standard_kinetic_estimation <- function(obj, features = NULL,
                                        quant_name = NULL){


  `.` <- list


  ### Figure out which fraction new estimates to use

  fnames <- names(obj[['fractions']])

  if(is.null(features)){

    ### Use the only available fractions; or throw an error if there
    ### is more than one fractions and thus auto-detection is impossible

    fractions_name <- names(fnames)

    if(length(fractions_name) > 1){

      stop("There is more than one fractions estimate data frame; therefore,
           you need to explicit specify a `features` vector to let EZbakR
           know which of these you would like to use!")

    }

    features_to_analyze <- unname(unlist(strsplit(fractions_name, "_")))

    fractions_name <- paste(features_to_analyze, collapse = "_")


  }else{

    ### Look for features in the names of the fractions data frames,

    supposed_fractions_name <- paste(gsub("_","",features), collapse = "_")

    if(!(supposed_fractions_name %in% names(fnames))){

      ### If that didn't work, check to see if one or more of the fractions
      ### names contain the expected feature vector as part of an
      ### isoform-specific fractions object

      fractions_name <- fnames[grepl(supposed_fractions_name, fnames) & grepl("isoforms_", fnames)]

      isoform_specific <- TRUE

      if(length(fractions_name) > 1){

        ### If there is more than one feature that matches the bill,
        ### narrow down by grepping for the provided quant_name

        if(is.null(quant_name)){

          stop("You appear to be requesting isoform level kinetic analyses, but
               have multiple isoform-level fraction estimates. Specify `quant_name`
               to ")

        }

        fractions_name <- fractions_name[grepl(paste0("_", quant_name), fractions_name)]

      }else if(length(fractions_name) == 0){

        stop("features do not have an associated fractions data frame!")


      }


    }else{

      features_to_analyze <- features

      fractions_name <- paste(features_to_analyze, collapse = "_")

    }

  }

  # Get fractions
  kinetics <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- get_features(kinetics, objtype = "fractions")



  ### Estimate kdegs

  # Determine which column to use for kinetic parameter estimation
  fraction_cols <- colnames(kinetics)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  setDT(kinetics)

  # Add label time info
  metadf <- obj$metadf

  setkey(metadf, sample)
  setkey(kinetics, sample)

  kinetics <- kinetics[metadf[,c("sample", "tl")], nomatch = NULL]

  # Make sure no -s4U controls made it through
  kinetics <- kinetics[tl > 0]


  kinetics[, kdeg := -log(1 - inv_logit(get(fraction_of_interest)))/tl]
  kinetics[, log_kdeg := log(kdeg)]


  ### Get normalized read counts

  if(isoform_specific){

    reads_norm <- get_normalized_read_counts(obj = obj,
                                             features_to_analyze = features_to_analyze,
                                             isoform_fraction_name = fractions_name)


  }else{

    reads_norm <- get_normalized_read_counts(obj = obj,
                                             features_to_analyze = features_to_analyze)


  }


  ### Estimate ksyn

  # Merge with kinetics
  setkeyv(reads_norm, c("sample", features_to_analyze))
  setkeyv(kinetics, c("sample", features_to_analyze))

  cols_to_keep <- c("sample", features_to_analyze, "normalized_reads")

  kinetics <- kinetics[reads_norm[,..cols_to_keep], nomatch = NULL]

  # Estimate ksyn
  kinetics[, ksyn := normalized_reads*kdeg]


  # Figure out what to name output
  reads_name <- paste0("normalized_readcounts_", fractions_name)

  obj[["kinetics"]][[fractions_name]] <- kinetics %>%
    dplyr::select(-!!fraction_of_interest, dplyr::everything(), n, normalized_reads, tl)

  # Eventually want to add count matrix output
  obj[["readcounts"]][[reads_name]] <- reads_norm %>%
    dplyr::select(sample, !!features_to_analyze, n, normalized_reads, geom_mean, scale_factor)

  if(!is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)


}


# Estimate TILAC ratio
tilac_ratio_estimation <- function(obj,
                                   features){


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

  setDT(kinetics)

  # Add label time info
  metadf <- obj$metadf
  metacols <- colnames(metadf)

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

  obj[['kinetics']][[kinetics_name]] <- kinetics

  # Eventually want to add count matrix output
  obj[['readcounts']][[reads_name]] <- list(reads_df = reads_norm)

  if(!is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)



}




