#' Generic function for estimating kinetic parameters
#'
#' @export
EstimateKinetics <- function(obj,
                             features = "all",
                             strategy = "standard"){

  ### Check that input is valid

  strategy <- match.arg(strategy)


  ### "Method dispatch"

  if(strategy == "standard"){

    Standard_kinetic_estimation(obj,
                                features = features)


  }


}


# Kinetic parameters = -log(1 - fn)/tl
Standard_kinetic_estimation <- function(obj, features = "all"){

  ### Figure out which fraction new estimates to use

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


  # Name of fractions table to use
  fractions_table_name <- paste(c("fractions", features_to_analyze), collapse = "_")

  # Get fractions
  kinetics <- obj[[fractions_table_name]]



  ### Estimate kdegs

  # Determine which column to use for kinetic parameter estimation
  fraction_cols <- colnames(kinetics)

  fraction_of_interest <- fraction_cols[grepl("high", fraction_cols)]

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

  cB <- obj$cB

  # Calc read counts for each feature
  reads <- cB[,.(reads = sum(n)), by = c("sample", features_to_analyze)]

  # Median of ratios normalization
  reads[, geom_mean := exp(mean(log(reads))), by = features_to_analyze]
  scales <- reads[, .(scale_factor =  median(reads/geom_mean)), by = .(sample)]

  setkey(scales, sample)
  setkey(reads, sample)

  # Normalize read counts
  reads_norm <- reads[scales, nomatch = NULL]

  reads_norm[,normalized_reads := reads/scale_factor]


  ### Estimate ksyn

  # Merge with kinetics
  setkeyv(reads_norm, c("sample", features_to_analyze))
  setkeyv(kinetics, c("sample", features_to_analyze))

  cols_to_keep <- c("sample", features_to_analyze, "normalized_reads")

  kinetics <- kinetics[reads_norm[,..cols_to_keep], nomatch = NULL]

  # Estimate ksyn
  kinetics[, ksyn := normalized_reads*kdeg]


  # Figure out what to name output
  kinetics_name <- paste(c("kinetics", features_to_analyze), collapse = "_")
  reads_name <- paste(c("readcounts", features_to_analyze), collapse = "_")

  obj[[kinetics_name]] <- kinetics

  # Eventually want to add count matrix output
  obj[[reads_name]] <- list(reads_df = reads_norm)

  if(!is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)


}

