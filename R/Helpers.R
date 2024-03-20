############################
### READ COUNT NORMALIZATION
############################

get_normalized_read_counts <- function(obj,
                                       features_to_analyze){

  UseMethod("get_normalized_read_counts")

}

get_normalized_read_counts.EZbakRFractions <- function(obj,
                                                       features_to_analyze){

  fraction_name <- paste0("fractions_", paste(features_to_analyze, collapse = "_"))

  reads <- data.table::setDT(obj[[fraction_name]])

  reads <- normalize_reads(reads, features_to_analyze)

  return(reads)

}

get_normalized_read_counts.default <- function(obj,
                                               features_to_analyze){

  ### Get normalized read counts

  cB <- obj$cB

  # Calc read counts for each feature
  reads <- cB[,.(n = sum(n)), by = c("sample", features_to_analyze)]


  # Normalize
  reads <- normalize_reads(reads, features_to_analyze)

  return(reads)

}

normalize_reads <- function(reads, features_to_analyze){

  # Median of ratios normalization
  reads[, geom_mean := exp(mean(log(n))), by = features_to_analyze]
  scales <- reads[, .(scale_factor =  median(n/geom_mean)), by = .(sample)]

  setkey(scales, sample)
  setkey(reads, sample)

  # Normalize read counts
  reads_norm <- reads[scales, nomatch = NULL]

  reads_norm[,normalized_reads := n/scale_factor]

  return(reads_norm)

}



##################################
# DETERMINE FEATURES BEING TRACKED
##################################


# Infer features from EZbakR objects
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

