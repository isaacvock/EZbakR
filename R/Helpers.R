############################
### READ COUNT NORMALIZATION
############################

get_normalized_read_counts <- function(obj,
                                       features_to_analyze,
                                       isoform_fraction_name = NULL){

  UseMethod("get_normalized_read_counts")

}

get_normalized_read_counts.EZbakRFractions <- function(obj,
                                                       features_to_analyze,
                                                       isoform_fraction_name = NULL){

  if(!is.null(isoform_fraction_name)){

    reads <- data.table::copy(data.table::setDT(obj[['fractions']][[isoform_fraction_name]]))

    reads <- normalize_reads(reads, features_to_analyze)

    return(reads)

  }else{

    fraction_name <-  paste(gsub("_","",features_to_analyze), collapse = "_")

    reads <- data.table::copy(data.table::setDT(obj[['fractions']][[fraction_name]]))

    reads <- normalize_reads(reads, features_to_analyze)

    return(reads)

  }



}

get_normalized_read_counts.default <- function(obj,
                                               features_to_analyze,
                                               isoform_fraction_name = NULL){

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

    features <- fraction_cols[!(fraction_cols %in% c(substrings, "n", "sample",
                                                     "expected_count", "effective_length",
                                                     "RPK"))]

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


# What is the name of the table of fractions/kinetics/etc. to analyze?
get_table_name <- function(obj, tabletype,
                           features){

  fnames <- names(obj[[tabletype]])
  isoform_specific <- FALSE

  if(is.null(features)){

    ### Use the only available fractions; or throw an error if there
    ### is more than one fractions and thus auto-detection is impossible

    table_name <- names(fnames)

    if(length(table_name) > 1){

      stop(paste0("There is more than one ", tabletype," data frame; therefore,
           you need to explicit specify a `features` vector to let EZbakR
           know which of these you would like to use!"))

    }

    features_to_analyze <- unname(unlist(strsplit(table_name, "_")))

    table_name <- paste(features_to_analyze, collapse = "_")


  }else{

    ### Look for features in the names of the fractions data frames,

    supposed_table_name <- paste(gsub("_","",features), collapse = "_")

    if(!(supposed_table_name %in% fnames)){

      isoform_specific <- TRUE

      ### If that didn't work, check to see if one or more of the fractions
      ### names contain the expected feature vector as part of an
      ### isoform-specific fractions object

      table_name <- fnames[grepl(supposed_table_name, fnames)]


      if(length(table_name) > 1){

        ### If there is more than one feature that matches the bill,
        ### narrow down by grepping for the provided quant_name

        if(is.null(quant_name)){

          stop(paste0("You appear to be requesting isoform level kinetic analyses, but
               have multiple isoform-level ", tabletype," tables! Specify `quant_name`
               to tell EZbakR which of these to use."))

        }

        table_name <- table_name[grepl(paste0("_", quant_name), table_name)]

      }else if(length(table_name) == 0){

        stop(paste0("features do not have an associated ", tabletype," data frame!"))


      }


    }else{

      table_name <- paste(features, collapse = "_")

    }

  }

  return(list(table_name = table_name,
              isoform_specific = isoform_specific))

}


