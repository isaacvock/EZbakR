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
                           features, quant_name = NULL){


  fnames <- names(obj[[tabletype]])



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


    ### Look for features in the names of the fractions data frames

    # Feature names will show up with '_'s removed
    features <- gsub("_","",features)


    # Don't search for fullfit output
    if(tabletype = 'averages'){

      fnames <- fnames[!grepl('fullfit', fnames)]

    }

    features_in_fnames <- strsplit(fnames, split = "_")




    # Determine whether or not to search for quantification tool name
    if(is.null(quant_name)){

      # Remove quantification tool name if present
      features_in_fnames <- features_in_fnames[!(features_in_fnames %in%
                                                 c("custom", "salmon", "sailfish",
                                                   "alevin", "piscem", "kallisto",
                                                   "rsem", "stringtie"))]

      present <- unlist(lapply(features_in_fnames, function(x) {
        return(all(features %in% x ) & all(x %in% features))
      })
      )


    }else{

      present <- unlist(lapply(features_in_fnames, function(x) {
        return(all(c(features, quant_name) %in% x ) & all(x %in% c(features, quant_name)))
      })
      )


    }




    if(length(present) == 0){


      # Remove quantification tool name if present
      features_in_fnames <- features_in_fnames[!(features_in_fnames %in%
                                                   c("custom", "salmon", "sailfish",
                                                     "alevin", "piscem", "kallisto",
                                                     "rsem", "stringtie"))]

      # Print out the allowed features
      features_allowed <- lapply(features_in_fnames, function(x){
              return(paste0('c(',paste(paste0("'",x,"'"), collapse = ', '), ')'))
            })

      message(paste(c("`features` should be one of the following (unordered) vectors: ",
                      unlist(features_allowed)),
                    collapse = '\n'))

      stop(paste0("features do not match any of the ", tabletype," tables present!"))

    }

    if(length(present) > 1){

      if(is.null(quant_name)){

        stop(paste0("Features match more than one ", tabletype, " table! If you
                    are trying to specify an isoform-specific analysis table,
                    and you have multiple such tables, make sure to specify
                    `quant_name`"))


      }

      stop(paste0("Features match more than one ", tabletype, " table!
                  Alert the developer by posting an Issues on the EZbakR,
                  Github as this edge case should not ever happen."))


    }


  }

  table_name <- fnames[present]

  if("isoforms" %in% features){
    isoform_specific <- TRUE
  }else{
    isoform_specific <- FALSE
  }


  return(list(table_name = table_name,
              isoform_specific = isoform_specific))

}


