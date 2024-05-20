
############################
### MISCELLANEOUS
############################



# Logit and sigmoid functions that I use a lot
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))




############################
### READ COUNT NORMALIZATION
############################



#' @export
get_normalized_read_counts <- function(obj,
                                       features_to_analyze,
                                       isoform_fraction_name = NULL){

  UseMethod("get_normalized_read_counts")

}

#' @export
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

#' @export
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

  }else if(objtype == 'averages'){

    averages <- obj

    ### What are the features?

    avg_cols <- colnames(averages)

    substrings <- avg_cols[grepl("^mean_", avg_cols) |
                             grepl("^sd_", avg_cols) |
                             grepl("^coverage_", avg_cols)]

    features <- avg_cols[!(avg_cols %in% substrings)]

    if(!(length(features) > 0)){

      rlang::abort(
        "No features are defined in your averages data frame!",
        class = "no_averages_features"
      )

    }


    return(features)


  }else{

    stop("objtype is undefined!")

  }



}


# What is the name of the table of fractions/kinetics/etc. to analyze?
get_table_name <- function(obj, tabletype,
                           features, quant_name = NULL,
                           parameter = NULL){



  fnames <- names(obj[[tabletype]])

  # Don't search for fullfit output
  if(tabletype == 'averages'){

    fnames <- fnames[!grepl('fullfit', fnames)]

  }


  if(is.null(features)){

    ### Use the only available fractions; or throw an error if there
    ### is more than one fractions and thus auto-detection is impossible


    if(length(fnames) > 1){

      stop(paste0("There is more than one ", tabletype," data frame; therefore,
           you need to explicit specify a `features` vector to let EZbakR
           know which of these you would like to use!"))

    }

    features_to_analyze <- unname(unlist(strsplit(fnames, "_")))

    table_name <- paste(features_to_analyze, collapse = "_")


  }else{


    ### Look for features in the names of the fractions data frames

    # Feature names will show up with '_'s removed
    features <- gsub("_","",features)



    features_in_fnames <- strsplit(fnames, split = "_")


    if(!is.null(quant_name)){

      features <- c(features, quant_name)

    }else{

      # Remove quantification tool name if present
      features_in_fnames <- lapply(features_in_fnames,
                                   function(x){
                                     return(x[!(x %in% c("custom", "salmon", "sailfish",
                                                         "alevin", "piscem", "kallisto",
                                                         "rsem", "stringtie"))])
                                   })

    }

    if(!is.null(parameter)){

      features <- c(features, gsub("_", "",parameter))

    }else{

      # Remove quantification tool name if present
      features_in_fnames <- lapply(features_in_fnames,
                                   function(x){
                                     return(x[!(x %in% gsub("_", "",parameter))])
                                   })

    }






    # Look for table of interest
    present <- unlist(lapply(features_in_fnames, function(x) {
      return(all(features %in% x ) & all(x %in% features))
    })
    )




    if(sum(present) == 0){


      if(!is.null(parameter)){

        # Remove quantification tool name if present
        features_in_fnames <- lapply(features_in_fnames,
                                     function(x){
                                       return(x[!(x %in% c("custom", "salmon", "sailfish",
                                                           "alevin", "piscem", "kallisto",
                                                           "rsem", "stringtie", gsub("_","",parameter)))])
                                     })

      }else{

        # Remove quantification tool name if present
        features_in_fnames <- lapply(features_in_fnames,
                                     function(x){
                                       return(x[!(x %in% c("custom", "salmon", "sailfish",
                                                           "alevin", "piscem", "kallisto",
                                                           "rsem", "stringtie"))])
                                     })

      }


      # Print out the allowed features
      features_allowed <- lapply(features_in_fnames, function(x){
              return(paste0('c(',paste(paste0("'",x,"'"), collapse = ', '), ')'))
            })

      message(paste(c("`features` should be one of the following (unordered) vectors: ",
                      unlist(features_allowed)),
                    collapse = '\n'))

      stop(paste0("features do not match any of the ", tabletype," tables present!"))

    }

    if(sum(present) > 1){


      if(is.null(quant_name)){

        stop(paste0("Features match more than one ", tabletype, " table! If you
                    are trying to specify an isoform-specific analysis table,
                    and you have multiple such tables, make sure to specify
                    `quant_name`"))


      }

      stop(paste0("Features match more than one ", tabletype, " table!
                  Alert the developer by posting an Issues on the EZbakR,
                  Github as this edge case should never happen."))


    }


    table_name <- fnames[present]


  }


  if(grepl("^isoforms", table_name)){
    isoform_specific <- TRUE
  }else{
    isoform_specific <- FALSE
  }


  return(list(table_name = table_name,
              isoform_specific = isoform_specific))

}



# What comparative analysis table do you want?
get_comparison <- function(obj,
               condition,
               reference,
               experimental,
               features){

  # Function to get rid of the comparison data frame
  remove_comparison <- function(input_list) {

    input_list <- input_list[names(input_list) != "comparison"]

    return(input_list)
  }


  metadata <- lapply(obj, remove_comparison)

  if(length(metadata) == 1){

    # Index of list of interest is the only list of interest
    index <- 1

  }else{

    # Function to find particular metadata
    check_meta <- function(input_list, element, value){

      if(input_list[[element]] == value){

        return(TRUE)

      }else{

        return(FALSE)

      }

    }


    if(!is.null(condition)){

      condition_check <- check_meta(metadata, "condition", condition)

    }else{

      condition_check <- rep(TRUE, times = length(metadata))

    }

    if(!is.null(reference)){

      reference_check <- check_meta(metadata, "reference", reference)

    }else{

      reference_check <- rep(TRUE, times = length(metadata))

    }

    if(!is.null(experimental)){

      experimental_check <- check_meta(metadata, "experimental", experimental)

    }else{

      experimental_check <- rep(TRUE, times = length(metadata))

    }

    if(!is.null(feature)){

      features_check <- check_meta(metadata, "features", features)

    }else{

      features_check <- rep(TRUE, times = length(metadata))

    }


    index <- condition_check &
      reference_check &
      experimental_check &
      features_check



  }



  return(obj[[index]])


}
