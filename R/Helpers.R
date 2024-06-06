

############################
### DATA GETTER
############################

#' Easily get EZbakR table of estimates of interest
#'
#' @param obj EZbakRData object
#' @param type The class of EZbakR outputs would you like to search through.
#' Equivalent to the name of the list in the EZbakRData object that contains
#' the tables of interest.
#' @param features Features that must be present in the table of interest.
#' If `exactMatch` is TRUE, then these features must also be the only features
#' present in the table.
#' @param populations Only relevant if `type` == "fractions". Mutational
#' populations that must have been analyzed to generate the table of interest.
#' @param fraction_design Only relevant if `type` == "fractions". Fraction design
#' table used to generate the table of interest.
#' @param parameter Only relevant if `type` == "averages" or "comparisons". Which
#' parameter was being averaged or compared?
#' @param returnNameOnly If TRUE, then only the names of tables that passed your
#' search criteria will be returned. Else, the single table passing your search
#' criteria will be returned. If there is more than one table that passes your
#' search criteria and `returnNameOnly` == `FALSE`, an error will be thrown.
#' @param exactMatch If TRUE, then for `features` and `populations`, which can be vectors,
#' ensure that provided vectors of features and populations exactly match the relevant metadata
#' vectors.
#' @param alwaysCheck If TRUE, then even if there is only a single table for the `type`
#' of interest, still run all checks against queries.
#' @export
EZget <- function(obj,
                  type = c("fractions", "kinetics",
                           "averages", "comparisons"),
                  features = NULL,
                  populations = NULL,
                  fraction_design = NULL,
                  parameter = NULL,
                  returnNameOnly = FALSE,
                  exactMatch = FALSE,
                  alwaysCheck = FALSE){


  type = match.arg(type)

  metadata <- obj[['metadata']][[type]]


  table_of_interest <- NULL


  # If only one table is present, that's the one you want
  if(length(metadata) == 1 & !alwaysCheck){

    table_of_interest <- names(metadata)

    if(returnNameOnly){

      return(table_of_interest)

    }else{

      return(obj[[type]][[table_of_interest]])

    }

    break

  }

  ### DESIGN IDEA
  # For each searchable metadata element, there are different
  # search strategies:
  # 1) features and populations: Either check for exact vector
  # equality (if exactMatch == TRUE) or check that all elements
  # in provided vector are in relevant metadata vector.
  # 2) fraction_design: make sure data frame is identical through anti_join
  # 3) parameters: exact string match

  possible_tables <- names(metadata)

  if(!is.null(features)){

    possible_tables_f <- vector_ezsearch(metadata,
                                         features,
                                         "features",
                                         exactMatch = exactMatch)

    possible_tables <- intersect(possible_tables, possible_tables_f)


  }


  if(!is.null(populations)){

    possible_tables_p <- vector_ezsearch(metadata,
                                         populations,
                                         "populations",
                                         exactMatch = exactMatch)

    possible_tables <- intersect(possible_tables, possible_tables_p)



  }


  if(!is.null(fraction_design)){

    lm <- length(metadata)
    possible_tables_fd <- c()

    for(m in 1:lm){

      fd_subject <- metadata[[m]][['fraction_design']]


      cnames <- colnames(fd_subject)

      if(all(cnames %in% colnames(fraction_design))){

        aj_test <- dplyr::anti_join(fd_subject,
                                    fraction_design,
                                    by = cnames)


      }else{

        aj_test <- data.frame(dummy = 1)

      }

      if(nrow(aj_test) == 0){

        possible_tables_fd <- c(possible_tables_fd, names(metadata)[m])

      }



    }


    possible_tables <- intersect(possible_tables, possible_tables_fd)



  }


  if(!is.null(parameter)){

    lm <- length(metadata)
    possible_tables_par <- c()


    for(m in 1:lm){

      par_subject <- metadata[[m]][['parameter']]

      if(par_subject == parameter){

        possible_tables_par <- c(possible_tables_par, names(metadata)[m])

      }

    }


    possible_tables <- intersect(possible_tables, possible_tables_par)


  }


  # Can return multiple names of candidate tables, but can
  # only return one table
  if(returnNameOnly){

    if(length(possible_tables) > 1){

      warning("More than one table fits your search criterion!")

    }

    return(possible_tables)

  }else{

    if(length(possible_tables) > 1){

      stop("More than one table fits your search criterion!")

    }

    return(obj[[type]][[possible_tables]])

  }



}



vector_ezsearch <- function(metadata,
                                    queries,
                                    object = c("features", "populations"),
                            exactMatch){


  potential_tables <- c()

  object <- match.arg(object)
  lmeta <- length(metadata)

  for(m in 1:lmeta){

    subjects <- metadata[[m]][[object]]

    if(exactMatch){

      if(all(queries %in% subjects) &
         all(subjects %in% queries)){

        potential_tables <- c(potential_tables, names(metadata)[m])

      }


    }else{

      if(all(queries %in% subjects)){

        potential_tables <- c(potential_tables, names(metadata)[m])

      }


    }

  }

  return(potential_tables)


}



###########################
### DETERMINE OUTPUT NAME
###########################

decide_output <- function(obj, proposed_name,
                          type = c("fractions"),
                          features = NULL,
                          populations = NULL,
                          fraction_design = NULL,
                          parameter = NULL,
                          overwrite = TRUE){

  ### Does same analysis output already exist?
  existing_output <- EZget(obj,
                             type = "fractions",
                             features = features,
                             populations = populations,
                             fraction_design = fraction_design,
                             parameter = parameter,
                             returnNameOnly = TRUE,
                             exactMatch = TRUE)

  if(is.null(existing_output)){

    # Have to find a name that doesn't exist
    if(proposed_name %in% names(obj[[type]])){

      need_new_name <- TRUE
      num_repeats <- 2

      while(need_new_name){

        test_proposed_name <- paste0(proposed_name, "_", num_repeats)

        if(test_proposed_name %in% names(obj[[type]])){

          num_repeats <- num_repeats + 1

        }else{

          proposed_name <- test_proposed_name
          need_new_name <- FALSE

        }

      }


    }


  }else if(overwrite){

    return(existing_output)


  }else{

    # Create overwrite name
    need_new_name <- TRUE
    num_repeats <- 2

    while(need_new_name){

      test_proposed_name <- paste0(existing_output, "_", num_repeats)

      if(test_proposed_name %in% names(obj[[type]])){

        num_repeats <- num_repeats + 1

      }else{

        proposed_name <- test_proposed_name
        need_new_name <- FALSE

      }

    }

  }

}





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
