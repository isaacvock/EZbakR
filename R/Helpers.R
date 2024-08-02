

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
#' @param isoforms If the relevant table is the result of isoform deconvolution
#' @param kstrat Only relevant if `type` == "kinetics". Short for "kinetics strategy";
#' the strategy used to infer kinetic parameters.
#' @param parameter Only relevant if `type` == "averages" or "comparisons". Which
#' parameter was being averaged or compared?
#' @param counttype String denoting what type of read count information you are looking
#' for. Current options are "TMM_normalized", "transcript", and "matrix". TO-DO:
#' Not sure this is being used in any way currently...
#' @param condition Condition specified in relevant run of `CompareParameters()`.
#' Therefore, only relevant if type == "comparisons".
#' @param experimental Experimental condition specified in relevant run of `CompareParameters()`.
#' Therefore, only relevant if type == "comparisons".
#' @param reference Reference condition specified in relevant run of `CompareParameters()`.
#' Therefore, only relevant if type == "comparisons".
#' @param mean_vars Sample features from metadf that were used in formula for
#' parameter average linear model. Only relevant if type == "averages".
#' @param sd_vars Sample features from metadf that were used in formula for
#' parameter standard deviation linear model. Only relevant if type == "averages".
#' @param sub_features Only relevant if type == "dynamics". Feature columns
#' that distinguished between the different measured species when running
#' `EZDynamics()`.
#' @param grouping_features Only relevant if type == "dynamics. Features
#' that were the overarching feature assignments by which `sub_features` were grouped
#' when running `EZDynamics()`.
#' @param sample_feature Only relevant if type == "dynamics". Name of the metadf
#' column that distinguished the different classes of samples when running
#' `EZDynamics()`.
#' @param graph Only relevant if type == "dynamics". NxN adjacency matrix,
#' where N represents the number of species modeled when running `EZDynamics()`.
#' @param repeatID Numerical ID for duplicate objects with same metadata.
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
                  type = c("fractions", "kinetics", "readcounts",
                           "averages", "comparisons", "dynamics"),
                  features = NULL,
                  populations = NULL,
                  fraction_design = NULL,
                  isoforms = NULL,
                  kstrat = NULL,
                  parameter = NULL,
                  counttype = NULL,
                  condition = NULL,
                  experimental = NULL,
                  reference = NULL,
                  mean_vars = NULL,
                  sd_vars = NULL,
                  repeatID = NULL,
                  sub_features = NULL,
                  grouping_features = NULL,
                  sample_feature = NULL,
                  graph = NULL,
                  returnNameOnly = FALSE,
                  exactMatch = FALSE,
                  alwaysCheck = FALSE){



  type <- match.arg(type)

  if(!is.null(counttype)){

    counttype <- match.arg(counttype, c("TMM_normalized", "transcript", "matrix"))

  }

  if(!is.null(kstrat)){

    kstrat <- match.arg(kstrat, c("standard", "tilac", "NSS",
                                  "shortfeed", "pulse-chase"))

  }

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

  if(!is.null(isoforms)){

    possible_tables_iso <- exact_ezsearch(metadata,
                                          query = isoforms,
                                          object = "isoforms")


    possible_tables <- intersect(possible_tables, possible_tables_iso)


  }


  if(!is.null(kstrat)){

    possible_tables_ks <- exact_ezsearch(metadata,
                                          query = kstrat,
                                          object = "kstrat")


    possible_tables <- intersect(possible_tables, possible_tables_ks)


  }


  if(!is.null(parameter)){

    possible_tables_par <- exact_ezsearch(metadata,
                                          query = parameter,
                                          object = "parameter")


    possible_tables <- intersect(possible_tables, possible_tables_par)


  }


  if(!is.null(counttype)){

    possible_tables_cnt <- exact_ezsearch(metadata,
                                          query = counttype,
                                          object = "counttype")


    possible_tables <- intersect(possible_tables, possible_tables_cnt)


  }

  if(!is.null(condition)){

    possible_tables_c <- exact_ezsearch(metadata,
                                          query = condition,
                                          object = "condition")


    possible_tables <- intersect(possible_tables, possible_tables_c)


  }

  if(!is.null(experimental)){

    possible_tables_e <- exact_ezsearch(metadata,
                                          query = experimental,
                                          object = "experimental")


    possible_tables <- intersect(possible_tables, possible_tables_e)


  }

  if(!is.null(reference)){

    possible_tables_r <- exact_ezsearch(metadata,
                                          query = reference,
                                          object = "reference")


    possible_tables <- intersect(possible_tables, possible_tables_r)


  }

  if(!is.null(repeatID)){

    possible_tables_rep <- exact_ezsearch(metadata,
                                          query = repeatID,
                                          object = "repeatID")

    possible_tables <- intersect(possible_tables, possible_tables_rep)

  }


  if(!is.null(mean_vars)){

    possible_tables_mv <- vector_ezsearch(metadata,
                                          mean_vars,
                                          "mean_vars",
                                          exactMatch = exactMatch)

    possible_tables <- intersect(possible_tables, possible_tables_mv)


  }


  if(!is.null(sd_vars)){

    possible_tables_sd <- vector_ezsearch(metadata,
                                          sd_vars,
                                          "sd_vars",
                                          exactMatch = exactMatch)

    possible_tables <- intersect(possible_tables, possible_tables_sd)


  }

  if(!is.null(sub_features)){

    possible_tables_sf <- vector_ezsearch(metadata,
                                          sub_features,
                                          "sub_features",
                                          exactMatch = exactMatch)

    possible_tables <- intersect(possible_tables, possible_tables_sf)

  }

  if(!is.null(grouping_features)){

    possible_tables_gf <- vector_ezsearch(metadata,
                                          grouping_features,
                                          "grouping_features",
                                          exactMatch = exactMatch)

    possible_tables <- intersect(possible_tables, possible_tables_gf)

  }

  if(!is.null(sample_feature)){

    possible_tables_sf <- exact_ezsearch(metadata,
                                         sample_feature,
                                         "sample_feature")

    possible_tables <- intersect(possible_tables, possible_tables_sf)

  }

  if(!is.null(graph)){

    possible_tables_g <- exact_ezsearch(metadata,
                                         graph,
                                         "graph")

    possible_tables <- intersect(possible_tables, possible_tables_g)


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

    if(length(possible_tables) == 0){

      stop("No tables fit your search criterion!")

    }

    return(obj[[type]][[possible_tables]])

  }




}

vector_ezsearch <- function(metadata,
                            queries,
                            object,
                            exactMatch){


  potential_tables <- c()

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


exact_ezsearch <- function(metadata,
                           query,
                           object){

  lm <- length(metadata)
  possible_tables <- c()


  for(m in 1:lm){

    subject <- metadata[[m]][[object]]

    if(identical(subject, query)){

      possible_tables <- c(possible_tables, names(metadata)[m])

    }

  }

  return(possible_tables)


}



###########################
### DETERMINE OUTPUT NAME
###########################

decide_output <- function(obj, proposed_name,
                          type = c("fractions", "kinetics",
                                   "readcounts", "averages",
                                   "comparisons", "dynamics"),
                          features = NULL,
                          populations = NULL,
                          fraction_design = NULL,
                          counttype = NULL,
                          kstrat = NULL,
                          parameter = NULL,
                          condition = NULL,
                          reference = NULL,
                          experimental = NULL,
                          mean_vars = NULL,
                          sd_vars = NULL,
                          sub_features = NULL,
                          grouping_features = NULL,
                          sample_feature = NULL,
                          graph = NULL,
                          overwrite = TRUE){

  type = match.arg(type)

  ### Does same analysis output already exist?
  existing_output <- EZget(obj,
                             type = type,
                             features = features,
                             populations = populations,
                             fraction_design = fraction_design,
                             parameter = parameter,
                             kstrat = kstrat,
                             counttype = counttype,
                           condition = condition,
                           reference = reference,
                           experimental = experimental,
                           mean_vars = mean_vars,
                           sd_vars = sd_vars,
                           sub_features = sub_features,
                           grouping_features = grouping_features,
                           sample_feature = sample_feature,
                           graph = graph,
                             returnNameOnly = TRUE,
                             exactMatch = TRUE,
                           alwaysCheck = TRUE)

  if(is.null(existing_output)){

    # Have to find a name that doesn't exist
    if(proposed_name %in% names(obj[[type]])){

      need_new_name <- TRUE
      num_repeats <- 2

      while(need_new_name){

        test_proposed_name <- paste0(proposed_name, num_repeats)

        if(test_proposed_name %in% names(obj[[type]])){

          num_repeats <- num_repeats + 1

        }else{

          proposed_name <- test_proposed_name
          need_new_name <- FALSE

        }

      }


    }


  }else if(overwrite){

    proposed_name <- existing_output


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

  return(proposed_name)

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


#' Get normalized read counts from either a cB table or `EZbakRFractions` object.
#'
#' Uses TMM normalization strategy, similar to that used by DESeq2 and edgeR.
#'
#' @param obj An `EZbakRData` or `EZbakRFractions` object.
#' @param features_to_analyze Features in relevant table
#' @param fractions_name Name of fractions table to use
#' @export
get_normalized_read_counts <- function(obj,
                                       features_to_analyze,
                                       fractions_name = NULL){

  UseMethod("get_normalized_read_counts")

}

#' Get normalized read counts from either an `EZbakRFractions` object.
#'
#' Uses TMM normalization strategy, similar to that used by DESeq2 and edgeR.
#'
#' @param obj An `EZbakRFractions` object.
#' @param features_to_analyze Features in relevant table
#' @param fractions_name Name of fractions table to use
#' @export
get_normalized_read_counts.EZbakRFractions <- function(obj,
                                                       features_to_analyze,
                                                       fractions_name = NULL){

  reads <- data.table::setDT(data.table::copy(obj[['fractions']][[fractions_name]]))

  reads <- normalize_reads(reads, features_to_analyze)

  return(reads)


}

#' Get normalized read counts from either a cB table in an `EZbakRData` object.
#'
#' Uses TMM normalization strategy, similar to that used by DESeq2 and edgeR.
#'
#' @param obj An `EZbakRData` object.
#' @param features_to_analyze Features in relevant table
#' @param fractions_name Name of fractions table to use
#' @export
get_normalized_read_counts.default <- function(obj,
                                               features_to_analyze,
                                               fractions_name = NULL){

  ### Hack to deal with devtools::check() NOTEs
  n <- NULL

  `.` <- list


  ### Get normalized read counts

  cB <- data.table::setDT(data.table::copy(obj$cB))

  # Calc read counts for each feature
  reads <- cB[,.(n = sum(n)), by = c("sample", features_to_analyze)]


  # Normalize
  reads <- normalize_reads(reads, features_to_analyze)

  return(reads)

}

normalize_reads <- function(reads, features_to_analyze){

  ### Hack to deal with devtools::check() NOTEs
  geom_mean <- n <- normalized_reads <- scale_factor <- NULL

  `.` <- list

  # Median of ratios normalization
  reads[, geom_mean := exp(mean(log(n))), by = features_to_analyze]
  scales <- reads[, .(scale_factor =  stats::median(n/geom_mean)), by = .(sample)]

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

