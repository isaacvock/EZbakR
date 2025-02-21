#' Generic function for estimating kinetic parameters.
#'
#' @param obj An `EZbakRFractions` object, which is an `EZbakRData` object on
#' which `EstimateFractions()` has been run.
#' @param strategy Kinetic parameter estimation strategy.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item NSS: Use strategy similar to that presented in Narain et al., 2021 that
#'  assumes you have data which provides a reference for how much RNA was present
#'  at the start of each labeling. In this case, `grouping_factor` and `reference_factor`
#'  must also be set.
#'  \item shortfeed: Estimate kinetic parameters assuming no degradation of labeled
#'  RNA, most appropriate if the metabolic label feed time is much shorter than
#'  the average half-life of an RNA in your system.
#'  \item tilac: Estimate TILAC-ratio as described in Courvan et al., 2022.
#'  \item pulse-chase (NOT YET IMPLEMENTED): Estimate kdeg for a pulse-chase experiment. A kdeg will be estimated
#'  for each time point at which label was present. This includes any pulse-only samples,
#'  as well as all samples including a chase after the pulse.
#'  \item custom (NOT YET IMPLEMENTED): Provide a custom function that takes
#'  fraction estimates as input and produces as output kinetic parameter estimates.
#' }
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param populations Mutational populations that were analyzed to generate the
#' fractions table to use. For example, this would be "TC" for a standard
#' s4U-based nucleotide recoding experiment.
#' @param fraction_design "Design matrix" specifying which RNA populations exist
#' in your samples. By default, this will be created automatically and will assume
#' that all combinations of the `mutrate_populations` you have requested to analyze are
#' present in your data. If this is not the case for your data, then you will have
#' to create one manually. See docs for `EstimateFractions` (run ?EstimateFractions()) for more details.
#' @param repeatID If multiple `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param grouping_factor Which sample-detail columns in the metadf should be used
#' to group samples by for calculating the average RPM at a particular time point?
#' Only relevant if `strategy = "NSS"`.
#' @param reference_factor Which sample-detail column in the metafd should be used to
#' determine which group of samples provide information about the starting levels of RNA
#' for each sample? This should have values that match those in `grouping_factor`.
#' #' Only relevant if `strategy = "NSS"`.
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param feature_lengths Table of effective lengths for each feature combination in your
#' data. For example, if your analysis includes features named GF and XF, this
#' should be a data frame with columns GF, XF, and length.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @return `EZbakRKinetics` object.
#' @import data.table
#' @export
EstimateKinetics <- function(obj,
                             strategy = c("standard", "tilac", "NSS",
                                          "shortfeed",
                                          "pulse-chase"),
                             features = NULL,
                             populations = NULL,
                             fraction_design = NULL,
                             repeatID = NULL,
                             exactMatch = TRUE,
                             grouping_factor = NULL,
                             reference_factor = NULL,
                             character_limit = 20,
                             feature_lengths = NULL,
                             overwrite = TRUE){

  ### Check that input is valid

  # EZbakRData object on which EstimateFractions() has been run?
  if(!methods::is(obj, "EZbakRFractions")){

    if(methods::is(obj, "EZbakRData")){

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

  if(strategy == "tilac"){

    obj <- tilac_ratio_estimation(obj,
                                  features = features,
                                  populations = populations,
                                  fraction_design = fraction_design,
                                  repeatID = repeatID,
                                  exactMatch = exactMatch,
                                  character_limit = character_limit,
                                  feature_lengths = feature_lengths,
                                  overwrite = overwrite)


  }else{

    obj <- Standard_kinetic_estimation(obj,
                                       strategy = strategy,
                                       features = features,
                                       populations = populations,
                                       fraction_design = fraction_design,
                                       reference_factor = reference_factor,
                                       grouping_factor = grouping_factor,
                                       repeatID = repeatID,
                                       exactMatch = exactMatch,
                                       character_limit = character_limit,
                                       feature_lengths = feature_lengths,
                                       overwrite = overwrite)

  }


  return(obj)

}


# kdeg = -log(1 - fn)/tl
# ksyn = (normalized read count)*kdeg
Standard_kinetic_estimation <- function(obj,
                                        strategy = c("standard", "NSS",
                                                     "shortfeed",
                                                     "pulse-chase"),
                                        features = NULL,
                                        populations = NULL,
                                        fraction_design = NULL,
                                        repeatID = NULL,
                                        exactMatch = TRUE,
                                        character_limit = 20,
                                        reference_factor = NULL,
                                        grouping_factor = NULL,
                                        feature_lengths = NULL,
                                        overwrite = TRUE){

  ### Hack to deal with devtools::check() NOTEs
  tl <- kdeg <- log_kdeg <- se_log_kdeg <- ..cols_to_keep <- ..kinetics_cols_to_keep <- NULL
  ksyn <- normalized_reads <- log_ksyn <- se_log_ksyn <- scale_factor <- n <- nolabel_rpm <- NULL
  old_rpm <- new_rpm <- geom_mean <- rpm <- nolabel_n <- nolabel_reps <- NULL
  reference_rpm <- reference_n <- reference_reps <- NULL

  rm(..cols_to_keep)
  rm(..kinetics_cols_to_keep)


  `.` <- list


  strategy <- match.arg(strategy)

  ### Figure out which fraction new estimates to use

  # Function is in Helpers.R
  fractions_name <- EZget(obj,
                          features = features,
                          populations = populations,
                          fraction_design = fraction_design,
                          repeatID = repeatID,
                          exactMatch = exactMatch,
                          returnNameOnly = TRUE)


  # Get fractions
  kinetics <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[["metadata"]][["fractions"]][[fractions_name]][["features"]]


  # Determine which column to use for kinetic parameter estimation
  fraction_cols <- colnames(kinetics)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]


  ### Normalize read counts


  # TO-DO: ALLOW USERS TO JUST USE THE TPM FROM ISOFORM QUANTIFICATION
  reads_norm <- get_normalized_read_counts(obj = obj,
                                           features_to_analyze = features_to_analyze,
                                           fractions_name = fractions_name,
                                           feature_lengths = feature_lengths)



  if(strategy == "standard"){

    ### Estimate kdegs

    kinetics <- setDT(data.table::copy(kinetics))

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


    ### Estimate uncertainty in log(kdeg)

    lkdeg_uncert <- function(fn, se_lfn){

      se_lfn <- ifelse(is.nan(se_lfn),
             1.5,
             se_lfn)

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



  }else if(strategy == "NSS"){

    ### IDEA
    # After working through the mathematical formalism, turns out
    # the ideas I had are identical to those presented inNarain et al., 2021.
    # Best you can do is fit a single exponential model using -s4U and +s4U data
    # to infer the "fraction of old reads that are still present by the end of
    # the label time". Turns out to be identical to estimating the average rate
    # constant for a non-homogeneous Poisson process during the label time. If
    # decay dynamics are non-Poissonian (time to decay is not ~ Exponential()),
    # then memorylessness is out the window and you are always going to be screwed.
    # I get the sense that exponential departure time is something of an approximation
    # or central limit theorem of almost any "standard" departure process. Pretty
    # much will be a good approximation as long as there is a predominantly rate-limiting
    # step in the RNA decay process.
    #
    # Analysis steps:
    # 1) Extract -s4U read counts from relevant source
    # Default is just from fractions object
    # Can also use a read count data frame stored in EZbakRData object
    # 2) Normalize read counts
    # Default is TMM
    # Can also provide a table of scale factors
    # Maybe can also choose not to normalize (technically would be same
    # as providing scale table of 1, but I can autmoate that for user)
    # 3) Estimate kdeg = -log(norm. old reads +s4U / norm. -s4U reads)/tl
    # 4) Estimate ksyn = (fn * norm. total reads +s4U * kdeg)/(1 - exp(-kdeg*tl))

    ### NOTE: a lot of steps 1 and 2 are duplicated from CorrectDropout, so
    ### really should refactor this in the future


    ##### STEP 1: GET -S4U READ COUNTS

    ### Which columns should -s4U samples be grouped by?

    metadf <- obj$metadf


    if(is.null(grouping_factor) | is.null(reference_factor)){

      stop("Grouping factors and reference factors must both be specified if
           using the NSS strategy!")

    }else if(length(grouping_factor) != 1 | length(reference_factor) != 1){

      stop("grouping_factor and reference_factor must each be a single
           column of your metadf!")

    }else if(!(grouping_factor %in% colnames(metadf)) | !(reference_factor %in% colnames(metadf))){

      stop("grouping_factor and reference_factor must both share a name
           with columns in your metadf!")

    }

    # Necessary generalizations:
    # 1) Metadf column used (e.g., pulse-chase)
    reference_data <- kinetics %>%
      dplyr::inner_join(metadf,
                        by = "sample") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(reference_rpm = n/(sum(n)/1000000),
                    reference_n = n) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(grouping_factor, features_to_analyze)))) %>%
      dplyr::summarise(reference_rpm = mean(reference_rpm),
                       reference_n = mean(reference_n),
                       reference_reps = dplyr::n()) %>%
      dplyr::select(!!grouping_factor, !!features_to_analyze,
                    reference_rpm, reference_n, reference_reps)


    ##### STEP 2: INTEGRATE WITH +S4U TO GET ADJUSTED FRACTION NEW
    kinetics <- kinetics %>%
      dplyr::inner_join(metadf %>% dplyr::filter(tl > 0) %>%
                          dplyr::select(sample, tl,
                                        !!reference_factor),
                        by = "sample") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(rpm = n/(sum(n)/1000000)) %>%
      dplyr::inner_join(reference_data %>%
                          dplyr::rename(!!reference_factor := !!dplyr::sym(grouping_factor)),
                        by = c(reference_factor, features_to_analyze)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(old_rpm = (1 - !!dplyr::sym(fraction_of_interest)) * rpm,
                    new_rpm = (!!dplyr::sym(fraction_of_interest)) * rpm) %>%
      dplyr::mutate(kdeg = dplyr::case_when(
        old_rpm >= reference_rpm ~ -log(1 - !!dplyr::sym(fraction_of_interest))/tl,
        .default = -log(old_rpm/reference_rpm)/tl
      ),
      ksyn = new_rpm * kdeg / (1 - exp(-kdeg*tl)),
      log_kdeg = log(kdeg),
      log_ksyn = log(ksyn))


    # Add normalized read counts
    kinetics <- kinetics %>%
      dplyr::inner_join(reads_norm %>%
                          dplyr::select(sample, !!features_to_analyze, "normalized_reads"),
                        by = c("sample", features_to_analyze))


    ##### STEP 3: ESTIMATE UNCERTAINTY
    ##### (higher than in standard due to -s4U read variance)

    ### TO-DO: Incorporate number of -s4U replicates used to infer nolabel_n

    ### Uncertainty approximated with delta method
    lkdeg_uncert_nss <- function(fn, se_lfn, Rs4U, Rctl, tl, kdeg){

      part_1 <- (fn^2)*(se_lfn^2)
      part_2 <- ((1/(Rs4U + 1))^2)*Rs4U/(tl^2)
      part_3 <- ((1/(Rctl + 1))^2)*Rctl/(tl^2)
      kdeg_var <- part_1 + part_2 + part_3
      uncert <- sqrt(((1/kdeg)^2)*kdeg_var)

      return(uncert)

    }

    lksyn_uncert_nss <- function(fn, se_lfn, Rs4U, Rctl, tl, kdeg){

      ### kdeg and log(kdeg) variances
      part_1 <- (fn^2)*(se_lfn^2)
      part_2 <- ((1/(Rs4U + 1))^2)*Rs4U/(tl^2)
      part_3 <- ((1/(Rctl + 1))^2)*Rctl/(tl^2)
      kdeg_var <- part_1 + part_2 + part_3
      lkdeg_var <- sqrt(((1/kdeg)^2)*kdeg_var)

      ### log(ksyn) variance
      part_1s <- ((1 - fn)^2)*(se_lfn^2)
      part_2s <- ((1/(Rs4U + 1))^2)*Rs4U
      part_3s <- lkdeg_var
      part_4s <- ((tl*exp(-tl*kdeg))/(1 - exp(-tl*kdeg)))*kdeg_var
      lksyn_var <- part_1s + part_2s + part_3s + part_4s

      return(sqrt(lksyn_var))


    }

    se_of_interest <- paste0("se_logit_", fraction_of_interest)

    kinetics <- setDT(kinetics)
    kinetics[, se_log_kdeg := lkdeg_uncert_nss(fn = get(fraction_of_interest),
                                               se_lfn = get(se_of_interest),
                                               Rs4U = n,
                                               Rctl = reference_n,
                                               tl = tl,
                                               kdeg = kdeg)]

    # Estimate uncertainty (assuming normalized_reads ~ Poisson(normalized_reads)/scale_factor)
    kinetics[, se_log_ksyn := lksyn_uncert_nss(fn = get(fraction_of_interest),
                                               se_lfn = get(se_of_interest),
                                               Rs4U = n,
                                               Rctl = reference_n,
                                               tl = tl,
                                               kdeg = kdeg)]



  }else if(strategy == "shortfeed"){

    ### IDEA
    # Take it to the limit where there is no degradation (so dR/dt = ks).
    # kdeg = fn / tl
    # ksyn = fn * (normalized read count) / tl; basically new RNA produced per unit time
    # Unclear how much this will matter, as it is literally just the limit of
    # the standard analysis, but worth implementing I think.


    ### NOTE: Currently a near perfect duplication of standard code, just
    ### refactor into a function for goodness sake!

    ### Estimate kdegs

    kinetics <- setDT(data.table::copy(kinetics))

    # Add label time info
    metadf <- obj$metadf

    metadf <- setDT(data.table::copy(metadf))
    setkey(metadf, sample)
    setkey(kinetics, sample)

    kinetics <- kinetics[metadf[,c("sample", "tl")], nomatch = NULL]

    # Make sure no -s4U controls made it through
    kinetics <- kinetics[tl > 0]


    kinetics[, kdeg := get(fraction_of_interest)/tl]
    kinetics[, log_kdeg := log(kdeg)]




    ### Estimate uncertainty in log(kdeg)

    # This needs to be a bit different in this case
    lkdeg_uncert_sf <- function(fn, se_lfn){

      deriv <- (1/fn)*(fn*(1-fn))

      uncert <- abs(deriv)*se_lfn

      return(uncert)

    }

    se_of_interest <- paste0("se_logit_", fraction_of_interest)

    kinetics[, se_log_kdeg := lkdeg_uncert_sf(fn = get(fraction_of_interest),
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


  }



  ##### PROCESS OUTPUT AND RETURN

  # What should output be named?
  kinetics_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")

  if(nchar(kinetics_vect) > character_limit){

    num_kinetics <- length(obj[['kinetics']])
    kinetics_vect <- paste0("kinetics", num_kinetics + 1)

  }


  readcount_vect <- paste0(paste(gsub("_", "", features_to_analyze), collapse = "_"),
                           "_readcount_df")

  if(nchar(readcount_vect) > character_limit){

    num_reads <- sum(grepl("^readcount_df", names(obj[['readcounts']])))
    readcount_vect <- paste0("readcount_df", num_reads + 1)

  }


  # Are there any metadata objects
  if(length(obj[['metadata']][['kinetics']]) > 0){

    kinetics_vect <- decide_output(obj,
                                   proposed_name = kinetics_vect,
                                   type = "kinetics",
                                   features = features_to_analyze,
                                   kstrat = strategy,
                                   overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'kinetics',
                               features = features_to_analyze,
                               kstrat = strategy,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }


  if(length(obj[['metadata']][['readcounts']]) > 0){

    readcount_vect <- decide_output(obj,
                                    proposed_name = readcount_vect,
                                    type = "readcounts",
                                    features = features_to_analyze,
                                    counttype = "TMM_normalized",
                                    overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      rc_repeatID <- 1

    }else{

      rc_repeatID <- length(EZget(obj,
                                  type = 'readcounts',
                                  features = features_to_analyze,
                                  counttype = "TMM_normalized",
                                  returnNameOnly = TRUE,
                                  exactMatch = TRUE,
                                  alwaysCheck = TRUE)) + 1
    }

  }else{

    rc_repeatID <- 1

  }

  obj[["kinetics"]][[kinetics_vect]] <- dplyr::as_tibble(kinetics) %>%
    dplyr::select(sample, !!features_to_analyze, kdeg, log_kdeg, se_log_kdeg, ksyn, log_ksyn, se_log_ksyn, normalized_reads, n)

  # Eventually want to add count matrix output
  obj[["readcounts"]][[readcount_vect]] <- dplyr::as_tibble(reads_norm) %>%
    dplyr::select(sample, !!features_to_analyze, n, normalized_reads, geom_mean, scale_factor)


  obj[["metadata"]][["kinetics"]][[kinetics_vect]] <- list(features = features_to_analyze,
                                                           kstrat = strategy,
                                                           repeatID = repeatID)
  obj[["metadata"]][["readcounts"]][[readcount_vect]] <- list(features = features_to_analyze,
                                                              counttype = "TMM_normalized",
                                                              repeatID = rc_repeatID)


  if(!methods::is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)


}


# Estimate TILAC ratio
tilac_ratio_estimation <- function(obj,
                                   features = NULL,
                                   populations = NULL,
                                   exactMatch = TRUE,
                                   fraction_design = NULL,
                                   grouping_factor = NULL,
                                   repeatID = NULL,
                                   character_limit = 20,
                                   feature_lengths = NULL,
                                   overwrite = TRUE){

  ### Hack to deal with devtools::check() NOTEs
  ..cols_to_keep <- TILAC_ratio <- log_TILAC_ratio <- ..cols_keep <- NULL

  rm(..cols_to_keep)
  rm(..cols_keep)

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

  kinetics <- setDT(data.table::copy(kinetics))

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
                                           features_to_analyze = features_to_analyze,
                                           feature_lengths = feature_lengths)


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

  if(!methods::is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)



}
