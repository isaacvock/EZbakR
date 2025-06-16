#' `EZbakRKinetics` object constructor
#'
#' \code{new_EZbakRKinetics} efficiently creates an object of class `EZbakRKinetics`.
#' It does not perform any rigorous checks of the legitimacy of this object.
#' @param kinetics Data frame containing information about the kinetic parameters
#' of interest for each set of features tracked.
#' @param features Features tracked in `kinetics` data frame. Needs to be specified
#' explicitly as it cannot be automatically inferred.
#' @param metadf Data frame describing each of the samples included
#' @param name Optional; name to give to fractions table.
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' in `kinetics`
new_EZbakRKinetics <- function(kinetics, features, metadf,
                               name = NULL,
                               character_limit = 20){
  stopifnot(is.data.frame(kinetics))
  stopifnot(is.data.frame(metadf))

  if(is.null(name)){

    # Can't have '_' in feature names when naming the kinetics object
    features_gsub <- gsub("_", "", features)

    # Name for output
    kinetics_name <- paste(features_gsub, collapse = "_")

    if(nchar(kinetics_name) > character_limit){

      kinetics_name <- "kinetics1"

    }

  }else{

    kinetics_name <- name

  }



  # Make output
  struct_list <- list(kinetics = list(kinetics),
                      metadf = metadf)
  names(struct_list[['kinetics']]) <- kinetics_name

  struct_list[['metadata']][['kinetics']][[kinetics_name]] <- list(features = features,
                                                                    kstrat = "custom")

  structure(struct_list, class = "EZbakRKinetics")
}

#' `EZbakRKinetics` object validator
#'
#' \code{validate_EZbakRKinetics} ensures that input for `EZbakRKinetics` construction
#' is valid.
#' @param obj An object of class `EZbakRKinetics`
#' @param features Features tracked in `kinetics` data frame. Needs to be specified
#' explicitly as it cannot be automatically inferred.
validate_EZbakRKinetics <- function(obj, features){


  ### Extract objects and info from objects

  vals <- unclass(obj)

  kinetics <- dplyr::as_tibble(vals[['kinetics']][[1]])
  metadf <- dplyr::as_tibble(vals[['metadf']])

  kinetics_cols <- colnames(kinetics)
  metadf_cols <- colnames(metadf)


  ### Which columns represent kinetic parameters?
  param_cols <- kinetics_cols[!(kinetics_cols %in% c("sample", "n", features,
                                                     "normalized_reads"))]


  ### Are all supposed features actually in the kinetics table?

  if(!all(features %in% kinetics_cols)){

    stop("Not all features are columns in kinetics!")

  }


  ### Are all kinetic parameter columns numeric?

  param_check <- sapply(kinetics %>% dplyr::select(!!param_cols),
                        function(x) all(is.numeric(x)))

  if(!all(param_check)){

    stop("Not all kinetic parameter columns are exclusively numeric!")

  }


  ### Are there uncertainties for each kinetic parameter?

  se_cols <- param_cols[grepl("^se_", param_cols)]
  est_cols <- param_cols[!(param_cols %in% se_cols)]

  expected_ests <- gsub("^se_", "", se_cols)

  if(length(expected_ests) > length(est_cols)){

    stop("There are more se_* columns than there are parameter estimation
         columns. Each parameter estimate column should have a corresponding
         se_* column, named identically to the parameter estimate column but
         with 'se_' appended to the front. Make sure to not start any parameter
         estimate columns with the substring 'se_'.")

  }

  missing_ses <- est_cols[!(est_cols %in% expected_ests)]

  if(length(missing_ses) > 0){

    message("Kinetic parameter estimates missing corresponding uncertainty
          quantification will have an imputed uncertainty of 0.")


    missing_ses <- paste0("se_", missing_ses)

    for(m in missing_ses){

      kinetics[[m]] <- 0

    }

  }


  ### Does kinetics data frame contain columns named "sample" and "n"?
  if(!("sample" %in% kinetics_cols & "n" %in% kinetics_cols)){

    rlang::abort(
      "kinetics must include columns named sample (representing the sample ID)
      and n (representing the number of reads)!",
      class = "kinetics_sample_n"
    )

  }



  ### Does metadf contain sample and correct tl(s)?

  if(!("sample" %in% metadf_cols)){
    rlang::abort(
      "metadf must include a column named sample (representing the sample ID)!",
      class = "metadf_sample"
    )
  }



  ### Check if all of the samples in kinetics are also in metadf

  samps_kinetics <- unique(kinetics$sample)

  if(!all(samps_kinetics %in% metadf$sample)){
    rlang::abort("Not all samples in kinetics are present in the metadf!",
                 class = "kinetics_metadf_samples")
  }


  ### Check that n in kinetics is numeric and > 0

  is_pos_whole2 <- function(x){

    if(all(is.numeric(x))){

      bool <- all((floor(x) == x) & (x > 0))

    }else{

      bool <- FALSE

    }

    return(bool)

  }


  if(!all(sapply(kinetics[,"n"], is_pos_whole2))){

    rlang::abort("Not all columns of kinetics tracking counts of reads contain positive whole numbers
         strictly greating than 0!",
                 class = "kinetics_n_pos_whole")

  }


  ### Check that label times are numeric and >= 0


  is_pos_num <- function(x){

    if(all(is.numeric(x))){

      bool <- all(x >= 0)


    }else{

      bool <- FALSE

    }


    return(bool)

  }


  ### Infer and check label time columns

  tl_cols <- metadf_cols[grepl("^tl$", metadf_cols) |
                           grepl("^tpulse", metadf_cols) |
                           grepl("^tchase", metadf_cols)]


  if(!all(sapply(metadf[,tl_cols], is_pos_num))){

    rlang::abort("Not all columns of metadf representing label times are numbers
         greater than or equal to 0.",
                 class = "tl_pos_num")

  }


  ### Convert all other metadf columns to factors

  cols_to_convert <- metadf_cols[!(metadf_cols %in% c(tl_cols, "sample"))]

  metadf <- as.data.frame(metadf)

  metadf[cols_to_convert] <- lapply(metadf[cols_to_convert], as.factor)

  metadf <- dplyr::as_tibble(metadf)

  obj$metadf <- metadf


  kinetics_name <- names(obj[['metadata']][['kinetics']])

  obj[['metadata']][['kinetics']][[kinetics_name]] <- list(features = features,
                                                           kstrat = "custom")

  obj[['kinetics']][[kinetics_name]] <- kinetics

  return(obj)




}


#' `EZbakRKinetics` helper function for users
#'
#' \code{EZbakRKinetics} creates an object of class `EZbakRKinetics` and checks the validity
#' of the provided input.
#' @param kinetics Data frame with the following columns:
#' \itemize{
#'  \item sample: Name given to particular sample from which data was collected.
#'  \item features: Any columns that cannot be interpreted as a mutation count
#'  or base nucleotide count (and that aren't named `sample` or `n`) will be
#'  interpreted as an ID for a genomic "feature" from which a read originated.
#'  Common examples of features and typical column names for said features include:
#'  \itemize{
#'    \item Genes; common column names: gene, gene_id, gene_name, GF
#'    \item Genes-exonic; common column names: gene_exon, gene_id_exon, gene_name_exon, XF
#'    \item Transcripts; common column names: transcripts, TF
#'    \item Exonic bins; common column names: exonic_bins, EF, EB
#'    \item Exons; common column names: exons, exon_ids
#'  }
#'  In some cases, a read will often map to multiple features (e.g., exons). Many
#'  functions in bakR expect each of the feature IDs in these cases to be separated
#'  by `+`. For example, if a read overlaps with two exons, with IDs exon_1 and exon_2,
#'  then the corresponding entry in a  column of exonic assignments would be "exon_1+exon_2".
#'  The default expectation can be overwritten though and is thus not strictly enforced.
#'  \item n: Number of reads with identical values for all other columns.
#'  \item kinetic parameter estimates: These can be named whatever you would like as long as
#'  they do not start with the string "se_". This should be reserved for kinetic parameter
#'  uncertainties, if provided.
#'  \item kinetic parameter uncertainties: Uncertainty in your kinetic parameter estimates.
#'  These should be named "se_" followed by the kinetic parameter as its name appears in
#'  the relevant column name of the `kinetics` table.
#' }
#' @param metadf Data frame detailing various aspects of each of the samples included
#' in the kinetics data frame. This includes:
#' \itemize{
#'  \item `sample`: The sample ID, which should correspond to a sample ID in the provided kinetics data frame.
#'  \item `tl`: Metabolic label time. There are several edge cases to be aware of:
#'  \itemize{
#'    \item If more than one metabolic label was used in the set of samples described
#'    by the metadf (e.g., s4U and s6G were used), then the `tl` column should be
#'    replaced by `tl_<muttype>`, where `<muttype>` represents the corresponding mutation
#'    type referenced in the fractions that the label whose incubation time will be listed
#'    in this column. For example, if feeding with s4U in some samples and s6G in others,
#'    then performing standard nucleotide recoding chemistry, you will include
#'    `tl_TC` and `tl_GA` columns corresponding to the s4U and s6G label times, respectively.
#'    \item If a pulse-chase experimental design was used (!!this is strongly discouraged
#'    unless you have a legitimate reason to prefer this design to a pulse-label
#'    design!!), then you should have columns named `tpulse` and `tchase`, corresponding
#'    to the pulse and chase times respectively. The same _<muttype> convention should
#'    be used in the case of multi-label pulse-chase designs.
#'  }
#' \item sample characteristics: The remaining columns can be named whatever you like
#' and should include distinguishing features of groups of samples. Common columns might
#' include:
#' \itemize{
#'  \item `treatment`: The experimental treatment applied to a set of samples.
#'  This could represent things like genetic knockouts or knockdowns, drug treatments, etc.
#'  \item `batch`: An ID for sets of samples that were collected and/or processed together.
#'  Useful for regressing out technical batch effects
#'  }
#' \item `assay`: This optional column should include a string that
#' describes the type of experiment that was done so as to influence
#' how EZbakR analyzes and interprets the data from those samples.
#' Possible values for `assay` currently include:
#' \itemize{
#'  \item standard: Refers to the "standard" nucleotide recoding RNA-seq methods
#'  (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.), in which cells are fed with a
#'  single metabolic label, RNA is extracted and sequenced, and mutations of a particular
#'  type are counted
#'  \item STL: Refers to Start-TimeLapse-seq, a method combining Start-seq (developed
#'  by Karen Adelman's lab) with TimeLapse-seq. Used to infer the kinetics of
#'  transcription initiation and promoter-proximal pause site departure.
#'  \item TT: Refers to Transient-Transcriptome NR-seq, a method combining TT-seq (developed
#'  by Patrick Cramer's lab) with NR-seq. TT-seq involves biochemically enriching
#'  for labeled RNA. By combining this method with nucleotide recoding chemistry (as was
#'  first done by the Simon lab with TT-TimeLapse-seq and has since been done with SLAM
#'  chemistry, often referred to as TTchem-seq), it is possible to bioinformatically
#'  filter out reads coming from unlabeled RNA background.
#'  \item TILAC: Refers to TILAC, a method developed by the Simon lab to achieve
#'  spike-in free normalization of RNA-seq data through the use of a dual labeling
#'  approach inspired by the proteomic method SILAC.
#'  \item subcellular: Refers to techniques such as subcellular TimeLapse-seq (developed
#'  by Stirling Churchman's lab) which combine subcellular fractionation with
#'  NR-seq to infer additional kinetic parameters.
#'  \item sc: Refers to single-cell RNA-seq implementations of NR-seq.
#'  }
#'
#' }
#' @param features Features tracked in `kinetics` data frame. Needs to be specified
#' explicitly as it cannot be automatically inferred.
#' @param name Optional; name to give to fractions table.
#' @param character_limit If name is chosen automatically, limit on the number of
#' characters in said name. If default name yields a string longer than this,
#' then kinetics table will be named `kinetics1`
#' @return An EZbakRKinetics object. This is simply a list of the provided `kinetics` and
#' `metadf` with class `EZbakRKinetics`
#' @examples
#' # Simulate data
#' simdata <- EZSimulate(30)
#'
#' # Get kinetics table by estimating (for demonstration)
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#' ezbdo <- EstimateFractions(ezbdo)
#' ezbdo <- EstimateKinetics(ezbdo)
#' kinetics <- EZget(ezbdo, type = "kinetics")
#'
#' # Create EZbakRKinetics object
#' ezbko <- EZbakRKinetics(kinetics, simdata$metadf, features = "feature")
#'
#' @export
EZbakRKinetics <- function(kinetics, metadf, features,
                            name = NULL,
                            character_limit = 20){

  kinetics <- dplyr::as_tibble(kinetics)
  metadf <- dplyr::as_tibble(metadf)

  validate_EZbakRKinetics(new_EZbakRKinetics(kinetics, metadf, features = features,
                                               name = name,
                                               character_limit = character_limit),
                           features = features)

}
