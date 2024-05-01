#' Import transcript isoform quantification into EZbakRData object
#'
#' A convenient wrapper to \code{tximport()} for importing isoform quantification
#' data into an EZbakRData object. This is necessary to run functions such as
#' \code{EstimateIsoformFractions}.
#'
#' @export
ImportIsoformQuant <- function(obj, files,
                               quant_tool = c("none", "salmon", "sailfish",
                                                   "alevin", "piscem", "kallisto",
                                                   "rsem", "stringtie"),
                               ...){

  quant_tool <- match.arg(quant_tool)

  txi <- tximport::tximport(files, type = quant_tool, txIn = TRUE, txOut = TRUE)

  counts_df <- dplyr::as_tibble(txi$counts, rownames = "transcript_id",
                                .name_repair = "unique") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "expected_count",
                        cols = !transcript_id)

  tpm_df <- dplyr::as_tibble(txi$abundance, rownames = "transcript_id",
                             .name_repair = "unique") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "TPM",
                        cols = !transcript_id)

  final_df <- dplyr::as_tibble(txi$length,
                               rownames = "transcript_id",
                               .name_repair = "unique") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "effective_length",
                        cols = !transcript_id) %>%
    dplyr::inner_join(counts_df,
                      by = c("sample", "transcript_id")) %>%
    dplyr::inner_join(tpm_df,
                      by = c("sample", "transcript_id"))



  if(quant_tool == "none"){

    table_name <- "isoform_quant_custom"

  }else{

    table_name <- paste0("isoform_quant_", quant_tool)

  }


  outlist <- list(final_df)
  names(outlist) <- table_name

  if(exists("readcounts", where = obj)){

    obj[['readcounts']] <- outlist

  }else{

    obj[['readcounts']] <- append(obj[['readcounts']], outlist)

  }


  return(obj)

}

#' Estimate isoform-specific fraction news (or more generally "fractions").
#'
#' Combines the output of \code{EstimateFractions} with transcript isoform
#' quantification performed by an outside tool (e.g., RSEM, kallisto, salmon, etc.)
#' to infer transcript isoform-specific fraction news (or more generally fraction
#' of reads coming from a particular mutation population).
#' @param obj An `EZbakRData` object
#' @param quantification Transcript isoform quantification information. This
#' can be provided in multiple formats:
#' \itemize{
#'  \item A data frame with the columns `sample`, `TPM`, `transcript_id`, and `gene_id`.
#'  The `sample` column should contain the names of all samples present in the cB in
#'  the `EZbakRData` object for which fractions have been estimated. `TPM` represents
#'  the transcripts per million estimated by the quantification tool of your choice.
#'  \item A data frame with two columns, one named `sample` and another named `path`.
#'  The `sample` column should contain the names of all samples present in the cB in
#'  the `EZbakRData` object for which fractions have been estimated. The `path`
#'  column should be a path to a transcript isoform quantification file. **If this
#'  form of `quantification` is provided, then `quant_tool` must be specified so
#'  that the quantification files can be properly interpreted
#' }
#' @export
EstimateIsoformFractions <- function(obj,
                                     fractions_identifier = "transcripts",
                                     quant_name = NULL,
                                     fraction_name = NULL,
                                     gene_to_transcript = NULL){


  if(is.null(quant_name)){

    possible_quant_names <- names(obj$readcounts)

    quant_name <- possible_quant_names[grepl("isoform_quant", possible_quant_names)]

    if(length(quant_name) > 1){

      stop("More than one isoform quantification table exists in your EZbakRData
           object! You need to tell EZbakR which one to use for the isoform-level
           analysis by specifying `quant_name`")

    }else if(length(quant_name) == 0){

      stop("There are no detected isoform quantification tables in your EZbakRData
           object! Create one by running `ImportIsoformQuant()`, providing
           paths to quantification files.")

    }

  }

  if(is.null(fraction_name)){

    possible_fraction_names <- names(obj$fractions)

    fraction_name <- possible_fraction_names[(grepl(paste0(gsub("_","",fractions_identifier), "_"),
                                                   possible_fraction_names) |
                                               grepl(paste0("_", gsub("_","",fractions_identifier)),
                                                     possible_fraction_names)) &
                                               !grepl("^isoforms_",
                                                      possible_fraction_names)]

    if(length(fraction_name) > 1){

      stop("More than one transcripts fraction estimates table exists in your
           EZbakRData object! You need to tell EZbakR which one to use for
           isoform-level analysis by specifying `fraction_name")

    }else if(length(fraction_name) == 0){

      stop("There are no detected transcripts fraction estimates table in your
           EZbakRData object! EZbakR tried to auto-detect such a table by looking
           for a table in the `fractions` list with a name that contained the
           string 'transcripts_' or '_transcripts'. If a similar table with a different
           naming convention exists, provide the identifier as the `fractions_identifier` parameter
           for this function. Alternatively, specify the full table name as the
           `fraction_name` parameter.")

    }

  }


  ### Figure out what the gene name is

  fraction <- obj$fractions[[fraction_name]]

  features <- get_features(fraction, objtype = "fractions")

  lens <- rep(0, times = length(features))
  for(f in seq_along(features)){

    lens[f] <- length(unique(fraction[[features[f]]]))

  }

  isoform_feature <- features[which(lens == max(lens))]

  gene_colnames <- features[features != isoform_feature]



  ### Estimate fractions
  metadf <- obj$metadf
  # TO-DO; GENERALIZE FOR PULSE-CHASE
  samp_names <- metadf[['sample']][metadf$tl > 0]



  isoform_fit <- purrr::map(.x = samp_names,
                            .f = Isoform_Fraction_Disambiguation,
                            obj = obj,
                            quant_name = quant_name,
                            fraction_name = fraction_name,
                            gene_colnames = gene_colnames,
                            transcript_colname = isoform_feature,
                            gene_to_transcript = gene_to_transcript)

  isoform_fit <- dplyr::bind_rows(isoform_fit)

  output_name <- paste0(paste(gsub("_","",gene_colnames), collapse = "_"), paste0('_isoforms_', unlist(strsplit(quant_name, "_"))[3]))

  output_name <- paste0("isoforms_",
                        paste(gsub("_", "", gene_colnames), collapse = "_"),
                        "_",
                        gsub("_", "", isoform_feature),
                        "_",
                        unlist(strsplit(quant_name, "_"))[3])

  obj[['fractions']][[output_name]] <- isoform_fit

  return(obj)


}

# Helper functions that I will use on multiple occasions
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))


# Likelihood calculation for beta regression
beta_r_likelihood <- function(data, design_matrix, v, par,
                              prior_a, prior_b){



  mus <- design_matrix %*% inv_logit(par)

  alphas <- mus*v
  betas <- v - alphas

  logl <- sum(stats::dbeta(data, alphas, betas, log = TRUE)) + sum(stats::dnorm(par, prior_a, prior_b, log = TRUE))

  return(-logl)


}

# Infer mixture of isoform fractions using beta regression
fit_beta_regression <- function(data){



  Fns_onegene <- data %>%
    dplyr::mutate(nreads = n) %>%
    dplyr::select(fn, group, transcript_id, p, nreads) %>%
    tidyr::pivot_wider(names_from = transcript_id,
                       values_from = p,
                       values_fill = 0)

  design_matrix <- as.matrix(Fns_onegene %>% dplyr::ungroup() %>%
                               dplyr::select(-group, -fn, -nreads))


  fns <- Fns_onegene$fn
  v <- Fns_onegene$nreads


  fit <- stats::optim(par = rep(0, times = ncol(design_matrix)),
               beta_r_likelihood,
               data = fns,
               design_matrix = design_matrix,
               v = v,
               method = "L-BFGS-B",
               upper = 5,
               lower = -5,
               prior_a = 0,
               prior_b = 1)


  return(tibble(transcript_id = colnames(design_matrix),
                logit_fn = fit$par))


}


# Process EZbakR data so as to fit beta regression model
Isoform_Fraction_Disambiguation <- function(obj, sample_name,
                                            quant_name,
                                            fraction_name,
                                            gene_colnames,
                                            transcript_colname,
                                            gene_to_transcript){


  ### THINGS THAT NEED TO BE INFERRED
  # 1) Which fractions to grab
  # 3) What the 'set of transcripts' column is called
  # 4) Which quantification to use


  message(paste0("Analyzing sample ", sample_name, "..."))

  # Determine which column is the fraction estimate of interest
  fraction_cols <- colnames(obj$fractions[[fraction_name]])

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]
  logit_fraction_of_interest <- paste0("logit_", fraction_of_interest)

  if(length(fraction_of_interest) > 1){
    stop("There is more than one non-redundant fraction estimate. Isoform
         deconvolution is currently only compatible with single-label data
         (e.g., s4U-labeled or s6G-labeled, not dually labeled).")
  }

  # Filter for data from sample of interest
  Fns <- obj$fractions[[fraction_name]] %>%
    filter(sample == sample_name)

  quant <- obj$readcounts[[quant_name]] %>%
    filter(sample == sample_name & expected_count > 10 & TPM > 1)

  # Need to have one row for each transcript ID from a group of
  # transcript IDs, and need to keep track of which reads came from
  # which groups of transcript IDs
  # Also includes a hack to fix issue where a read gets assigned to a transcript
  # isoform on the opposite strand from which the read's RNA originated.
  if(is.null(gene_to_transcript)){

    Fns <- Fns %>%
      dplyr::mutate(group = 1:n()) %>%
      tidyr::separate_rows(!!transcript_colname, sep = "\\+") %>%
      dplyr::mutate(transcript_id = !!sym(transcript_colname)) %>%
      dplyr::select(-!!transcript_colname) %>%
      dplyr::inner_join(quant,
                        by = c("transcript_id", "sample"))

  }else{

    Fns <- Fns %>%
      dplyr::mutate(group = 1:n()) %>%
      tidyr::separate_rows(!!transcript_colname, sep = "\\+") %>%
      dplyr::mutate(transcript_id = !!sym(transcript_colname)) %>%
      dplyr::select(-!!transcript_colname) %>%
      dplyr::inner_join(gene_to_transcript,
                        by = c("transcript_id", gene_colnames)) %>%
      dplyr::inner_join(quant,
                        by = c("transcript_id", "sample"))

  }


  # What fraction of the reads are expected to have come from each isoform?
  Fns <- Fns %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(p = expected_count/effective_length/sum(expected_count/effective_length))

  # Filter out genes that only have one isoform
  single_isoforms <- Fns %>%
    dplyr::ungroup() %>%
    dplyr::select(!!gene_colnames, transcript_id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(gene_colnames))) %>%
    dplyr::count() %>%
    dplyr::filter(n == 1)

  Fns_single <- Fns %>%
    dplyr::inner_join(single_isoforms %>% dplyr::select(-n),
                      by = gene_colnames) %>%
    dplyr::group_by(across(all_of(c("sample", gene_colnames, "transcript_id")))) %>%
    dplyr::summarise(!!fraction_of_interest := sum(n*(!!sym(fraction_of_interest)))/sum(n),
                     effective_length = mean(effective_length),
                     expected_count = mean(expected_count)) %>%
    dplyr::mutate(!!logit_fraction_of_interest := logit(!!sym(fraction_of_interest)))


  Fns_multi <- Fns %>%
    dplyr::left_join(single_isoforms %>%
                dplyr::mutate(isoform_cnt = n) %>%
                dplyr::select(-n),
              by = gene_colnames) %>%
    dplyr::filter(is.na(isoform_cnt) & !is.na(expected_count))

  # Run beta regression model to infer isoform-specific fractions
  logit_fraction_of_interest <- paste0("logit_", fraction_of_interest)


  # WAY I WANT TO DO IT
  Fns_multi <- Fns_multi %>%
    dplyr::mutate(fn = !!sym(fraction_of_interest)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(gene_colnames))) %>%
    tidyr::nest() %>%
    dplyr::mutate(fnest = lapply(data, fit_beta_regression)) %>%
    dplyr::select(!!gene_colnames, fnest) %>%
    tidyr::unnest(cols = c(fnest)) %>%
    dplyr::mutate(!!logit_fraction_of_interest := logit_fn,
                  !!fraction_of_interest := inv_logit(logit_fn)) %>%
    dplyr::select(-logit_fn) %>%
    dplyr::inner_join(quant,
                      by = c("transcript_id")) %>%
    dplyr::select(-TPM)

  # # CHESS-specific hack for right now
  # Fns_multi <- Fns_multi %>%
  #   dplyr::ungroup() %>%
  #   dplyr::rowwise() %>%
  #   dplyr::filter(grepl(paste0(GF, "\\."), transcript_id)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(fn = !!sym(fraction_of_interest)) %>%
  #   dplyr::group_by(dplyr::across(dplyr::all_of(gene_colnames))) %>%
  #   tidyr::nest() %>%
  #   dplyr::mutate(fnest = lapply(data, fit_beta_regression)) %>%
  #   dplyr::select(!!gene_colnames, fnest) %>%
  #   tidyr::unnest(cols = c(fnest)) %>%
  #   dplyr::mutate(!!logit_fraction_of_interest := logit_fn,
  #                 !!fraction_of_interest := inv_logit(logit_fn)) %>%
  #   dplyr::select(-logit_fn) %>%
  #   dplyr::inner_join(quant,
  #                     by = c("transcript_id"))


  # Combine single and multi-isoform gene estimates

  output <- Fns_single %>%
    ungroup() %>%
    dplyr::select(sample, !!gene_colnames, transcript_id, !!fraction_of_interest, !!logit_fraction_of_interest,
                  expected_count, effective_length) %>%
    dplyr::bind_rows(Fns_multi %>% dplyr::ungroup()) %>%
    dplyr::select(sample, !!gene_colnames, transcript_id, everything()) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n = expected_count,
                  RPK = expected_count/(effective_length/1000),
                  !!logit_fraction_of_interest := ifelse(abs(!!sym(logit_fraction_of_interest)) >= 5,
                                                  rnorm(1, !!sym(logit_fraction_of_interest), 0.25),
                                                  !!sym(logit_fraction_of_interest)),
                  !!fraction_of_interest := inv_logit(!!sym(logit_fraction_of_interest))) %>%
    dplyr::select(-expected_count, -effective_length)

  return(output)


}

