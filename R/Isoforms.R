#' Import transcript isoform quantification into EZbakRData object
#'
#' A convenient wrapper to \code{tximport()} for importing isoform quantification
#' data into an EZbakRData object. This is necessary to run functions such as
#' \code{EstimateIsoformFractions}.
#'
ImportIsoformQuant <- function(obj, files,
                               quant_tool = c("none", "salmon", "sailfish",
                                                   "alevin", "piscem", "kallisto",
                                                   "rsem", "stringtie"),
                               ...){

  quant_tool <- match.arg(quant_tool)

  txi <- tximport::tximport(files, type = quant_tool, txIn = TRUE, txOut = TRUE)

  counts_df <- dplyr::as_tibble(txi$counts, rownames = "transcript_id") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "expected_count",
                        cols = !transcript_id)

  final_df <- dplyr::as_tibble(txi$length,
                               rownames = "transcript_id") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "effective_length",
                        cols = !transcript_id) %>%
    dplyr::inner_join(rsem_counts_df,
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
EstimateIsoformFractions <- function(obj,
                                     quant_name = NULL,
                                     fraction_name = NULL){

  if(is.null(quant_name)){

    possible_quant_names <- names(ezbdo$readcounts)

    quant_name <- possible_quant_names[grepl("isoform_quant", possible_quant_names)]

    if(length(quant_name) > 1){

      stop("More than one isoform quantification table exists in your EZbakRData
           object! You need to tell EZbakR which one to use for the isoform-level
           analysis by specifying `quant_name`")

    }

  }

  if(is.null(fraction_name)){

    possible_fraction_names <- names(ezbdo$fractions)

    fraction_name <- possible_fraction_names[grepl("transcripts_",
                                                   possible_fraction_names) |
                                               grepl("_transcripts",
                                                     possible_fraction_names)]

    if(length(fraction_name) > 1){

      stop("More than one transcripts fraction estimates table exists in your
           EZbakRData object! You need to tell EZbakR which one to use for
           isoform-level analysis by specifying `fraction_name")

    }

  }


  ### Estimate fractions
  samp_names <- unique(obj$fractions_GF_transcripts$sample)

  isoform_fit <- purrr::map(.x = samp_names,
                            .f = Isoform_Fraction_Disambiguation,
                            obj = obj,
                            quant_name = quant_name,
                            fraction_name = fraction_name)

  isoform_fit <- dplyr::bind_rows(isoform_fit)

  output_name <-

  obj[['fractions']][['isoforms']] <- isoform_fit

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
    dplyr::mutate(nreads = n,
                  fn = fraction_highTC) %>%
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
                                            quant_name = NULL,
                                            fraction_name = NULL){

  message(paste0("Analyzing sample ", sample_name, "..."))

  ### THINGS THAT NEED TO BE INFERRED
  # 1) Which fractions to grab
  # 2) What the gene ID column is called
  # 3) What the 'set of transcripts' column is called
  # 4) Which quantification to use

  if(is.null(fraction_name)){




  }

  Fns <- obj$fractions[[fraction_name]] %>%
    filter(sample == sample_name)

  quant <- obj$readcounts[[quant_name]] %>%
    filter(sample == sample_name & expected_count > 0)


  Fns <- Fns %>%
    dplyr::mutate(group = 1:n()) %>%
    tidyr::separate_rows(transcripts, sep = "\\+") %>%
    dplyr::mutate(transcript_id = transcripts) %>%
    dplyr::select(-transcripts) %>%
    dplyr::inner_join(quant,
               by = c("transcript_id", "sample"))

  Fns <- Fns %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(p = expected_count/effective_length/sum(expected_count/effective_length))

  # Filter out genes that only have one isoform
  single_isoforms <- Fns %>%
    dplyr::ungroup() %>%
    dplyr::select(GF, transcript_id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(GF) %>%
    dplyr::count() %>%
    dplyr::filter(n == 1)

  Fns_single <- Fns %>%
    dplyr::inner_join(single_isoforms %>% dplyr::select(-n),
                      by = c("GF"))

  Fns_multi <- Fns %>%
    dplyr::left_join(single_isoforms %>%
                dplyr::mutate(isoform_cnt = n) %>%
                dplyr::select(-n),
              by = "GF") %>%
    dplyr::filter(is.na(isoform_cnt))


  Fns_multi <- Fns_multi %>%
    dplyr::group_by(GF) %>%
    tidyr::nest() %>%
    dplyr::mutate(fnest = lapply(data, fit_beta_regression)) %>%
    dplyr::select(GF, fnest) %>%
    tidyr::unnest(cols = c(fnest)) %>%
    dplyr::select(-logit_fn)


  output <- Fns_single %>%
    ungroup() %>%
    dplyr::select(sample, GF, transcript_id, fraction_highTC, logit_fraction_highTC,
                  expected_count, effective_length) %>%
    dplyr::bind_rows(Fns_multi %>% dplyr::ungroup()) %>%
    dplyr::select(sample, GF, transcript_id, everything())

  return(output)


}
