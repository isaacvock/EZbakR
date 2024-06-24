#' Import transcript isoform quantification into EZbakRData object
#'
#' A convenient wrapper to \code{tximport()} for importing isoform quantification
#' data into an EZbakRData object. This is necessary to run functions such as
#' \code{EstimateIsoformFractions}.
#'
#' @param obj An `EZbakRData` object.
#' @param files A named vector of paths to all transcript quantification files that you would like
#' to import. This will be passed as the first argument of `tximport::tximport()` (also named `files`).
#' The names of this vector should be the same as the sample names as they appear in the
#' metadf of the `EZbakRData` object.
#' @param quant_tool String denoting the type of software used to generate the abundances. Will
#' get passed to the `type` argument of `tximport::tximport()`. As described in the documentation for
#' `tximport` 'Options are "salmon", "sailfish", "alevin", "piscem", "kallisto", "rsem", "stringtie",
#' or "none". This argument is used to autofill the arguments below (geneIdCol, etc.) "none" means
#' that the user will specify these columns. Be aware that specifying type other than "none" will
#' ignore the arguments below (geneIdCol, etc.)'. Referenced 'arguments below' can be specified
#' as part of `...`.
#' @param ... Additional arguments to be passed to `tximport::tximport()`. Especially relevant if
#' you set `quant_tool` to "none".
#' @export
ImportIsoformQuant <- function(obj, files,
                               quant_tool = c("none", "salmon", "sailfish",
                                                   "alevin", "piscem", "kallisto",
                                                   "rsem", "stringtie"),
                               ...){

  ### Hack to deal with devtools::check() NOTEs
  trascript_id <- NULL

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
#' @param quant_name Name of transcript isoform quantification table to use. Should be stored
#' in the obj$readcounts list under this name. Use `ImportIsoformQuant()` to create
#' this table. If `quant_name` is `NULL`, it will search for tables containing the string
#' "isoform_quant" in their name, as that is the naming convention used by `ImportIsoformQuant()`.
#' If more than one such table exists, an error will be thrown and you will have to specify
#' the exact name in `quant_name`.
#' @param gene_to_transcript Table with columns `transcript_id` and all feature related
#' columns that appear in the relevant fractions table. This is only relevant as a hack to
#' to deal with the case where STAR includes in its transcriptome alignment transcripts
#' on the opposite strand from where the RNA actually originated. This table will be used
#' to filter out such transcript-feature combinations that should not exist.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @export
EstimateIsoformFractions <- function(obj,
                                     features = NULL,
                                     populations = NULL,
                                     fraction_design = NULL,
                                     quant_name = NULL,
                                     gene_to_transcript = NULL,
                                     overwrite = TRUE){


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

    fraction_name <- EZget(obj = obj,
                           features = features,
                           populations = populations,
                           fraction_design = fraction_design,
                           returnNameOnly = TRUE)
  }


  ### Figure out what the gene name is

  fraction <- obj$fractions[[fraction_name]]

  features <- obj[['metadata']][['fractions']][[fraction_name]][['features']]

  fraction <- setDT(fraction)

  # Cute strat: transcript set should have the most unique elements
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


  ##### DETERMINE OUTPUT NAME AND RETURN

  metadata_list <- obj[['metadata']][['fractions']][[fraction_name]]
  fd <- metadata_list[['fraction_design']]
  pops <- metadata_list[['populations']]

  output_name <- "isoforms"

  output_name <- decide_output(obj, output_name,
                                type = "fractions",
                                features = c(gene_colnames, "transcript_id"),
                                populations = pops,
                                fraction_design = fd,
                                overwrite = overwrite)


  obj[['fractions']][[output_name]] <- isoform_fit


  obj[['metadata']][['fractions']][[output_name]] <- list(features = c(gene_colnames, "transcript_id"),
                                                          populations = pops,
                                                          fraction_design = fd)


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

  ### Hack to deal with devtools::check() NOTEs
  n <- fn <- group <- transcript_id <- p <- nreads <- NULL


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
               upper = 9,
               lower = -9,
               prior_a = 0,
               prior_b = 1,
               hessian = TRUE)

  uncertainty <- sqrt(diag(solve(fit$hessian)))

  return(dplyr::tibble(transcript_id = colnames(design_matrix),
                logit_fn = fit$par,
                se_logit_fn = uncertainty))


}


# Process EZbakR data so as to fit beta regression model
Isoform_Fraction_Disambiguation <- function(obj, sample_name,
                                            quant_name,
                                            fraction_name,
                                            gene_colnames,
                                            transcript_colname,
                                            gene_to_transcript){

  ### Hack to deal with devtools::check() NOTEs
  expected_count <- TPM <- n <- group <- effective_length <- transcript_id <- n <- isoform_cnt <- data <- NULL
  fnest <- logit_fn <- se_logit_fn

  ### THINGS THAT NEED TO BE INFERRED
  # 1) Which fractions to grab
  # 2) What the 'set of transcripts' column is called
  # 3) Which quantification to use


  message(paste0("Analyzing sample ", sample_name, "..."))

  # Determine which column is the fraction estimate of interest
  fraction_cols <- colnames(obj$fractions[[fraction_name]])

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]
  logit_fraction_of_interest <- paste0("logit_", fraction_of_interest)
  logit_fraction_se <- paste0("se_logit_", fraction_of_interest)

  if(length(fraction_of_interest) > 1){
    stop("There is more than one non-redundant fraction estimate. Isoform
         deconvolution is currently only compatible with single-label data
         (e.g., s4U-labeled or s6G-labeled, not dually labeled).")
  }

  # Filter for data from sample of interest
  Fns <- obj$fractions[[fraction_name]] %>%
    dplyr::filter(sample == sample_name)

  quant <- obj$readcounts[[quant_name]] %>%
    dplyr::filter(sample == sample_name & expected_count > 10 & TPM > 1)

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
    dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", gene_colnames, "transcript_id")))) %>%
    dplyr::summarise(!!fraction_of_interest := sum(n*(!!sym(fraction_of_interest)))/sum(n),
                     effective_length = mean(effective_length),
                     expected_count = mean(expected_count),
                     !!logit_fraction_se := mean(!!sym(logit_fraction_se))) %>%
    dplyr::mutate(!!logit_fraction_of_interest := logit(!!sym(fraction_of_interest)))


  Fns_multi <- Fns %>%
    dplyr::left_join(single_isoforms %>%
                dplyr::mutate(isoform_cnt = n) %>%
                dplyr::select(-n),
              by = gene_colnames) %>%
    dplyr::filter(is.na(isoform_cnt) & !is.na(expected_count))

  # Run beta regression model to infer isoform-specific fractions

  logit_fraction_of_interest <- paste0("logit_", fraction_of_interest)

  Fns_multi <- Fns_multi %>%
    dplyr::mutate(fn = !!sym(fraction_of_interest)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(gene_colnames))) %>%
    tidyr::nest() %>%
    dplyr::mutate(fnest = lapply(data, fit_beta_regression)) %>%
    dplyr::select(!!gene_colnames, fnest) %>%
    tidyr::unnest(cols = c(fnest)) %>%
    dplyr::mutate(!!logit_fraction_of_interest := logit_fn,
                  !!logit_fraction_se := se_logit_fn,
                  !!fraction_of_interest := inv_logit(logit_fn)) %>%
    dplyr::select(-logit_fn, -se_logit_fn) %>%
    dplyr::inner_join(quant,
                      by = c("transcript_id")) %>%
    dplyr::select(-TPM)


  # Combine single and multi-isoform gene estimates

  output <- Fns_single %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, !!gene_colnames, transcript_id, !!fraction_of_interest, !!logit_fraction_of_interest,
                  !!logit_fraction_se, expected_count, effective_length) %>%
    dplyr::bind_rows(Fns_multi %>% dplyr::ungroup()) %>%
    dplyr::select(sample, !!gene_colnames, transcript_id, everything()) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n = expected_count,
                  RPK = expected_count/(effective_length/1000),
                  !!logit_fraction_of_interest := ifelse(abs(!!sym(logit_fraction_of_interest)) >= 9,
                                                  stats::rnorm(1, !!sym(logit_fraction_of_interest), 0.25),
                                                  !!sym(logit_fraction_of_interest)),
                  !!fraction_of_interest := inv_logit(!!sym(logit_fraction_of_interest))) %>%
    dplyr::select(-expected_count, -effective_length)

  return(output)


}

