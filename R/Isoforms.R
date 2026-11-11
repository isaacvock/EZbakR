#' Import transcript isoform quantification into EZbakRData object
#'
#' A convenient wrapper to [\code{tximport()}](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) for importing isoform quantification
#' data into an EZbakRData object. You need to run this before running
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
#' @param txIn Whether or now you are providing isoform level quantification files.
#' Alternative (`txIn = FALSE`) is gene-level quantification. In `ImportIsoformQuant`,
#' `txIn` gets passed to BOTH the `txIn` and `txOut` parameters in `tximport()`.
#' @param ... Additional arguments to be passed to `tximport::tximport()`. Especially relevant if
#' you set `quant_tool` to "none".
#' @return An `EZbakRData` object with an additional element in the `readcounts`
#' list named "isform_quant_<quant_tool>". It contains TPM, expected_count,
#' and effective length information for each transcript_id and each sample.
#' @examples
#'
#' # Dependencies for example
#' library(dplyr)
#' library(data.table)
#'
#' # Simulate and analyze data
#' simdata <- EZSimulate(30, nreps = 2)
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' # Hack to generate example quantification files
#' savedir <- tempdir()
#' rsem_data <- tibble(
#'   transcript_id = paste0("tscript_feature", 1:30),
#'   gene_id = paste0("feature", 1:30),
#'   length = 1000,
#'   effective_length = 1000,
#'   expected_count = 1000,
#'   TPM = 10,
#'   FPKM = 10,
#'   IsoPct = 1
#' )
#'
#' fwrite(rsem_data, file.path(savedir, "Sample_1.isoforms.results"), sep = '\t')
#'
#' files <- file.path(savedir,"Sample_1.isoforms.results")
#' names(files) <- "Sample_1"
#'
#' # Read in file
#' ezbdo <- ImportIsoformQuant(ezbdo, files, quant_tool = "rsem")
#'
#'
#'
#' @export
ImportIsoformQuant <- function(obj, files,
                               quant_tool = c("none", "salmon", "sailfish",
                                                   "alevin", "piscem", "kallisto",
                                                   "rsem", "stringtie"),
                               txIn = TRUE,
                               ...){

  ### Hack to deal with devtools::check() NOTEs
  transcript_id <- NULL

  quant_tool <- match.arg(quant_tool)

  txi <- tximport::tximport(files, type = quant_tool, txIn = txIn, txOut = txIn,
                            ...)

  quant_feature <- ifelse(txIn,
                          "transcript_id",
                          "gene_id")

  counts_df <- dplyr::as_tibble(txi$counts, rownames = quant_feature,
                                .name_repair = "unique") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "expected_count",
                        cols = !dplyr::sym(quant_feature))

  tpm_df <- dplyr::as_tibble(txi$abundance, rownames = quant_feature,
                             .name_repair = "unique") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "TPM",
                        cols = !dplyr::sym(quant_feature))

  final_df <- dplyr::as_tibble(txi$length,
                               rownames = quant_feature,
                               .name_repair = "unique") %>%
    tidyr::pivot_longer(names_to = "sample",
                        values_to = "effective_length",
                        cols = !dplyr::sym(quant_feature)) %>%
    dplyr::inner_join(counts_df,
                      by = c("sample", quant_feature)) %>%
    dplyr::inner_join(tpm_df,
                      by = c("sample", quant_feature))

  if(txIn){

    if(quant_tool == "none"){

      table_name <- "isoform_quant_custom"

    }else{

      table_name <- paste0("isoform_quant_", quant_tool)

    }


  }else{

    if(quant_tool == "none"){

      table_name <- "gene_quant_custom"

    }else{

      table_name <- paste0("gene_quant_", quant_tool)

    }


  }




  outlist <- list(final_df)
  names(outlist) <- table_name

  if(exists("readcounts", where = obj)){

    obj[['readcounts']] <- append(obj[['readcounts']], outlist)



  }else{

    obj[['readcounts']] <- outlist

  }


  return(obj)

}

#' Estimate isoform-specific fraction news (or more generally "fractions").
#'
#' Combines the output of \code{EstimateFractions} with transcript isoform
#' quantification performed by an outside tool (e.g., RSEM, kallisto, salmon, etc.)
#' to infer transcript isoform-specific fraction news (or more generally fraction
#' of reads coming from a particular mutation population).
#'
#' \code{EstimateIsoformFractions} expects as input a "fractions" table with estimates
#' for transcript equivalence class (TEC) fraction news. A transcript equivalence class
#' is the set of transcript isoforms with which a sequencing read is compatible.
#' [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR) is able to assign
#' reads to these equivalence classes so that EZbakR can estimate the fraction of
#' reads in each TEC that are from labeled RNA.
#'
#' \code{EstimateIsoformFractions} estimates transcript isoform fraction news
#' by fitting a linear mixing model to the TEC fraction new estimates + transcript
#' isoform abundance estimates. In other words, each TEC fraction new (data) is modeled
#' as a weighted average of transcript isoform fraction news (parameters to estimate),
#' with the weights determined by the relative abundances of the transcript isoforms
#' in the TEC (data). The TEC fraction new is modeled as a Beta distribution with mean
#' equal to the weighted transcript isoform fraction new average and variance related
#' to the number of reads in the TEC.
#'
#' Transcript isoforms with estimated TPMs less than with an estimated
#' TPM greater than `TPM_min` (1 by default) or an estimated number of read
#' counts less than `count_min` (10 by default) are removed from TECs and will
#' not have their fraction news estimated.
#'
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
#' @param repeatID If multiple `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param fraction_name Name of fraction estimate table to use. Should be stored in the
#' `obj$fractions` list under this name. Can also rely on specifying `features` and/or `populations`
#' and having `EZget()` find it.
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
#' @param TPM_min Minimum TPM for a transcript to be kept in analysis.
#' @param count_min Minimum expected_count for a transcript to be kept in analysis.
#' @return An `EZbakRData` object with an additional table under the "fractions"
#' list. Has the same form as the output of `EstimateFractions()`, and will have the
#' feature column "transcript_id".
#' @examples
#'
#' # Load dependencies
#' library(dplyr)
#'
#' # Simulates a single sample worth of data
#' simdata_iso <- SimulateIsoforms(nfeatures = 25)
#'
#' # We have to manually create the metadf in this case
#' metadf <- tibble(sample = 'sampleA',
#'                      tl = 4,
#'                      condition = 'A')
#'
#' ezbdo <- EZbakRData(simdata_iso$cB,
#'                     metadf)
#'
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' ### Hack in the true, simulated isoform levels
#' reads <- simdata_iso$ground_truth %>%
#'   dplyr::select(transcript_id, true_count, true_TPM) %>%
#'   dplyr::mutate(sample = 'sampleA',
#'                 effective_length = 10000) %>%
#'   dplyr::rename(expected_count = true_count,
#'                 TPM = true_TPM)
#'
#' # Name of table needs to have "isoform_quant" in it
#' ezbdo[['readcounts']][['simulated_isoform_quant']] <- reads
#'
#' ### Perform deconvolution
#' ezbdo <- EstimateIsoformFractions(ezbdo)
#'
#' @export
EstimateIsoformFractions <- function(obj,
                                     features = NULL,
                                     populations = NULL,
                                     fraction_design = NULL,
                                     repeatID = NULL,
                                     exactMatch = TRUE,
                                     fraction_name = NULL,
                                     quant_name = NULL,
                                     gene_to_transcript = NULL,
                                     overwrite = TRUE,
                                     TPM_min = 1,
                                     count_min = 10){


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
                           exactMatch = exactMatch,
                           repeatID = repeatID,
                           returnNameOnly = TRUE)

    if(length(fraction_name) > 1){

      stop("More than one fractions table fits search criterion!
           Either specify additional criterion (i.e., features, populations,
           or fraction_design) or explicitly specify the table name (fraction_name).")
    }
  }


  ### Figure out what the gene name is

  fraction <- obj$fractions[[fraction_name]]

  features <- obj[['metadata']][['fractions']][[fraction_name]][['features']]

  fraction <- data.table::as.data.table(fraction)

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
                            gene_to_transcript = gene_to_transcript,
                            TPM_min = TPM_min,
                            count_min = count_min)

  isoform_fit <- dplyr::bind_rows(isoform_fit)


  ##### DETERMINE OUTPUT NAME AND RETURN

  metadata_list <- obj[['metadata']][['fractions']][[fraction_name]]
  fd <- metadata_list[['fraction_design']]
  pops <- metadata_list[['populations']]

  output_name <- "isoforms"

  # How many identical tables already exist?
  if(overwrite){

    repeatID <- 1

  }else{

    repeatID <- length(EZget(obj,
                                   type = 'fractions',
                                   features = c(gene_colnames, "transcript_id"),
                                   populations = pops,
                                   fraction_design = fd,
                                   returnNameOnly = TRUE,
                                   exactMatch = TRUE,
                                   alwaysCheck = TRUE)) + 1
  }

  # Get name for output
  output_name <- decide_output(obj, output_name,
                                type = "fractions",
                                features = c(gene_colnames, "transcript_id"),
                                populations = pops,
                                fraction_design = fd,
                                overwrite = overwrite)


  obj[['fractions']][[output_name]] <- isoform_fit


  obj[['metadata']][['fractions']][[output_name]] <- list(features = c(gene_colnames, "transcript_id"),
                                                          populations = pops,
                                                          fraction_design = fd,
                                                          repeatID = repeatID)


  return(obj)


}

# Helper functions that I will use on multiple occasions
logit <- function(x) stats::qlogis(x)
inv_logit <- function(x) stats::plogis(x)


# Likelihood calculation for beta regression
beta_r_likelihood <- function(data, design_matrix, v, par,
                              prior_a, prior_b){

  data_clamped <- pmin(pmax(data, inv_logit(-9)), inv_logit(9))

  mus <- as.vector(design_matrix %*% inv_logit(par))
  mus <- pmin(pmax(mus, 1e-12), 1 - 1e-12)

  alphas <- mus*v
  betas <- v - alphas

  logl <- sum(stats::dbeta(data_clamped, alphas, betas, log = TRUE)) + sum(stats::dnorm(par, prior_a, prior_b, log = TRUE))

  return(-logl)


}

# Infer mixture of isoform fractions using beta regression
fit_beta_regression <- function(data){

  cols <- colnames(data)
  deconvolved_feature <- if("transcript_id" %in% cols) "transcript_id" else "gene_id"

  grp_vals <- data$group
  u_grp <- unique(grp_vals)
  feat_vals <- as.character(data[[deconvolved_feature]])
  u_feat <- unique(feat_vals)
  n_g <- length(u_grp)
  n_f <- length(u_feat)

  design_matrix <- matrix(0, nrow = n_g, ncol = n_f, dimnames = list(NULL, u_feat))
  r_idx <- match(grp_vals, u_grp)
  c_idx <- match(feat_vals, u_feat)
  design_matrix[cbind(r_idx, c_idx)] <- as.numeric(data$p)

  first_row <- match(u_grp, grp_vals)
  fns <- as.numeric(data$fn[first_row])
  v <- as.numeric(data$n[first_row])

  prior_a <- 0
  prior_b <- 1.25
  inv_pb2 <- 1 / (prior_b * prior_b)
  prior_const <- -0.5 * n_f * log(2 * pi) - n_f * log(prior_b)

  y <- pmin(pmax(fns, inv_logit(-9)), inv_logit(9))
  log_y <- log(y)
  log_1my <- log1p(-y)
  logit_y <- log_y - log_1my
  lgamma_v <- lgamma(v)

  fn_beta <- function(par){
    theta <- stats::plogis(par)
    mus_raw <- as.vector(design_matrix %*% theta)
    mus <- pmin(pmax(mus_raw, 1e-12), 1 - 1e-12)
    alphas <- mus * v
    betas <- v - alphas
    ll_data <- sum(lgamma_v - lgamma(alphas) - lgamma(betas) +
                     (alphas - 1) * log_y + (betas - 1) * log_1my)
    ll_prior <- prior_const - 0.5 * sum((par - prior_a)^2) * inv_pb2
    -(ll_data + ll_prior)
  }

  gr_beta <- function(par){
    theta <- stats::plogis(par)
    mus_raw <- as.vector(design_matrix %*% theta)
    unclamped <- (mus_raw > 1e-12) & (mus_raw < (1 - 1e-12))
    mus <- pmin(pmax(mus_raw, 1e-12), 1 - 1e-12)
    alphas <- mus * v
    betas <- v - alphas
    d_mu <- ifelse(unclamped,
                   v * (logit_y - digamma(alphas) + digamma(betas)),
                   0)
    xt_dmu <- as.vector(crossprod(design_matrix, d_mu))
    -theta * (1 - theta) * xt_dmu + (par - prior_a) * inv_pb2
  }

  tryCatch(
    {
      fit <- stats::optim(par = rep(0, times = n_f),
                          fn = fn_beta,
                          gr = gr_beta,
                          method = "L-BFGS-B",
                          upper = 9,
                          lower = -9,
                          hessian = FALSE)
    },
    error = function(e){
      message(paste0("Estimation of isoform fraction news failed for
                     features: [",
                     paste0(data[[deconvolved_feature]], collapse = ", "), "]"))
      stop("Error from optim: ", e)
    }
  )

  # Exact analytical Hessian at converged par (fallback to numerical if non-invertible)
  par_hat <- fit$par
  theta <- stats::plogis(par_hat)
  s_vec <- theta * (1 - theta)
  mus_raw <- as.vector(design_matrix %*% theta)
  unclamped <- (mus_raw > 1e-12) & (mus_raw < (1 - 1e-12))
  mus <- pmin(pmax(mus_raw, 1e-12), 1 - 1e-12)
  alphas <- mus * v
  betas <- v - alphas
  d_mu <- ifelse(unclamped,
                 v * (logit_y - digamma(alphas) + digamma(betas)),
                 0)
  w_mu <- ifelse(unclamped,
                 (v * v) * (trigamma(alphas) + trigamma(betas)),
                 0)
  xt_dmu <- as.vector(crossprod(design_matrix, d_mu))
  xtwx <- crossprod(design_matrix, w_mu * design_matrix)
  H_mat <- (s_vec %o% s_vec) * xtwx
  diag(H_mat) <- diag(H_mat) - s_vec * (1 - 2 * theta) * xt_dmu + inv_pb2

  uncertainty <- tryCatch(
    sqrt(abs(diag(solve(H_mat)))),
    error = function(e){
      fit_num <- stats::optim(par = par_hat,
                              fn = beta_r_likelihood,
                              data = fns,
                              design_matrix = design_matrix,
                              v = v,
                              method = "L-BFGS-B",
                              upper = 9,
                              lower = -9,
                              prior_a = prior_a,
                              prior_b = prior_b,
                              hessian = TRUE)
      sqrt(abs(diag(solve(fit_num$hessian))))
    }
  )

  if(deconvolved_feature == "transcript_id"){
    return(data.table::data.table(transcript_id = u_feat,
                                  logit_fn = par_hat,
                                  se_logit_fn = uncertainty))
  }else{
    return(data.table::data.table(gene_id = u_feat,
                                  logit_fn = par_hat,
                                  se_logit_fn = uncertainty))
  }

}


# Fast batch wrapper avoiding tidyr::nest() + tidyr::unnest() overhead across thousands of genes
fit_beta_regression_batch <- function(df_multi, group_cols, keep_cols, deconvolved_feature){
  if(nrow(df_multi) == 0L){
    empty_out <- df_multi[0L, keep_cols, drop = FALSE]
    empty_out[[deconvolved_feature]] <- character(0)
    empty_out$logit_fn <- numeric(0)
    empty_out$se_logit_fn <- numeric(0)
    return(dplyr::as_tibble(empty_out))
  }
  grouped <- dplyr::group_by(df_multi, dplyr::across(dplyr::all_of(group_cols)))
  g_keys <- dplyr::group_keys(grouped)
  g_rows <- dplyr::group_rows(grouped)
  G <- nrow(g_keys)

  res_list <- vector("list", G)
  n_per_g <- integer(G)
  for(g in seq_len(G)){
    sub_df <- df_multi[g_rows[[g]], , drop = FALSE]
    fit_dt <- fit_beta_regression(sub_df)
    res_list[[g]] <- fit_dt
    n_per_g[g] <- nrow(fit_dt)
  }

  bound_fits <- data.table::rbindlist(res_list)
  rep_idx <- rep.int(seq_len(G), n_per_g)
  keys_expanded <- g_keys[rep_idx, keep_cols, drop = FALSE]
  dplyr::bind_cols(keys_expanded, dplyr::as_tibble(bound_fits))
}


# Process EZbakR data so as to fit beta regression model
Isoform_Fraction_Disambiguation <- function(obj, sample_name,
                                            quant_name,
                                            fraction_name,
                                            gene_colnames,
                                            transcript_colname,
                                            gene_to_transcript,
                                            count_min = 10,
                                            TPM_min = 1){

  ### Hack to deal with devtools::check() NOTEs
  expected_count <- TPM <- n <- group <- effective_length <- transcript_id <- n <- isoform_cnt <- data <- NULL
  fnest <- logit_fn <- se_logit_fn <- NULL

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
    dplyr::filter(sample == sample_name)

  quant_cols <- colnames(quant)

  if("expected_count" %in% quant_cols){

    quant <- quant %>%
      dplyr::filter(expected_count > count_min)

  }

  if("TPM" %in% quant_cols){

    quant <- quant %>%
      dplyr::filter(TPM > TPM_min)

  }


  # Need to have one row for each transcript ID from a group of
  # transcript IDs, and need to keep track of which reads came from
  # which groups of transcript IDs
  # Also includes a hack to fix issue where a read gets assigned to a transcript
  # isoform on the opposite strand from which the read's RNA originated.
  if(is.null(gene_to_transcript)){

    Fns <- Fns %>%
      dplyr::mutate(group = 1:dplyr::n()) %>%
      tidyr::separate_rows(!!transcript_colname, sep = "\\+") %>%
      dplyr::mutate(transcript_id = !!sym(transcript_colname)) %>%
      dplyr::select(-!!transcript_colname) %>%
      dplyr::inner_join(quant,
                        by = c("transcript_id", "sample"))

  }else{

    Fns <- Fns %>%
      dplyr::mutate(group = 1:dplyr::n()) %>%
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
    dplyr::ungroup() %>%
    dplyr::mutate(fn = !!sym(fraction_of_interest)) %>%
    fit_beta_regression_batch(
      group_cols = gene_colnames,
      keep_cols = gene_colnames,
      deconvolved_feature = "transcript_id"
    ) %>%
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
    dplyr::mutate(n = expected_count,
                  RPK = expected_count/(effective_length/1000),
                  !!logit_fraction_of_interest := pmin(pmax(!!sym(logit_fraction_of_interest), -9), 9),
                  !!fraction_of_interest := inv_logit(!!sym(logit_fraction_of_interest))) %>%
    dplyr::select(-expected_count, -effective_length)

  return(output)


}


#' Deconvolve multi-feature fractions.
#'
#' Combines the output of \code{EstimateFractions} with feature
#' quantification performed by an outside tool (e.g., RSEM, kallisto, salmon, etc.)
#' to infer fraction news for features that reads cannot always be assigned to
#' unambiguously. This is a generalization of \code{EstimateIsoformFractions}
#' which performs this deconvolution for transcript isoforms.
#'
#' \code{DeconvolveFractions} expects as input a "fractions" table with estimates
#' for fraction news of at least one convolved feature set. A convolved feature
#' set is one where some reads cannot be unambiguously assigned to one instance
#' of that feature type. For example, it is often impossible to assign short
#' reads to a single transcript isoform. Thus, something like the "TEC" assignment
#' provided by fastq2EZbakR is an instance of a convolved feature set, as it
#' is an assignment of reads to transcript isoforms with which they are compatible.
#' Another example is assignment to the exonic regions of genes, for fusion genes
#' (where a read may be consistent with both the putative fusion gene as well
#' as one of the fusion components).
#'
#' \code{DeconvolveFractions} deconvolves fraction news
#' by fitting a linear mixing model to the convolved fraction new estimates +
#' feature abundance estimates. In other words, each convolved fraction new (data) is modeled
#' as a weighted average of single feature fraction news (parameters to estimate),
#' with the weights determined by the relative abundances of the features
#' in the convolved set (data). The convolved fraction new is modeled as a beta distribution with mean
#' equal to the weighted feature fraction new average and variance related
#' to the number of reads in the convolved feature set.
#'
#' Features with estimated TPMs less than `TPM_min` (1 by default) or an estimated number of read
#' counts less than `count_min` (10 by default) are removed from convolved feature sets and will
#' not have their fraction news estimated.
#'
#' @param obj An `EZbakRData` object
#' @param feature_type Either `"gene"` (if deconvolving gene-level fraction news,
#' i.e., for resolving fusion gene and component gene fraction news) or
#' `"isoform"` (if deconvolving transcript isoform fraction news).
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
#' @param fraction_name Name of fraction estimate table to use. Should be stored in the
#' `obj$fractions` list under this name. Can also rely on specifying `features` and/or `populations`
#' and having `EZget()` find it.
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
#' @param TPM_min Minimum TPM for a transcript to be kept in analysis.
#' @param count_min Minimum expected_count for a transcript to be kept in analysis.
#' @return An `EZbakRData` object with an additional table under the "fractions"
#' list. Has the same form as the output of `EstimateFractions()`, and will have the
#' feature column "transcript_id".
#' @examples
#'
#' # Load dependencies
#' library(dplyr)
#'
#' # Simulates a single sample worth of data
#' simdata_iso <- SimulateIsoforms(nfeatures = 30)
#'
#' # We have to manually create the metadf in this case
#' metadf <- tibble(sample = 'sampleA',
#'                      tl = 4,
#'                      condition = 'A')
#'
#' ezbdo <- EZbakRData(simdata_iso$cB,
#'                     metadf)
#'
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' ### Hack in the true, simulated isoform levels
#' reads <- simdata_iso$ground_truth %>%
#'   dplyr::select(transcript_id, true_count, true_TPM) %>%
#'   dplyr::mutate(sample = 'sampleA',
#'                 effective_length = 10000) %>%
#'   dplyr::rename(expected_count = true_count,
#'                 TPM = true_TPM)
#'
#' # Name of table needs to have "isoform_quant" in it
#' ezbdo[['readcounts']][['simulated_isoform_quant']] <- reads
#'
#' ### Perform deconvolution
#' ezbdo <- DeconvolveFractions(ezbdo, feature_type = "isoform")
#'
#' @export
DeconvolveFractions <- function(obj,
                                feature_type = c("gene", "isoform"),
                                features = NULL,
                                populations = NULL,
                                fraction_design = NULL,
                                repeatID = NULL,
                                exactMatch = TRUE,
                                fraction_name = NULL,
                                quant_name = NULL,
                                gene_to_transcript = NULL,
                                overwrite = TRUE,
                                TPM_min = 1,
                                count_min = 10){

  feature_type <- match.arg(feature_type)

  if(is.null(quant_name)){

    search_string <- paste0(feature_type, "_quant")

    possible_quant_names <- names(obj$readcounts)

    quant_name <- possible_quant_names[grepl(search_string, possible_quant_names)]

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
                           exactMatch = exactMatch,
                           repeatID = repeatID,
                           returnNameOnly = TRUE)

    if(length(fraction_name) > 1){

      stop("More than one fractions table fits search criterion!
           Either specify additional criterion (i.e., features, populations,
           or fraction_design) or explicitly specify the table name (fraction_name).")
    }

    if(is.null(fraction_name)){
      stop("No tables from your EZbakR analysis match your search criteria!")
    }

  }


  ### Figure out what the gene name is

  fraction <- obj$fractions[[fraction_name]]

  features <- obj[['metadata']][['fractions']][[fraction_name]][['features']]

  fraction <- data.table::as.data.table(fraction)

  # Cute strat: transcript set should have the most unique elements
  lens <- rep(0, times = length(features))
  ambigs <- rep(FALSE, times = length(features))
  for(f in seq_along(features)){

    uniq_fxns <- unique(fraction[[features[f]]])
    lens[f] <- length(uniq_fxns)
    ambigs[f] <- any(grepl("\\+", uniq_fxns))

  }

  if(all(!ambigs)){

    stop("None of the features have multiple assignments to deconvolve!
         There is thus no reason to run DeconvolveFractions() with the
         chosen fractions table.")

  }

  convolved_feature <- features[which(lens == max(lens) & ambigs)]

  other_features <- features[features != convolved_feature]



  ### Estimate fractions
  metadf <- obj$metadf
  # TO-DO; GENERALIZE FOR PULSE-CHASE
  samp_names <- metadf[['sample']][metadf$tl > 0]



  isoform_fit <- purrr::map(.x = samp_names,
                            .f = General_Fraction_Disambiguation,
                            obj = obj,
                            quant_name = quant_name,
                            fraction_name = fraction_name,
                            other_features = other_features,
                            convolved_feature = convolved_feature,
                            gene_to_transcript = gene_to_transcript,
                            feature_type = feature_type,
                            TPM_min = TPM_min,
                            count_min = count_min)

  isoform_fit <- dplyr::bind_rows(isoform_fit)


  ##### DETERMINE OUTPUT NAME AND RETURN

  metadata_list <- obj[['metadata']][['fractions']][[fraction_name]]
  fd <- metadata_list[['fraction_design']]
  pops <- metadata_list[['populations']]

  output_name <- paste0("deconvolved_", feature_type)

  new_feature <- ifelse(
    feature_type == "gene",
    "gene_id",
    "transcript_id"
  )

  # How many identical tables already exist?
  if(overwrite){

    repeatID <- 1

  }else{

    repeatID <- length(EZget(obj,
                             type = 'fractions',
                             features = c(other_features, new_feature),
                             populations = pops,
                             fraction_design = fd,
                             deconvolved = TRUE,
                             returnNameOnly = TRUE,
                             exactMatch = TRUE,
                             alwaysCheck = TRUE)) + 1
  }

  # Get name for output
  output_name <- decide_output(obj, output_name,
                               type = "fractions",
                               features = c(other_features, new_feature),
                               populations = pops,
                               fraction_design = fd,
                               deconvolved = TRUE,
                               overwrite = overwrite)


  obj[['fractions']][[output_name]] <- isoform_fit


  obj[['metadata']][['fractions']][[output_name]] <- list(features = c(other_features, new_feature),
                                                          populations = pops,
                                                          fraction_design = fd,
                                                          deconvolved = TRUE,
                                                          repeatID = repeatID)


  return(obj)


}



# Process EZbakR data so as to fit beta regression model
General_Fraction_Disambiguation <- function(obj, sample_name,
                                            quant_name,
                                            fraction_name,
                                            other_features,
                                            convolved_feature,
                                            feature_type,
                                            gene_to_transcript,
                                            count_min = 10,
                                            TPM_min = 1){

  ### Hack to deal with devtools::check() NOTEs
  expected_count <- TPM <- n <- group <- effective_length <- transcript_id <- n <- isoform_cnt <- data <- NULL
  fnest <- logit_fn <- se_logit_fn <- feature_group_id <- NULL


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
    dplyr::filter(sample == sample_name)

  quant_cols <- colnames(quant)

  if("expected_count" %in% quant_cols){

    quant <- quant %>%
      dplyr::filter(expected_count > count_min)

  }

  if("TPM" %in% quant_cols){

    quant <- quant %>%
      dplyr::filter(TPM > TPM_min)

  }




  # Need to have one row for each transcript ID from a group of
  # transcript IDs, and need to keep track of which reads came from
  # which groups of transcript IDs
  # Also includes a hack to fix issue where a read gets assigned to a transcript
  # isoform on the opposite strand from which the read's RNA originated.

  quant_feature <- ifelse(
    feature_type == "gene",
    "gene_id",
    "transcript_id"
  )


  # Strategy to note which features are ever convolved with other features.
  group_df <- identify_groups(Fns %>%
                                dplyr::select(!!convolved_feature) %>%
                                dplyr::distinct(),
                              quant_feature)



  if(is.null(gene_to_transcript) | feature_type == "gene"){

    Fns <- Fns %>%
      dplyr::mutate(group = 1:dplyr::n()) %>%
      tidyr::separate_rows(!!convolved_feature, sep = "\\+") %>%
      dplyr::mutate(!!quant_feature := !!sym(convolved_feature)) %>%
      dplyr::select(-!!convolved_feature) %>%
      dplyr::inner_join(quant,
                        by = c(quant_feature, "sample"))

  }else{

    tx2g_cols <- colnames(gene_to_transcript)
    gene_colnames <- intersect(other_features,
                               tx2g_cols[tx2g_cols != "transcript_id"])

    if(length(gene_colnames) == 0){
      stop("None of the feature columns match the gene ID column of the
           gene_to_transcript table you provided!")
    }

    Fns <- Fns %>%
      dplyr::mutate(group = 1:dplyr::n()) %>%
      tidyr::separate_rows(!!convolved_feature, sep = "\\+") %>%
      dplyr::mutate(transcript_id = !!sym(convolved_feature)) %>%
      dplyr::select(-!!convolved_feature, -!!gene_colnames) %>%
      dplyr::inner_join(gene_to_transcript,
                        by = c("transcript_id")) %>%
      dplyr::inner_join(quant,
                        by = c("transcript_id", "sample"))

  }


  # What fraction of the reads are expected to have come from each feature?
  Fns <- Fns %>%
    dplyr::inner_join(
      group_df,
      by = quant_feature
    ) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(p = expected_count/effective_length/sum(expected_count/effective_length))



  grouping_cols <- c(other_features, "feature_group_id")


  # Separate features with only one instance that don't need deconvolving
  single_isoforms <- Fns %>%
    dplyr::ungroup() %>%
    dplyr::select(!!grouping_cols, !!quant_feature) %>%
    dplyr::distinct() %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
    dplyr::count() %>%
    dplyr::filter(n == 1)

  Fns_single <- Fns %>%
    dplyr::inner_join(single_isoforms %>% dplyr::select(-n),
                      by = grouping_cols) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", grouping_cols, quant_feature)))) %>%
    dplyr::summarise(!!fraction_of_interest := sum(n*(!!sym(fraction_of_interest)))/sum(n),
                     effective_length = mean(effective_length),
                     expected_count = mean(expected_count),
                     !!logit_fraction_se := mean(!!sym(logit_fraction_se))) %>%
    dplyr::mutate(!!logit_fraction_of_interest := logit(!!sym(fraction_of_interest))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-feature_group_id)


  Fns_multi <- Fns %>%
    dplyr::left_join(single_isoforms %>%
                       dplyr::mutate(isoform_cnt = n) %>%
                       dplyr::select(-n),
                     by = grouping_cols) %>%
    dplyr::filter(is.na(isoform_cnt) & !is.na(expected_count))

  # Run beta regression model to infer isoform-specific fractions

  logit_fraction_of_interest <- paste0("logit_", fraction_of_interest)

  Fns_multi <- Fns_multi %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fn = !!sym(fraction_of_interest)) %>%
    fit_beta_regression_batch(
      group_cols = grouping_cols,
      keep_cols = other_features,
      deconvolved_feature = quant_feature
    ) %>%
    dplyr::mutate(!!logit_fraction_of_interest := logit_fn,
                  !!logit_fraction_se := se_logit_fn,
                  !!fraction_of_interest := inv_logit(logit_fn)) %>%
    dplyr::select(-logit_fn, -se_logit_fn) %>%
    dplyr::inner_join(quant,
                      by = quant_feature) %>%
    dplyr::select(-TPM)



  # Combine single and multi-isoform gene estimates

  output <- Fns_single %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, !!other_features, !!quant_feature, !!fraction_of_interest, !!logit_fraction_of_interest,
                  !!logit_fraction_se, expected_count, effective_length) %>%
    dplyr::bind_rows(Fns_multi %>% dplyr::ungroup()) %>%
    dplyr::select(sample, !!other_features, !!quant_feature, everything()) %>%
    dplyr::mutate(n = expected_count,
                  RPK = expected_count/(effective_length/1000),
                  !!logit_fraction_of_interest := pmin(pmax(!!sym(logit_fraction_of_interest), -9), 9),
                  !!fraction_of_interest := inv_logit(!!sym(logit_fraction_of_interest))) %>%
    dplyr::select(-expected_count, -effective_length)

  return(output)


}


# qf = quant_feature; linear-time Disjoint Set Union (Union-Find) for connected components
identify_groups <- function(df, qf){

  convolved_col <- colnames(df)[1L]
  vals <- unique(as.character(df[[convolved_col]]))
  if(length(vals) == 0L){
    return(dplyr::tibble(
      !!qf := character(0),
      feature_group_id = integer(0)
    ))
  }

  subs_list <- strsplit(vals, "\\+")
  all_nodes <- unique(unlist(subs_list, use.names = FALSE))
  n_nodes <- length(all_nodes)
  parent <- seq_len(n_nodes)

  find_root <- function(x){
    r <- x
    while(parent[r] != r) r <- parent[r]
    while(x != r){
      nxt <- parent[x]
      parent[x] <<- r
      x <- nxt
    }
    r
  }

  for(subs in subs_list){
    if(length(subs) > 1L){
      ids <- match(subs, all_nodes)
      r1 <- find_root(ids[1L])
      for(k in 2L:length(ids)){
        rk <- find_root(ids[k])
        if(rk != r1) parent[rk] <- r1
      }
    }
  }

  roots <- vapply(seq_len(n_nodes), find_root, integer(1))
  comp_id <- match(roots, unique(roots))

  dplyr::tibble(
    !!qf := all_nodes,
    feature_group_id = comp_id
  )

}

