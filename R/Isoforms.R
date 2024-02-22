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
EstimateIsoformFractions <- function(obj, quantification,
                                     quant_tool = NULL,
                                     cores = 1){

  if(cores > 1){

    # Check if furrr is available
    if (!rlang::is_installed("furrr")){
      stop("You have specified a value of cores > 1. Isoform fraction new estimation
           is only parallelized if the package `furrr` is installed. Either
           install furrr (e.g., with `install.packages('furrr')` or set cores to 1.")

    }

    #
    isoform_fit <- furrr::future_map(.x = samp_names,
                              .f = Isoform_Fraction_Disambiguation,
                              obj = obj,
                              quantification = quantification)

  }else{

    isoform_fit <- purrr::map(.x = samp_names,
                              .f = Isoform_Fraction_Disambiguation,
                              obj = obj,
                              quantification = quantification)

  }





}

# Calculate beta mean
beta_mu <- function(plist, ILparam){

  # Dot product of parameter vector (estimated transcript fns) with probabilities of
  # each isoform being in the
  mus <- sapply(plist, function(v) sum(v * ILparam))

  return(mus)

}

# Calculate parameter prior penalization
calc_prior <- function(ILparam, prior_a, prior_b){

  priors <- dnorm(ILparam, prior_a, prior_b, log = TRUE)

  return(sum(priors))

}


# Beta regression likelihood
beta_lik <- function(param, plist, v, fns, prior_a, prior_b){


  mus <- beta_mu(plist, inv_logit(param))

  alphas <- mus*v
  betas <- v - alphas

  logl <- sum((alphas - 1)*log(fns) + (betas - 1)*log(1 - fns) - lbeta(alphas, betas) )

  logl <- logl + calc_prior(param, prior_a, prior_b)

  return(-logl)

}

Isoform_Fraction_Disambiguation <- function(obj, sample_name,
                                            quantification){

  Fns <- bfo$Fast_Fit$Fn_Estimates %>%
    filter(sample == sample_name) %>%
    dplyr::select(logit_fn, logit_fn_se, XF, nreads) %>%
    dplyr::mutate(fn = inv_logit(logit_fn))

  # RSEM quantification
  setwd(rsem_path)
  rsem_file <- paste0(sample_name, ".isoforms.results")
  rsem <- as_tibble(fread(rsem_file)) %>%
    filter(IsoPct > 0 & TPM >= 1)


  # Add gene name to bakR output
  Fns <- Fns %>%
    mutate(group = 1:n()) %>%
    separate_rows(XF, sep = "\\+") %>%
    mutate(transcript_id = XF) %>%
    dplyr::select(-XF) %>%
    inner_join(rsem,
               by = "transcript_id")

  # Add coverage weights
  Fns <- Fns %>%
    group_by(group) %>%
    mutate(p = expected_count/effective_length/sum(expected_count/effective_length))


  ### Beta regression

  if(debug){
    browser()
  }

  # Filter out genes that only have one isoform
  single_isoforms <- Fns %>%
    ungroup() %>%
    dplyr::select(gene_id, transcript_id) %>%
    dplyr::distinct() %>%
    group_by(gene_id) %>%
    dplyr::count() %>%
    filter(n == 1)

  Fns_single <- Fns %>%
    filter(gene_id %in% single_isoforms$gene_id)

  Fns_multi <- Fns %>%
    filter(!(gene_id %in% single_isoforms$gene_id))


  # Ugly for loop over each gene to estimate fraction news of all its components
  genes <- unique(Fns_multi$gene_id)



  # Effective prior reads
  prior_reads <- 20

  for(g in seq_along(genes)){


    ### Construct list of vectors of proportions of reads coming from each isoform

    Fns_gene <- Fns_multi %>%
      filter(gene_id == genes[g])

    transcripts <- unique(Fns_gene$transcript_id)
    sets <- unique(Fns_gene$group)

    plist <- vector(mode = "list", length = length(sets))

    for(s in seq_along(sets)){

      Fns_set <- Fns_gene %>% filter(group == sets[s])


      tset <- Fns_set$transcript_id

      pvect <- rep(0, times = length(transcripts))
      ps <- Fns_set$p
      count <- 1
      for(t in seq_along(transcripts)){

        if(transcripts[t] %in% tset){
          pvect[t] <- Fns_set$p[Fns_set$transcript_id == transcripts[t]]
          count <- count + 1
        }else{
          pvect[t] <- 0
        }


      }
      plist[[s]] <- pvect



    }

    nt <- length(transcripts)

    low_ps <- rep(-7, times = nt)
    high_ps <- rep(7, times = nt)

    Fns_group <- Fns_gene %>%
      dplyr::select(fn, nreads, group) %>%
      dplyr::distinct()

    ### Estimate prior
    prior <- Fns_group %>%
      ungroup() %>%
      summarise(fn_gene = sum(fn*nreads)/sum(nreads)) %>%
      unlist() %>%
      unname()

    # prior_a <- prior*(prior_reads)
    # prior_b <- prior_reads - prior_a

    # logit_normal prior instead
    prior_a <- logit(prior)
    prior_b <- 1


    ### Fit model

    fit <- stats::optim(par=rep(0, times = nt),
                        fn = beta_lik,
                        plist = plist, v = Fns_group$nreads, fns = Fns_group$fn,
                        prior_a = prior_a, prior_b = prior_b,
                        method = "L-BFGS-B",
                        lower = low_ps, upper = high_ps)

    fns <- fit$par

    if(g == 1){

      final_df <- data.frame(transcript_id = transcripts,
                             logit_fn = fit$par) %>%
        mutate(fn = inv_logit(logit_fn))

    }else{

      estimate_df <- tibble(transcript_id = transcripts,
                            logit_fn = fit$par) %>%
        mutate(fn = inv_logit(logit_fn))


      final_df <- rbind(final_df, estimate_df)


    }



  }


  return(list(multi_isoforms = final_df,
              single_isoforms = Fns_single))


}
