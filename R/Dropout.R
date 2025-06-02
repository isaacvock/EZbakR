#' Correct for experimental/bioinformatic dropout of labeled RNA.
#'
#' Uses the strategy described [here](https://simonlabcode.github.io/bakR/articles/Dropout.html), and similar to that originally presented
#' in [Berg et al. 2024](https://academic.oup.com/nar/article/52/7/e35/7612100).
#'
#' Dropout is the disproportionate loss of labeled RNA/reads from said RNA
#' described independently [here](https://academic.oup.com/nar/article/52/7/e35/7612100)
#' and [here](https://www.biorxiv.org/content/10.1101/2023.05.24.542133v1). It can originate from a combination of
#' bioinformatic (loss of high mutation content reads due to alignment problems),
#' technical (loss of labeled RNA during RNA extraction), and biological (transcriptional
#' shutoff in rare cases caused by metabolic label toxicity) sources.
#' `CorrectDropout()` compares label-fed and label-free controls from the same
#' experimental conditions to estimate and correct for this dropout. It assumes
#' that there is a single number (referred to as the dropout rate, or pdo) which
#' describes the rate at which labeled RNA is lost (relative to unlabeled RNA).
#' pdo ranges from 0 (no dropout) to 1 (complete loss of all labeled RNA), and
#' is thus interpreted as the percentage of labeled RNA/reads from labeled RNA
#' disproportionately lost, relative to the equivalent unlabeled species.
#'
#' @param obj An EZbakRFractions object, which is an EZbakRData object on which
#' you have run `EstimateFractions()`.
#' @param grouping_factors Which sample-detail columns in the metadf should be used
#' to group -s4U samples by for calculating the average -s4U RPM? The default value of
#' `NULL` will cause all sample-detail columns to be used.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of `NULL`
#' will expect there to be only one fractions table in the EZbakRFractions object.
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
#' @param exactMatch If TRUE, then `features` must exactly match the `features`
#' metadata for a given fractions table for it to be used. Means that you cannot
#' specify a subset of features by default. Set this to FALSE if you would like
#' to specify a feature subset.
#' @param read_cutoff Minimum number of reads for a feature to be used to fit
#' the dropout model.
#' @param dropout_cutoff Maximum ratio of -s4U:+s4U RPMs for a feature to be
#' used to fit the dropout model (i.e., simple outlier filtering cutoff).
#' @return An `EZbakRData` object with the specified "fractions" table replaced
#' with a dropout corrected table.
#' @importFrom magrittr %>%
#' @export
CorrectDropout <- function(obj,
                           grouping_factors = NULL,
                           features = NULL,
                           populations = NULL,
                           fraction_design = NULL,
                           repeatID = NULL,
                           exactMatch = TRUE,
                           read_cutoff = 25,
                           dropout_cutoff = 5){


  ### Hack to deal with devtools::check() NOTEs

  tl <- n <- nolabel_rpm <- rpm <- pdo <- fit <- global_fraction <- corrected_gf <- corrected_fraction <- corrected_n <- NULL
  nolabel_n <- nolabel_reps <- dropout <- sig <- NULL


  ##### GENERAL STEPS:
  # 1) Get -s4U RPMs
  # 2) Get +s4U RPMs
  # 3) Calculate dropout (+s4U RPM)/(-s4U RPM)
  # 4) Fit dropout vs. estimated fraction new trend
  # 5) Correct fraction news and read counts accordingly.
  # 6) Return corrected EZbakRFractions object


  ### Figure out which fraction new estimates to use

  # Function in Helpers.R
  fractions_name <- EZget(obj, type = "fractions",
                          features = features,
                          exactMatch = exactMatch,
                          populations = populations,
                          fraction_design = fraction_design,
                          returnNameOnly = TRUE)

  # Get fractions
  fractions <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[['metadata']][['fractions']][[fractions_name]][['features']]


  ### Figure out column names to operate on

  fraction_cols <- colnames(fractions)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  # Other columns I will have to reference
  logit_fraction <- paste0("logit_", fraction_of_interest)
  logit_se <- paste0("se_", logit_fraction)


  ### Quantify and correct for dropout

  corrected <- calculate_dropout(obj = obj,
                                 grouping_factors = grouping_factors,
                                 features = features,
                                 populations = populations,
                                 fraction_design = fraction_design,
                                 repeatID = repeatID,
                                 exactMatch = exactMatch,
                                 read_cutoff = read_cutoff,
                                 dropout_cutoff = dropout_cutoff) %>%

    ### ESTIMATE DROPOUT PARAMETERS:
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      fit = I(list(stats::optim(
        par = c(0.2, 2),
        dropout_likelihood,
        theta = !!dplyr::sym(fraction_of_interest),
        dropout = dropout,
        sig = sig,
        method = "L-BFGS-B",
        upper = c(0.99, 100),
        lower = c(0.01, 1) # Scale factor can't be < 1 because -s4U molecule count > +s4U molecule count
      )))
    ) %>%
    dplyr::mutate(
      pdo := purrr::map_dbl(fit, ~ .x$par[1]),
      scale := purrr::map_dbl(fit, ~ .x$par[2]),
    ) %>%
    dplyr::select(-fit) %>%

    ### CORRECT DROPOUT:
    dplyr::mutate(corrected_fraction =
                    (!!dplyr::sym(fraction_of_interest))/((1 - pdo) + !!dplyr::sym(fraction_of_interest) * pdo)) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(global_fraction = sum(!!dplyr::sym(fraction_of_interest)*n)/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(corrected_gf = global_fraction/((1 - pdo) + global_fraction * pdo)) %>%
    dplyr::mutate(corrected_n = ceiling(n * (corrected_gf*(1 - pdo) + (1 - corrected_gf) /
                    (corrected_fraction*(1 - pdo) + (1 - corrected_fraction))))
    ) %>%
    dplyr::mutate(!!fraction_of_interest := corrected_fraction,
                  n = corrected_n,
                  !!logit_fraction := logit(corrected_fraction)) %>%
    dplyr::select(sample, !!features_to_analyze, !!fraction_of_interest,
                  !!logit_fraction, !!logit_se, n, pdo, scale)


  ### Get dropout parameters
  dropout_params <- corrected %>%
    dplyr::select(sample, pdo, scale) %>%
    dplyr::distinct()


  message(paste0(c("Estimated rates of dropout are:",
                   utils::capture.output(as.data.frame(dropout_params[,c("sample", "pdo")]))),
                 collapse = "\n"))

  ### Add back -s4U data
  corrected <- corrected %>%
    dplyr::bind_rows(fractions %>%
                       dplyr::filter(!!dplyr::sym(fraction_of_interest) == 0))


  ### Overwrite uncorrected data
  obj[["fractions"]][[fractions_name]] <- corrected %>%
    dplyr::select(-pdo, -scale)

  obj[["dropout_params"]][[fractions_name]] <- dropout_params


  return(obj)

}



#' Make plots to visually assess dropout trends
#'
#' Plots a measure of dropout (the ratio of -label to +label RPM) as a function
#' of feature fraction new, with the model fit depicted. Use this function to
#' qualitatively assess model fit and whether the modeling assumptions are met.
#'
#' @param obj An EZbakRFractions object, which is an EZbakRData object on which
#' you have run `EstimateFractions()`.
#' @param grouping_factors Which sample-detail columns in the metadf should be used
#' to group -s4U samples by for calculating the average -s4U RPM? The default value of
#' `NULL` will cause all sample-detail columns to be used.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of `NULL`
#' will expect there to be only one fractions table in the EZbakRFractions object.
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
#' @param exactMatch If TRUE, then `features` must exactly match the `features`
#' metadata for a given fractions table for it to be used. Means that you cannot
#' specify a subset of features by default. Set this to FALSE if you would like
#' to specify a feature subset.
#' @param n_min Minimum raw number of reads to make it to plot
#' @importFrom magrittr %>%
#' @return A list of `ggplot2` objects, one for each +label sample.
#' @export
VisualizeDropout <- function(obj,
                             grouping_factors = NULL,
                             features = NULL,
                             populations = NULL,
                             fraction_design = NULL,
                             repeatID = NULL,
                             exactMatch = TRUE,
                             n_min = 50){


  ### Hack to deal with devtools::check() NOTEs

  tl <- n <- nolabel_rpm <- rpm <- pdo <- fit <- global_fraction <- corrected_gf <- corrected_fraction <- corrected_n <- NULL
  nolabel_n <- nolabel_reps <- dropout <- sig <- point_density <- NULL


  ##### GENERAL STEPS:
  # 1) Get -s4U RPMs
  # 2) Get +s4U RPMs
  # 3) Calculate dropout (+s4U RPM)/(-s4U RPM)
  # 4) Fit dropout vs. estimated fraction new trend
  # 5) Correct fraction news and read counts accordingly.
  # 6) Return corrected EZbakRFractions object


  ### Figure out which fraction new estimates to use

  # Function in Helpers.R
  fractions_name <- EZget(obj, type = "fractions",
                          features = features,
                          exactMatch = exactMatch,
                          populations = populations,
                          fraction_design = fraction_design,
                          returnNameOnly = TRUE)

  # Get fractions
  fractions <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[['metadata']][['fractions']][[fractions_name]][['features']]


  ### Figure out column names to operate on

  fraction_cols <- colnames(fractions)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  # Other columns I will have to reference
  logit_fraction <- paste0("logit_", fraction_of_interest)
  logit_se <- paste0("se_", logit_fraction)


  ### Calculate dropout

  dropout_df <- calculate_dropout(obj = obj,
                                  grouping_factors = grouping_factors,
                                  features = features,
                                  populations = populations,
                                  fraction_design = fraction_design,
                                  repeatID = repeatID,
                                  exactMatch = exactMatch)

  ### Make dropout plots

  samps <- unique(dropout_df[['sample']])

  plot_list <- vector(mode = "list", length = length(samps))

  for(s in seq_along(samps)){

    message(paste0("Making plot for sample ", samps[s], "..."))

    plot_list[[s]] <- dropout_df %>%
      dplyr::filter(sample == samps[s] &
                      n > n_min) %>%
      dplyr::mutate(
        point_density = get_density(
          x = !!dplyr::sym(fraction_of_interest),
          y = dropout,
          n = 200
        )
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(x = !!dplyr::sym(fraction_of_interest),
            y = dropout,
            color = point_density)
      ) +
      ggplot2::geom_point() +
      ggplot2::theme_classic() +
      ggplot2::scale_color_viridis_c() +
      ggplot2::xlab("fraction labeled") +
      ggplot2::ylab("Dropout") +
      ggplot2::geom_hline(
        yintercept = 1,
        color = 'darkred',
        linewidth = 0.75,
        linetype = 'dotted'
      )

  }

  names(plot_list) <- samps


  return(plot_list)

}


# Calculate dropout for each feature in each +label sample
calculate_dropout <- function(obj,
                              grouping_factors = NULL,
                              features = NULL,
                              populations = NULL,
                              fraction_design = NULL,
                              repeatID = NULL,
                              exactMatch = TRUE,
                              read_cutoff = 25,
                              dropout_cutoff = 5){

  ### Hack to deal with devtools::check() NOTEs

  tl <- n <- nolabel_rpm <- rpm <- pdo <- fit <- global_fraction <- corrected_gf <- corrected_fraction <- corrected_n <- NULL
  nolabel_n <- nolabel_reps <- dropout <- sig <- NULL


  ##### GENERAL STEPS:
  # 1) Get -s4U RPMs
  # 2) Get +s4U RPMs
  # 3) Calculate dropout (+s4U RPM)/(-s4U RPM)
  # 4) Fit dropout vs. estimated fraction new trend
  # 5) Correct fraction news and read counts accordingly.
  # 6) Return corrected EZbakRFractions object


  ### Figure out which fraction new estimates to use

  # Function in Helpers.R
  fractions_name <- EZget(obj, type = "fractions",
                          features = features,
                          exactMatch = exactMatch,
                          populations = populations,
                          fraction_design = fraction_design,
                          returnNameOnly = TRUE)

  # Get fractions
  fractions <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[['metadata']][['fractions']][[fractions_name]][['features']]


  ### Figure out column names to operate on

  fraction_cols <- colnames(fractions)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  # Other columns I will have to reference
  logit_fraction <- paste0("logit_", fraction_of_interest)
  logit_se <- paste0("se_", logit_fraction)


  ### Which columns should -s4U samples be grouped by?

  metadf <- obj$metadf


  if(is.null(grouping_factors)){

    grouping_factors <- colnames(metadf)[!grepl("^tl", colnames(metadf)) &
                                           (colnames(metadf) != "sample") &
                                           !grepl("^tpulse", colnames(metadf)) &
                                           !grepl("^tchase", colnames(metadf))]


  }


  ### Quantify and correct for dropout


  # Necessary generalizations:
  # 1) Metadf column used (e.g., pulse-chase)
  nolabel_data <- fractions %>%
    dplyr::inner_join(metadf %>% dplyr::filter(tl == 0),
                      by = "sample") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(nolabel_rpm = n/(sum(n)/1000000),
                  nolabel_n = n) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(grouping_factors, features_to_analyze)))) %>%
    dplyr::summarise(nolabel_rpm = mean(nolabel_rpm),
                     nolabel_n = mean(nolabel_n),
                     nolabel_reps = dplyr::n()) %>%
    dplyr::select(!!grouping_factors, !!features_to_analyze,
                  nolabel_rpm, nolabel_n, nolabel_reps)

  if(nrow(nolabel_data) == 0){

    stop("No -s4U samples found!!")

  }

  ### CALCULATE DROPOUT:
  dropout_df <- fractions %>%
    dplyr::inner_join(metadf %>% dplyr::filter(tl > 0) %>%
                        dplyr::select(sample, !!grouping_factors),
                      by = "sample") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(rpm = n/(sum(n)/1000000)) %>%
    dplyr::inner_join(nolabel_data,
                      by = c(grouping_factors, features_to_analyze)) %>%
    dplyr::mutate(dropout = rpm / nolabel_rpm) %>%
    dplyr::mutate(sig = sqrt((!!dplyr::sym(logit_se))^2 )) %>%
    dplyr::filter(n > read_cutoff & dropout < dropout_cutoff) # Get rid of outliers


  return(dropout_df)

}


#' Normalize for experimental/bioinformatic dropout of labeled RNA.
#'
#' Uses the strategy described [here](https://simonlabcode.github.io/bakR/articles/Dropout.html), and similar to that originally presented
#' in [Berg et al. 2024](https://academic.oup.com/nar/article/52/7/e35/7612100),
#' to normalize for dropout. Normalizing for dropout means identifying a reference
#' sample with low dropout, and estimating dropout in each sample relative to
#' that sample.
#'
#' `NormalizeForDropout()` has a number of unique advantages relative to
#' `CorrectDropout()`:
#'
#'  - `NormalizeForDropout()` doesn't require -label control data.
#'  - `NormalizeForDropout()` compares an internally normalized quantity
#'  (fraction new) across samples, which has some advantages over the
#'  absolute dropout estimates derived from comparisons of normalized read
#'  counts in `CorrectDropout()`.
#'  - `NormalizeForDropout()` may be used to normalize half-life estimates
#'  across very different biological contexts (e.g., different cell types).
#'
#'  There are also some caveats to be aware of when using `NormalizeForDropout()`:
#'
#'  - Be careful using `NormalizeForDropout()` when you have multiple different
#'  label times. Dropout normalization requires each sample be compared to a reference
#'  sample with the same label time. Thus, normalization will be performed
#'  separately for groups of samples with different label times. If the extent
#'  of dropout in the references with different label times is different, there
#'  will still be unaccounted for dropout biases between some of the samples.
#'  - `NormalizeForDropout()` effectively assumes that there are no true global
#'  differences in turnover kinetics of RNA. If such differences actually exist
#'  (e.g., half-lives in one context are on average truly lower than those in
#'  another), `NormalizeForDropout()` risks normalizing away these real
#'  differences. This is similar to how statistical normalization strategies
#'  implemented in differential expression analysis software like DESeq2 assumes
#'  that there are no global differences in RNA levels.
#'
#'  By default, all samples with same label time are normalized with respect
#'  to a reference sample chosen from among them. If you want to further separate
#'  the groups of samples that are normalized together, specify the columns of
#'  your metadf by which you want to additionally group factors in the `grouping_factors`
#'  parameter. This behavior can be changed by setting `normalize_across_tls` to
#'  `TRUE`, which will
#'
#' @param obj An EZbakRFractions object, which is an EZbakRData object on which
#' you have run `EstimateFractions()`.
#' @param normalize_across_tls If TRUE, samples from different label times will
#' be normalized by finding a max inferred degradation rate constant (kdeg)
#' sample and using that as a reference. Degradation kinetics with this max
#' will be assumed to infer reference fraction news at different label times
#' @param grouping_factors Which sample-detail columns in the metadf should be used
#' to group -s4U samples by for calculating the average -s4U RPM? The default value of
#' `NULL` will cause no sample-detail columns to be used.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of `NULL`
#' will expect there to be only one fractions table in the EZbakRFractions object.
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
#' @param exactMatch If TRUE, then `features` must exactly match the `features`
#' metadata for a given fractions table for it to be used. Means that you cannot
#' specify a subset of features by default. Set this to FALSE if you would like
#' to specify a feature subset.
#' @param read_cutoff Minimum number of reads for a feature to be used to fit
#' dropout model.
#' @return An `EZbakRData` object with the specified "fractions" table replaced
#' with a dropout corrected table.
#' @importFrom magrittr %>%
#' @export
NormalizeForDropout <- function(obj,
                                 normalize_across_tls = FALSE,
                                 grouping_factors = NULL,
                                 features = NULL,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 repeatID = NULL,
                                 exactMatch = TRUE,
                                read_cutoff = 25){


  ### Hack to deal with devtools::check() NOTEs

  tl <- n <- nolabel_rpm <- rpm <- pdo <- fit <- global_fraction <- corrected_gf <- corrected_fraction <- corrected_n <- NULL
  nolabel_n <- nolabel_reps <- dropout <- sig <- `.` <- NULL
  med_fn <- normalization_reference <- reference_samp <- NULL
  nsamps_in_group <- num <- den <- NULL

  ##### GENERAL STEPS:
  # 1) Get -s4U RPMs
  # 2) Get +s4U RPMs
  # 3) Calculate dropout (+s4U RPM)/(-s4U RPM)
  # 4) Fit dropout vs. estimated fraction new trend
  # 5) Correct fraction news and read counts accordingly.
  # 6) Return corrected EZbakRFractions object


  ### Figure out which fraction new estimates to use

  # Function in Helpers.R
  fractions_name <- EZget(obj, type = "fractions",
                          features = features,
                          exactMatch = exactMatch,
                          populations = populations,
                          fraction_design = fraction_design,
                          returnNameOnly = TRUE)

  # Get fractions
  fractions <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[['metadata']][['fractions']][[fractions_name]][['features']]


  ### Figure out column names to operate on

  fraction_cols <- colnames(fractions)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  # Other columns I will have to reference
  logit_fraction <- paste0("logit_", fraction_of_interest)
  logit_se <- paste0("se_", logit_fraction)


  ### Which columns should -s4U samples be grouped by?

  metadf <- obj$metadf

  # Which column is label time column?
    # Could generalize to pulse-chases at some point
  tl_col <- "tl"


  if(!normalize_across_tls){
    grouping_factors <- unique(c(grouping_factors, tl_col))
  }



  ### Figure out which sample likely has lowest dropout in each group

  # Reference sample in each group is that which has the lowest
  # dropout and thus that other samples in that group will be
  # normalized with respect to
  reference_sample_info <- fractions %>%
    dplyr::inner_join(
      metadf,
      by = "sample"
    ) %>%
    dplyr::filter(tl > 0 & n > read_cutoff) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", grouping_factors)))) %>%
    dplyr::summarise(
      med_fn = ifelse(
        normalize_across_tls,
        stats::median(-log(1 - !!dplyr::sym(fraction_of_interest))/tl ), # Inferred kdeg
        stats::median(!!dplyr::sym(logit_fraction))
        )
    ) %>%
    dplyr::ungroup() %>%
    {if (length(grouping_factors) > 0) dplyr::group_by(., dplyr::across(dplyr::all_of(grouping_factors))) else .} %>%
    dplyr::mutate(
      normalization_reference = dplyr::case_when(
        med_fn == max(med_fn) ~ TRUE,
        .default = FALSE
      ),
      reference_samp = sample[med_fn == max(med_fn)],
      nsamps_in_group = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      sample, !!grouping_factors, normalization_reference, reference_samp, nsamps_in_group
    )


  ### Quantify dropout with respect to reference in each group


  fxn_wide <- fractions %>%
    dplyr::filter(n > read_cutoff & !!dplyr::sym(logit_fraction) > -Inf) %>%
    dplyr::select(sample, !!features_to_analyze, !!fraction_of_interest) %>%
    tidyr::pivot_wider(
      names_from = "sample",
      values_from = fraction_of_interest
    )

  pdos <- metadf %>%
    dplyr::filter(tl > 0) %>%
    dplyr::mutate(
      pdo = 0
    )

  for(i in seq_along(pdos$sample)){

    samp <- pdos$sample[i]

    if(reference_sample_info$normalization_reference[reference_sample_info$sample == samp]) next

    reference <- reference_sample_info$reference_samp[reference_sample_info$sample == samp]


    fxn_subset <- fxn_wide %>%
      dplyr::select(
        !!features_to_analyze,
        !!samp,
        !!reference
      )

    if(normalize_across_tls){

      ref_tl <- metadf[[tl_col]][metadf$sample == reference]
      samp_tl <- metadf[[tl_col]][metadf$sample == samp]

      # Infer reference fraction new assuming exponential decay kinetics
      fxn_subset <- fxn_subset %>%
        dplyr::mutate(
          # kdeg = -log(1 - fn)/tl
          # fn = 1 - exp(-kdeg*tl)
          !!reference := 1 - exp((log(1 - !!dplyr::sym(reference))/ref_tl) * samp_tl)
        )

    }

    pdos$pdo[i] <- fxn_subset %>%
      stats::na.omit() %>%
      dplyr::summarise(
        pdo = stats::optim(
          par = c(-2),
          fn = do_norm_ll,
          reffn = !!dplyr::sym(reference),
          fndo = !!dplyr::sym(samp),
          sig = 0.2, # Should be some function of the uncertainties
          method = "L-BFGS-B",
          upper = 6,
          lower = -6
        )$par[1]
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(pdo) %>%
      unlist() %>%
      unname() %>%
      inv_logit()


  }

  ### Normalize for dropout
  fractions_normalized <- fractions %>%
    dplyr::left_join(
      pdos %>%
        dplyr::select(sample, pdo),
      by = "sample"
    ) %>%
    dplyr::mutate(
      pdo = ifelse(is.na(pdo), 0, pdo)
    ) %>%
    dplyr::mutate(
      !!fraction_of_interest := (!!dplyr::sym(fraction_of_interest))/((1 - pdo) + !!dplyr::sym(fraction_of_interest) * pdo),
    ) %>%
    dplyr::mutate(
      !!logit_fraction := logit(!!dplyr::sym(fraction_of_interest))
    ) %>%
    dplyr::group_by(
      sample
    ) %>%
    dplyr::mutate(
      num = mean(!!dplyr::sym(fraction_of_interest))*(1-pdo) + (1-mean(!!dplyr::sym(fraction_of_interest))),
      den = !!dplyr::sym(fraction_of_interest)*(1-pdo) + (1 - !!dplyr::sym(fraction_of_interest)),
      n = ceiling(n*(num/den))
    ) %>%
    dplyr::select(
      -num, -den, -pdo
    )

  obj[["fractions"]][[fractions_name]] <- fractions_normalized
  obj[['metadata']][['fractions']][[fractions_name]][['dropout']] <- pdos

  return(obj)


}


do_norm_ll <- function(param,
                       reffn,
                       fndo,
                       sig){

  pdo <- inv_logit(param[1])

  Efn <- fndo / ((1 - pdo) + (fndo*pdo))

  ll <- stats::dnorm(logit(reffn),
                     logit(Efn),
                     sig,
                     log = TRUE)

  return(-sum(ll))

}



# Dropout likelihood function
# log(dropout) ~ Normal(f(pdo, theta, scale), sig)
dropout_likelihood <- function(param, dropout, theta, sig){

  pdo <- param[1]
  scale <- param[2]

  Edropout <- log((-(scale*pdo)*theta)/((1-pdo) + theta*pdo) + scale)

  ll <- stats::dnorm(dropout,
                     Edropout,
                     sig,
                     log = TRUE)

  return(-sum(ll))


}


# Estimating standard deviation of log(reads)
read_sig <- function(reads, reps){

  return(sqrt(1/reads)/sqrt(reps))

}


# Likelihood function for dropout normalization
do_norm_ll <- function(param,
                       reffn,
                       fndo,
                       sig){

  pdo <- inv_logit(param[1])

  Efn <- fndo / ((1 - pdo) + (fndo*pdo))

  ll <- stats::dnorm(logit(reffn),
                     logit(Efn),
                     sig,
                     log = TRUE)

  return(-sum(ll))

}


# Make sure dropout normalization is possible
check_donorm_validity <- function(
    fractions,
    metadf,
    features_to_analyze,
    logit_fraction,
    logit_se,
    grouping_factors
){

  ### Checklist:
  # 1) Need 2 or more replicates of each grouping factor




}
