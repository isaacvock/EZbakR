#' Correct for experimental/bioinformatic dropout of labeled RNA.
#'
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
#' @export
CorrectDropout <- function(obj,
                           grouping_factors = NULL,
                           features = NULL,
                           populations = NULL,
                           fraction_design = NULL,
                           repeatID = NULL,
                           exactMatch = TRUE){


  ### Hack to deal with devtools::check() NOTEs

  tl <- n <- nolabel_rpm <- rpm <- pdo <- fit <- global_fraction <- corrected_gf <- corrected_fraction <- corrected_n <- NULL
  nolabel_n <- nolabel_reps <- dropout <- sig <- density <- NULL


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
                                 exactMatch = exactMatch) %>%

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



# VisualizeDropout <- function(obj,
#                              grouping_factors = NULL,
#                              features = NULL,
#                              exactMatch = NULL)


#' Make plots to visually assess dropout trends
#'
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
        density = get_density(
          x = !!dplyr::sym(fraction_of_interest),
          y = dropout,
          n = 200
        )
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(x = !!dplyr::sym(fraction_of_interest),
            y = dropout,
            color = density)
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
                              exactMatch = TRUE){

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
    dplyr::filter(n > 25 & dropout < 5) # Get rid of outliers


  return(dropout_df)

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
