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
#' @param exactMatch If TRUE, then `features` must exactly match the `features`
#' metadata for a given fractions table for it to be used. Means that you cannot
#' specify a subset of features by default. Set this to FALSE if you would like
#' to specify a feature subset.
#' @export
CorrectDropout <- function(obj,
                           grouping_factors = NULL,
                           features = NULL,
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


  corrected <- fractions %>%

    ### CALCULATE DROPOUT:
    dplyr::inner_join(metadf %>% dplyr::filter(tl > 0) %>%
                        dplyr::select(sample, !!grouping_factors),
                      by = "sample") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(rpm = n/(sum(n)/1000000)) %>%
    dplyr::inner_join(nolabel_data,
                      by = c(grouping_factors, features_to_analyze)) %>%
    dplyr::mutate(dropout = rpm / nolabel_rpm) %>%
    dplyr::mutate(sig = sqrt((!!dplyr::sym(logit_se))^2 )) %>%

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
