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
#' @param quant_name Name of quantification tool appended to table name of interest.
#' This is only relevant if you are providing isoform-specific estimates, in which
#' case the isoform quantification tool's name may need to be provided in order
#' for EZbakR to uniquely identify the table of interest. Even in that case though,
#' this should only have to be non-null in the case where you have performed isoform-specific
#' fraction estimation with more than one quantification tool's output.
CorrectDropout <- function(obj,
                           grouping_factors = NULL,
                           features = NULL,
                           quant_name = NULL){

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

    grouping_factors <- colnames(metadf)[!grepl("tl", colnames(metadf)) &
                                           (colnames(metadf) != "sample")]


  }


  ### Quantify and correct for dropout


  # Necessary generalizations:
  # 1) Metadf column used (e.g., pulse-chase)
  nolabel_data <- fractions %>%
    dplyr::inner_join(metadf %>% dplyr::filter(tl == 0),
                      by = "sample") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(nolabel_rpm = n/(sum(n)/1000000)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(grouping_factors, features_to_analyze)))) %>%
    dplyr::summarise(nolabel_rpm = mean(nolabel_rpm)) %>%
    dplyr::select(!!grouping_factors, !!features_to_analyze, nolabel_rpm)


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

    ### ESTIMATE DROPOUT PARAMETERS:
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      fit = I(list(summary(stats::nls(log(dropout) ~ log((-(scale*pdo)*fraction_highTC)/((1-pdo) + fraction_highTC*pdo) + scale),
                                      start = list(scale = 1, pdo = 0.2)))))
    ) %>%
    dplyr::mutate(
      pdo := purrr::map_dbl(fit, ~ .x$coefficients[2,1]),
      scale := purrr::map_dbl(fit, ~ .x$coefficients[1,1]),
    ) %>%
    dplyr::mutate(
      pdo = ifelse(pdo < 0, 0, pdo)
    ) %>%
    dplyr::select(-fit) %>%

    ### CORRECT DROPOUT:
    dplyr::mutate(corrected_fraction =
                    (!!dplyr::sym(fraction_of_interest))/((1 - pdo) + !!dplyr::sym(fraction_of_interest) * pdo)) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(global_fraction = sum(!!dplyr::sym(fraction_of_interest)*n)/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(corrected_gf = global_fraction/((1 - pdo) + global_fraction * pdo)) %>%
    dplyr::mutate(corrected_n = n * (corrected_fraction*(1 - pdo) + (1 - corrected_fraction)) /
                    (corrected_gf*(1 - pdo) + (1 - corrected_gf)),
                  correction_factor = corrected_fraction / !!dplyr::sym(fraction_of_interest)
    ) %>%
    dplyr::mutate(!!fraction_of_interest := corrected_fraction,
                  n = corrected_n,
                  !!logit_fraction := logit(corrected_fraction),
                  !!logit_se := correction_factor*!!dplyr::sym(logit_se)) %>%
    dplyr::select(sample, !!features_to_analyze, !!fraction_of_interest,
                  !!logit_fraction, !!logit_se, n, pdo, scale)


  ### Get dropout parameters
  dropout_params <- corrected %>%
    dplyr::select(sample, pdo, scale) %>%
    dplyr::distinct()


  if(any(dropout_params$pdo > 1)){
    stop("Dropout was estimated to be 100%, in one of your samples, which is likely an estimation error.
              Is your data of unusually low coverage or metabolic label incorporation rate? This can lead to estimation problems.")
  }

  message(paste0(c("Estimated rates of dropout are:",
                   utils::capture.output(as.data.frame(dropout_params[,c("sample", "pdo")]))),
                 collapse = "\n"))

  ### Add back -s4U data
  corrected <- corrected %>%
    dplyr::bind_rows(fractions %>%
                       filter(fraction_highTC == 0))


  ### Overwrite uncorrected data
  obj[["fractions"]][[fractions_name]] <- corrected %>%
    dplyr::select(-pdo, -scale)

  obj[["dropout_params"]][[fractions_name]] <- dropout_params


  return(obj)

}
