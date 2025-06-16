#' Make a VolcanoPlot from EZbakR comparison
#'
#' Make a plot of effect size (x-axis) vs. multiple-test adjusted p-value (y-axis),
#' coloring points by position relative to user-defined decision cutoffs.
#'
#' `EZVolcanoPlot()` accepts as input the output of `CompareParameters()`, i.e.,
#' an `EZbakRData` object with at least one "comparisons" table. It will plot
#' the "difference" column in this table versus -log10 of the "padj" column.
#' In the simplest case, "difference" represents a log-fold change in a kinetic
#' parameter (e.g., kdeg) estimate. More complicated linear model fits and
#' comparisons can yield different parameter estimates.
#'
#' NOTE: some outputs of `CompareParameters()` are not meant for visualization
#' via a volcano plot. For example, when fitting certain interaction models,
#' some of the parameter estimates may represent average log(kinetic paramter)
#' in one condition. See discussion of one example of this [here](https://github.com/simonlabcode/bakR/issues/7#issuecomment-2371431127).
#'
#' EZbakR estimates kinetic parameters in `EstimateKinetics()` and `EZDynamics()`
#' on a log-scale. By default, since log2-fold changes are a bit easier to interpret
#' and more common for these kind of visualizations, `EZVolcanoPlot()` multiplies
#' the x-axis value by log2(exp(1)), which is the factor required to convert from
#' a log to a log2 scale. You can turn this off by setting `plotlog2` to `FALSE`.
#'
#' @param obj An object of class `EZbakRCompare`, which is an `EZbakRData` object
#' on which you have run `CompareParameters`
#' @param parameter Name of parameter whose comparison you want to plot.
#' @param design_factor Name of factor from `metadf` whose parameter estimates at
#' different factor values you would like to compare.
#' @param reference Name of reference `condition` factor level value.
#' @param experimental Name of `condition` factor level value to compare to reference.
#' @param param_name If you want to assess the significance of a single parameter,
#' rather than the comparison of two parameters, specify that one parameter's name
#' here.
#' @param param_function NOT YET IMPLEMENTED. Will allow you to specify more complicated
#' functions of parameters when hypotheses you need to test are combinations of parameters
#' rather than individual parameters or simple differences in two parameters.
#' @param features Character vector of feature names for which comparisons were made.
#' @param condition Defunct parameter that has been replaced with `design_factor`. If provided
#' gets passed to `design_factor` if `design_factor` is not already specified.
#' @param normalize_by_median Whether or not the median was subtracted from the estimated
#' parameter differences.
#' @param repeatID If multiple `kinetics` or `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` has to exactly match
#' those for a given comparisons table for that table to be used. Means that you can't
#' specify a subset of features by default, since this is TRUE
#' by default.
#' @param plotlog2 If TRUE, assume that log(parameter) difference is passed in and that
#' you want to plot log2(parameter) difference. TO-DO: probably best to change this
#' to a more general scale parameter by which the parameter is multiplied. Default
#' would be log2(exp(1)) to convert log() to log2().
#' @param FDR_cutoff False discovery cutoff by which to color points.
#' @param difference_cutoff Minimum absolute difference cutoff by which to color points.
#' @param size Size of points, passed to `geom_point()` size parameter. If not specified,
#' a point size is automatically chosen.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @examples
#'
#' # Simulate data to analyze
#' simdata <- EZSimulate(30)
#'
#' # Create EZbakR input
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#'
#' # Estimate Fractions
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' # Estimate Kinetics
#' ezbdo <- EstimateKinetics(ezbdo)
#'
#' # Average estimates across replicate
#' ezbdo <- AverageAndRegularize(ezbdo)
#'
#' # Compare parameters across conditions
#' ezbdo <- CompareParameters(
#' ezbdo,
#' design_factor = "treatment",
#' reference = "treatment1",
#' experimental = "treatment2"
#' )
#'
#' # Make volcano plot (ggplot object that you can save and add/modify layers)
#' EZVolcanoPlot(ezbdo)
#'
#' @return A `ggplot2` object. X-axis = log2(estimate of interest (e.g., fold-change
#' in degradation rate constant); Y-axis = -log10(multiple test adjusted p-value);
#' points colored by location relative to FDR and effect size cutoffs.
#' @export
EZVolcanoPlot <- function(obj,
                          parameter = "log_kdeg",
                          design_factor = NULL,
                          reference = NULL,
                          experimental = NULL,
                          param_name = NULL,
                          param_function = NULL,
                          features = NULL,
                          condition = NULL,
                          normalize_by_median = NULL,
                          repeatID = NULL,
                          exactMatch = TRUE,
                          plotlog2 = TRUE,
                          FDR_cutoff = 0.05,
                          difference_cutoff = log(2),
                          size = NULL){

  # Check for backwards compatibility
  if(!is.null(condition) & is.null(design_factor)){
    design_factor <- condition
  }

  ### Hack to deal with devtools::check() NOTEs
  comparison_name <- difference <- padj <- conclusion <- NULL


  ### Find the object you want to get

  # Function in Helpers.R
  comparison_name <- EZget(obj, type = "comparisons",
                           parameter = parameter,
                           features = features,
                          design_factor = design_factor,
                          reference = reference,
                          experimental = experimental,
                          normalize_by_median = normalize_by_median,
                          param_name = param_name,
                          param_function = param_function,
                          repeatID = repeatID,
                          exactMatch = exactMatch,
                          returnNameOnly = TRUE)

  comparison <- obj[['comparisons']][[comparison_name]]
  metadata <- obj[['metadata']][['comparisons']][[comparison_name]]

  # Infer x-axis scale and label
  if(plotlog2){

    scale_factor <- log2(exp(1))

    if(parameter == "log_kdeg"){

      axislabel <- "L2FC(kdeg)"

    }else if(parameter == "log_ksyn"){

      axislabel <- "L2FC(ksyn)"

    }else{

      axislabel <- "Difference"

    }
  }else{

    scale_factor <- 1

    if(parameter == "log_kdeg"){

      axislabel <- "LFC(kdeg)"

    }else if(parameter == "log_ksyn"){

      axislabel <- "LFC(ksyn)"

    }else{

      axislabel <- "Difference"

    }
  }


  # Infer point size if necessary
  if(is.null(size)){

    nr <- nrow(comparison)

    if(nr < 1000){
      size <- 1
    }else if (nr > 10000){
      size <- 0.5
    }else{

      size <- (0.5/(9000))*nr + (4/9)

    }

  }

  # Make Volcano plot
  colors <- c("Decreased" = "deepskyblue4",
              "Increased" = "darkorange",
              "Not Sig." = "darkgray")
  ggv <- comparison %>%
    dplyr::mutate(conclusion = factor(dplyr::case_when(
      difference < -abs(difference_cutoff) & padj < FDR_cutoff ~ "Decreased",
      difference > abs(difference_cutoff) & padj < FDR_cutoff ~ "Increased",
      .default = "Not Sig."
    ), levels = c("Decreased", "Increased", "Not Sig."))) %>%
    ggplot(aes(x = difference*scale_factor, y = -log10(padj), color = conclusion)) +
    geom_point(size = size) +
    scale_color_manual(values = colors) +
    theme_classic() +
    xlab(axislabel) +
    ylab("-log10(padj)") +
    geom_hline(yintercept = -log10(FDR_cutoff),
               color = "darkred",
               linetype = "dotted",
               linewidth = 0.75)

  if(abs(difference_cutoff) > 0){

    intercept <- difference_cutoff * ifelse(plotlog2,
                                            log2(exp(1)),
                                            1)

    ggv <- ggv +
      geom_vline(xintercept = -intercept,
                 color = 'darkred',
                 linetype = "dotted",
                 linewidth = 0.75) +
      geom_vline(xintercept = intercept,
                 color = 'darkred',
                 linetype = "dotted",
                 linewidth = 0.75)

  }


  return(ggv)
}

#' Make an MAPlot from EZbakR comparison
#'
#' Make a plot of effect size (y-axis) vs. log10(read coverage) (x-axis),
#' coloring points by position relative to user-defined decision cutoffs.
#'
#' `EZMAPlot()` accepts as input the output of `CompareParameters()`, i.e.,
#' an `EZbakRData` object with at least one "comparisons" table. It will plot
#' the "avg_coverage" column in this table vs. the "difference" column.
#' In the simplest case, "difference" represents a log-fold change in a kinetic
#' parameter (e.g., kdeg) estimate. More complicated linear model fits and
#' comparisons can yield different parameter estimates.
#'
#' NOTE: some outputs of `CompareParameters()` are not meant for visualization
#' via an MA plot. For example, when fitting certain interaction models,
#' some of the parameter estimates may represent average log(kinetic paramter)
#' in one condition. See discussion of one example of this [here](https://github.com/simonlabcode/bakR/issues/7#issuecomment-2371431127).
#'
#' EZbakR estimates kinetic parameters in `EstimateKinetics()` and `EZDynamics()`
#' on a log-scale. By default, since log2-fold changes are a bit easier to interpret
#' and more common for these kind of visualizations, `EZMAPlot()` multiplies
#' the y-axis value by log2(exp(1)), which is the factor required to convert from
#' a log to a log2 scale. You can turn this off by setting `plotlog2` to `FALSE`.
#'
#' @param obj An object of class `EZbakRCompare`, which is an `EZbakRData` object
#' on which you have run `CompareParameters`
#' @param parameter Name of parameter whose comparison you want to plot.
#' @param design_factor Name of factor from `metadf` whose parameter estimates at
#' different factor values you would like to compare.
#' @param reference Name of reference `condition` factor level value.
#' @param experimental Name of `condition` factor level value to compare to reference.
#' @param param_name If you want to assess the significance of a single parameter,
#' rather than the comparison of two parameters, specify that one parameter's name
#' here.
#' @param param_function NOT YET IMPLEMENTED. Will allow you to specify more complicated
#' functions of parameters when hypotheses you need to test are combinations of parameters
#' rather than individual parameters or simple differences in two parameters.
#' @param features Character vector of feature names for which comparisons were made.
#' @param condition Defunct parameter that has been replaced with `design_factor`. If provided
#' gets passed to `design_factor` if `design_factor` is not already specified.
#' @param repeatID If multiple `kinetics` or `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param plotlog2 If TRUE, assume that log(parameter) difference is passed in and that
#' you want to plot log2(parameter) difference.
#' @param FDR_cutoff False discovery cutoff by which to color points.
#' @param difference_cutoff Minimum absolute difference cutoff by which to color points.
#' @param size Size of points, passed to `geom_point()` size parameter. If not specified,
#' a point size is automatically chosen.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @examples
#'
#' # Simulate data to analyze
#' simdata <- EZSimulate(30)
#'
#' # Create EZbakR input
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#'
#' # Estimate Fractions
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' # Estimate Kinetics
#' ezbdo <- EstimateKinetics(ezbdo)
#'
#' # Average estimates across replicate
#' ezbdo <- AverageAndRegularize(ezbdo)
#'
#' # Compare parameters across conditions
#' ezbdo <- CompareParameters(
#' ezbdo,
#' design_factor = "treatment",
#' reference = "treatment1",
#' experimental = "treatment2"
#' )
#'
#' # Make MA plot (ggplot object that you can save and add/modify layers)
#' EZMAPlot(ezbdo)
#'
#' @return A `ggplot2` object. Y-axis = log2(estimate of interest (e.g., fold-change
#' in degradation rate constant); X-axis = log10(average normalized read coverage);
#' points colored by location relative to FDR and effect size cutoffs.
#' @export
EZMAPlot <- function(obj,
                         parameter = "log_kdeg",
                     design_factor = NULL,
                     reference = NULL,
                     experimental = NULL,
                     param_name = NULL,
                     param_function = NULL,
                         features = NULL,
                     condition = NULL,
                     repeatID = NULL,
                     exactMatch = TRUE,
                         plotlog2 = TRUE,
                         FDR_cutoff = 0.05,
                     difference_cutoff = log(2),
                         size = NULL){


  # Check for backwards compatibility
  if(!is.null(condition) & is.null(design_factor)){
    design_factor <- condition
  }

  ### Hack to deal with devtools::check() NOTEs
  comparison_name <- difference <- avg_coverage <- conclusion <- NULL


  ### Find the object you want to get

  # Function in Helpers.R
  comparison_name <- EZget(obj, type = "comparisons",
                           parameter = parameter,
                           features = features,
                           design_factor = design_factor,
                           reference = reference,
                           experimental = experimental,
                           param_name = param_name,
                           param_function = param_function,
                           repeatID = repeatID,
                           exactMatch = exactMatch,
                           returnNameOnly = TRUE)

  comparison <- obj[['comparisons']][[comparison_name]]
  metadata <- obj[['metadata']][['comparisons']][[comparison_name]]

  if(!("avg_coverage" %in% colnames(comparison))){

    stop("Can't make MA plot if you have no coverage information! This likely means
         that you are trying to compare parameters estimated by EZDynamics, where
         average coverage is not included as output due to it lacking an interpretable
         meaning.")

  }

  # Infer x-axis scale and label
  if(plotlog2){

    scale_factor <- log2(exp(1))

    if(parameter == "log_kdeg"){

      axislabel <- "L2FC(kdeg)"

    }else if(parameter == "log_ksyn"){

      axislabel <- "L2FC(ksyn)"

    }else{

      axislabel <- "Difference"

    }
  }else{

    scale_factor <- 1

    if(parameter == "log_kdeg"){

      axislabel <- "LFC(kdeg)"

    }else if(parameter == "log_ksyn"){

      axislabel <- "LFC(ksyn)"

    }else{

      axislabel <- "Difference"

    }
  }


  # Infer point size if necessary
  if(is.null(size)){

    nr <- nrow(comparison)

    if(nr < 1000){
      size <- 1
    }else if (nr > 10000){
      size <- 0.5
    }else{

      size <- (0.5/(9000))*nr + (4/9)

    }

  }

  # Make Volcano plot
  colors <- c("Decreased" = "deepskyblue4",
              "Increased" = "darkorange",
              "Not Sig." = "darkgray")
  ggma <- comparison %>%
    dplyr::mutate(conclusion = factor(dplyr::case_when(
      difference < -abs(difference_cutoff) & padj < FDR_cutoff ~ "Decreased",
      difference > abs(difference_cutoff) & padj < FDR_cutoff ~ "Increased",
      .default = "Not Sig."
    ), levels = c("Decreased", "Increased", "Not Sig."))) %>%
    ggplot(aes(y = difference*scale_factor, x = avg_coverage, color = conclusion)) +
    geom_point(size = size) +
    scale_color_manual(values = colors) +
    theme_classic() +
    xlab("Avgerage coverage") +
    ylab(axislabel)


  return(ggma)
}


# EZbakRPCAPlot <- function(obj){
#
#   return(ggpca)
# }
