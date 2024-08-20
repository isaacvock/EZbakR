#' Make a VolcanoPlot from EZbakR comparison
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
                          repeatID = NULL,
                          exactMatch = TRUE,
                          plotlog2 = TRUE,
                          FDR_cutoff = 0.05,
                          difference_cutoff = 0,
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


  return(ggv)
}

#' Make an MAPlot from EZbakR comparison
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
                         difference_cutoff = 0,
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
