#' Make a VolcanoPlot from EZbakR comparison
#'
#' @param obj An object of class `EZbakRCompare`, which is an `EZbakRData` object
#' on which you have run `CompareParameters`
#' @param parameter Name of parameter whose comparison you want to plot.
#' @param condition Name of condition for which you want to plot the comparison of
#' two of its levels. Condition should be the name of a column that appears in
#' your metadf table.
#' @param reference Name of reference level for comparison.
#' @param experimental Name of experimental level for comparison.
#' @param features Character vector of feature names for which comparisons were made.
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
                              condition = NULL,
                              reference = NULL,
                              experimental = NULL,
                              features = NULL,
                              plotlog2 = TRUE,
                              FDR_cutoff = 0.05,
                              difference_cutoff = 0,
                              size = NULL){


  ### Hack to deal with devtools::check() NOTEs
  comparison_name <- difference <- padj <- conclusion <- NULL


  ### Find the object you want to get

  # Function in Helpers.R
  comparison_name <- EZget(obj, type = "comparisons",
                          condition = condition,
                          reference = reference,
                          experimental = experimental,
                          features = features,
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
  ggv <- comparison %>%
    dplyr::mutate(conclusion = factor(dplyr::case_when(
      difference < -abs(difference_cutoff) & padj < FDR_cutoff ~ "Decreased",
      difference > abs(difference_cutoff) & padj < FDR_cutoff ~ "Increased",
      .default = "Not Sig."
    ), levels = c("Decreased", "Increased", "Not Sig."))) %>%
    ggplot(aes(x = difference*scale_factor, y = -log10(padj), color = conclusion)) +
    geom_point(size = size) +
    scale_color_manual(values = c("deepskyblue4",
                                  "darkorange",
                                  "darkgray")) +
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
#' @param condition Name of condition for which you want to plot the comparison of
#' two of its levels. Condition should be the name of a column that appears in
#' your metadf table.
#' @param reference Name of reference level for comparison.
#' @param experimental Name of experimental level for comparison.
#' @param features Character vector of feature names for which comparisons were made.
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
                         condition = NULL,
                         reference = NULL,
                         experimental = NULL,
                         features = NULL,
                         plotlog2 = TRUE,
                         FDR_cutoff = 0.05,
                         difference_cutoff = 0,
                         size = NULL){

  ### Hack to deal with devtools::check() NOTEs
  comparison_name <- difference <- avg_coverage <- conclusion <- NULL


  ### Find the object you want to get

  # Function in Helpers.R
  comparison_name <- EZget(obj, type = "comparisons",
                          condition = condition,
                          reference = reference,
                          experimental = experimental,
                          features = features,
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
  ggma <- comparison %>%
    dplyr::mutate(conclusion = factor(dplyr::case_when(
      difference < -abs(difference_cutoff) & padj < FDR_cutoff ~ "Decreased",
      difference > abs(difference_cutoff) & padj < FDR_cutoff ~ "Increased",
      .default = "Not Sig."
    ), levels = c("Decreased", "Increased", "Not Sig."))) %>%
    ggplot(aes(y = difference*scale_factor, x = avg_coverage, color = conclusion)) +
    geom_point(size = size) +
    scale_color_manual(values = c("deepskyblue4",
                                  "darkorange",
                                  "darkgray")) +
    theme_classic() +
    xlab("Avgerage coverage") +
    ylab(axislabel)


  return(ggma)
}


# EZbakRPCAPlot <- function(obj){
#
#   return(ggpca)
# }
