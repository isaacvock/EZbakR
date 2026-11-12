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
#' @param features_to_highlight Features you want to highlight in the plot (black circle
#' will be drawn around them). This can either be a data frame with one column per
#' feature type in the comparison table you are visualizing, or a vector of feature names
#' if the relevant comparison table will only have one feature type noted.
#' @param highlight_shape Shape of points overlayed on highlighted features. Defaults
#' to an open circle
#' @param highlight_size_diff Sets how much larger should the points overlayed on the highlighted
#' features be than the original plot points.
#' @param highlight_stroke Stroke width of the points overlayed on the highlighted
#' features.
#' @param highlight_fill Fill color of the points overlayed on the highlighted
#' features. Default is for them to be fill-less (`highlight_fill == NA`).
#' @param highlight_color Stroke color of the points overlayed on the highlighted points.
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
                          size = NULL,
                          features_to_highlight = NULL,
                          highlight_shape = 21,
                          highlight_size_diff = 1,
                          highlight_stroke = 0.7,
                          highlight_fill = NA,
                          highlight_color = "black"){

  # Check for the sake of backwards compatibility
  if(!is.null(condition) & is.null(design_factor)){
    design_factor <- condition
  }

  ### Hack to deal with devtools::check() NOTEs
  comparison_name <- difference <- padj <- conclusion <- highlight <- NULL


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

  if(is.null(comparison_name)){
    stop("No tables from your EZbakR analysis match your search criteria!")
  }

  comparison <- obj[['comparisons']][[comparison_name]]
  metadata <- obj[['metadata']][['comparisons']][[comparison_name]]


  # Track which features user wants highlighted
  if(!is.null(features_to_highlight)){

    features <- obj[['metadata']][['comparisons']][[comparison_name]][['features']]

    if(is.null(features)){
      # Mainly for backwards compatibility where a bug prevented the features
      # entry of the comparisons metadata from getting updated. Infer feature
      # columns from non-statistic columns of the comparison table.
      features <- setdiff(
        colnames(comparison),
        c("difference", "uncertainty", "stat", "pval", "padj", "avg_coverage")
      )
    }

    if(!is.data.frame(features_to_highlight)){
      # Make into one column data frame if user provided a single vector of
      # feature names

      if(length(features) > 1){

        stop("You provided a vector to features_to_highlight but your comparison table has more than one feature!
            features_to_highlight must be a data frame with columns for each of the feature
            types tracked in the relevant comparison table.")

      }

      features_to_highlight <- dplyr::tibble(
        !!features := as.character(features_to_highlight)
      )

    }

    if(!(all(features %in% colnames(features_to_highlight)))){
      stop("Not all features in comparison table are columns in features_to_highlight!")
    }


    comparison <- comparison %>%
      dplyr::left_join(
        features_to_highlight %>%
          dplyr::mutate(
            highlight = TRUE
          ),
        by = features
      ) %>%
      tidyr::replace_na(
        list(
          highlight = FALSE
        )
      )

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

  if(!is.null(features_to_highlight)){

    ggv <- ggv +
      geom_point(
        data = dplyr::filter(comparison, highlight),
        aes(x = difference*scale_factor, y = -log10(padj)),
        inherit.aes = FALSE,
        shape = highlight_shape,
        size  = size + highlight_size_diff,
        stroke = highlight_stroke,
        fill = highlight_fill,
        color = highlight_color
      )

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
#' @param features_to_highlight Features you want to highlight in the plot (black circle
#' will be drawn around them). This can either be a data frame with one column per
#' feature type in the comparison table you are visualizing, or a vector of feature names
#' if the relevant comparison table will only have one feature type noted.
#' @param highlight_shape Shape of points overlayed on highlighted features. Defaults
#' to an open circle
#' @param highlight_size_diff Sets how much larger should the points overlayed on the highlighted
#' features be than the original plot points.
#' @param highlight_fill Fill color of the points overlayed on the highlighted
#' features. Default is for them to be fill-less (`highlight_fill == NA`).
#' @param highlight_stroke Stroke width of the points overlayed on the highlighted
#' features.
#' @param highlight_color Color of the points overlayed on the highlighted points.
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
                     size = NULL,
                     features_to_highlight = NULL,
                     highlight_shape = 21,
                     highlight_size_diff = 1,
                     highlight_stroke = 0.7,
                     highlight_fill = NA,
                     highlight_color = "black"){


  # Check for backwards compatibility
  if(!is.null(condition) & is.null(design_factor)){
    design_factor <- condition
  }

  ### Hack to deal with devtools::check() NOTEs
  comparison_name <- difference <- avg_coverage <- conclusion <- highlight <- NULL


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

  if(is.null(comparison_name)){
    stop("No tables from your EZbakR analysis match your search criteria!")
  }

  comparison <- obj[['comparisons']][[comparison_name]]
  metadata <- obj[['metadata']][['comparisons']][[comparison_name]]

  if(!("avg_coverage" %in% colnames(comparison))){

    stop("Can't make MA plot if you have no coverage information! This likely means
         that you are trying to compare parameters estimated by EZDynamics, where
         average coverage is not included as output due to it lacking an interpretable
         meaning.")

  }


  # Track which features user wants highlighted
  if(!is.null(features_to_highlight)){

    features <- obj[['metadata']][['comparisons']][[comparison_name]][['features']]

    if(is.null(features)){
      # Mainly for backwards compatibility where a bug prevented the features
      # entry of the comparisons metadata from getting updated. Infer feature
      # columns from non-statistic columns of the comparison table.
      features <- setdiff(
        colnames(comparison),
        c("difference", "uncertainty", "stat", "pval", "padj", "avg_coverage")
      )
    }

    if(!is.data.frame(features_to_highlight)){
      # Make into one column data frame if user provided a single vector of
      # feature names

      if(length(features) > 1){

        stop("You provided a vector to features_to_highlight but your comparison table has more than one feature!
            features_to_highlight must be a data frame with columns for each of the feature
            types tracked in the relevant comparison table.")

      }

      features_to_highlight <- dplyr::tibble(
        !!features := as.character(features_to_highlight)
      )

    }

    if(!(all(features %in% colnames(features_to_highlight)))){
      stop("Not all features in comparison table are columns in features_to_highlight!")
    }


    comparison <- comparison %>%
      dplyr::left_join(
        features_to_highlight %>%
          dplyr::mutate(
            highlight = TRUE
          ),
        by = features
      ) %>%
      tidyr::replace_na(
        list(
          highlight = FALSE
        )
      )

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


  if(abs(difference_cutoff) > 0){

    intercept <- difference_cutoff * ifelse(plotlog2,
                                            log2(exp(1)),
                                            1)

    ggma <- ggma +
      geom_hline(yintercept = -intercept,
                 color = 'darkred',
                 linetype = "dotted",
                 linewidth = 0.75) +
      geom_hline(yintercept = intercept,
                 color = 'darkred',
                 linetype = "dotted",
                 linewidth = 0.75)

  }


  if(!is.null(features_to_highlight)){

    ggma <- ggma +
      geom_point(
        data = dplyr::filter(comparison, highlight),
        aes(y = difference*scale_factor, x = avg_coverage),
        inherit.aes = FALSE,
        shape = highlight_shape,
        size  = size + highlight_size_diff,
        stroke = highlight_stroke,
        fill = highlight_fill,
        color = highlight_color
      )

  }


  return(ggma)
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
#' @param data_type Specifies what data to use for the PCA. Options are "fraction_labeled" (default;
#' means using fraction high T-to-C or other mutation type estimate from EZbakR) or "reads" (means
#' using log10(read counts + 1)).
#' @param features Character vector of feature names for which comparisons were made.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param variance_decile Integer between (inclusive) 1 and 9. Features with sample-to-sample
#' variance greater than the nth decile (n = `variance_decile`) will go into PCA.
#' @param center From `prcomp()`: a logical value indicating whether the variables should be shifted to be zero centered.
#' Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.
#' @param scale From `prcomp()`: a logical value indicating whether the variables should be scaled to have unit variance
#'  before the analysis takes place. Alternatively, a vector of length equal the number of columns of x can be supplied.
#'  The value is passed to scale.
#' @param point_size Size of points in PCA plot
#' @param metadf_cols_to_use Columns in the EZbakR metadf that will be used to
#' color points in the PCA plot. Points will be colored by the interaction between
#' all of these columns (i.e., samples with unique combinations of values of these columns
#' will get unique colors). Default is to use all columns (except "sample"), specified
#' as "all".
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
#'
#' # Make MA plot (ggplot object that you can save and add/modify layers)
#' EZpcaPlot(ezbdo)
#'
#' @return A `ggplot2` object.
#' @export
EZpcaPlot <- function(obj, data_type = c("fraction_labeled", "reads"),  features = NULL, exactMatch = TRUE,
                      variance_decile = 7, center = TRUE, scale = TRUE,
                      point_size = 3, metadf_cols_to_use = "all"){

  repeatID <- tl <- tpulse <- n <- variance <- combo <- NULL

  data_type <- match.arg(data_type)

  # Bind variables locally to resolve devtools::check() Notes
  PC1 <- PC2 <- NULL

  ### Extract logit(fn)

  fractions_name <- EZget(obj,
                          features = features,
                          repeatID = repeatID,
                          exactMatch = exactMatch,
                          returnNameOnly = TRUE)

  if(is.null(fractions_name)){
    stop("No tables from your EZbakR analysis match your search criteria!")
  }

  fraction_table <- obj[["fractions"]][[fractions_name]]


  ### Get useful metadata
  features_to_analyze <- obj[["metadata"]][["fractions"]][[fractions_name]][["features"]]

  fraction_cols <- colnames(fraction_table)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]
  logit_col <- paste0("logit_", fraction_of_interest)


  ### Sample information table used in plot annotation

  sample_lookup <- obj$metadf


  if(length(metadf_cols_to_use) == 1 && metadf_cols_to_use == "all"){
    metadf_cols_to_use <- colnames(obj$metadf)[colnames(obj$metadf) != "sample"]
  }

  if(!(all(metadf_cols_to_use %in% colnames(obj$metadf)))){
    stop("Some of the metadf_cols_to_use specified are not in your metadf!")
  }




  if("tl" %in% colnames(sample_lookup)){

    nlabel_fed <- sample_lookup %>%
      dplyr::filter(
        tl > 0
      ) %>%
      nrow()

  }else if("tpulse" %in% colnames(sample_lookup)){

    nlabel_fed <- sample_lookup %>%
      dplyr::filter(
        tpulse > 0
      ) %>%
      nrow()

  }

  ### Make logit(fn) matrix


  data_col <- ifelse(
    data_type == "fraction_labeled",
    logit_col,
    "log10_n"
  )


  if(data_col == logit_col){

    logit_fn_df <- fraction_table %>%
      dplyr::filter(
        !!dplyr::sym(fraction_of_interest) > 0
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(features_to_analyze)
        )
      ) %>%
      dplyr::filter(
        dplyr::n() == nlabel_fed
      )

  }else{

    logit_fn_df <- fraction_table %>%
      dplyr::ungroup() %>%
      dplyr::select(
        sample, !!features_to_analyze, n
      ) %>%
      tidyr::pivot_wider(
        values_from = "n",
        names_from = "sample",
        values_fill = 0
      ) %>%
      tidyr::pivot_longer(
        cols = -features_to_analyze,
        names_to = "sample",
        values_to = "n"
      ) %>%
      dplyr::mutate(
        log10_n = log10(n + 1)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(features_to_analyze)
        )
      )

  }


  logit_fn_df <- logit_fn_df %>%
    dplyr::mutate(
      variance = stats::var(!!dplyr::sym(data_col))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      variance > stats::quantile(variance, probs = seq(0, 1, 0.1))[variance_decile]
    ) %>%
    dplyr::select(
      sample, !!features_to_analyze, !!data_col
    ) %>%
    tidyr::pivot_wider(
      values_from = !!data_col,
      names_from = "sample"
    )

  logit_fn_mat <- logit_fn_df %>%
    dplyr::select(
      -!!features_to_analyze
    )


  ### Perform PCA

  fn_pcs <- stats::prcomp(t(logit_fn_mat), center = center, scale. = scale)


  ### Extract loadings and PCs

  fn_eigenvect <- fn_pcs$x

  fn_PC1 <- fn_eigenvect[,c("PC1")]
  fn_PC2 <- fn_eigenvect[,c("PC2")]


  meta_cols <- colnames(obj$metadf)[colnames(obj$metadf) != "sample"]


  samps_to_keep <- colnames(logit_fn_mat)
  fn_pca_df <- dplyr::tibble(sample = rownames(fn_eigenvect),
                             PC1 = fn_PC1,
                             PC2 = fn_PC2) %>%
    dplyr::left_join(obj$metadf %>%
                       dplyr::filter(
                         sample %in% samps_to_keep
                       ) %>%
                       dplyr::select(
                         sample, !!metadf_cols_to_use
                       ),
                     by = "sample") %>%
    dplyr::mutate(
      combo = interaction(
        dplyr::across(
          dplyr::all_of(metadf_cols_to_use)
        ),
        drop = TRUE, lex.order = TRUE, sep = ":"
      )
    )

  fn_prop_var <- unclass(summary(fn_pcs))$importance[c("Proportion of Variance"),]

  ### Plot

  legend_title <- paste0(meta_cols[meta_cols %in% metadf_cols_to_use], collapse=":")

  g_pca <- ggplot(fn_pca_df, aes(x = PC1, y = PC2, color = combo)) +
    geom_point(size = point_size) +
    xlab(paste0("PC1 (", fn_prop_var[1]*100, "% of var.)")) +
    ylab(paste0("PC2 (", fn_prop_var[2]*100, "% of var.)")) +
    theme_classic() +
    scale_color_viridis_d() +
    guides(color = guide_legend(title=legend_title))




  return(g_pca)

}


#' Visualize dataset-wide and sample-level goodness of fit
#'
#' `EZGoFPlot()` provides four complementary diagnostic views of model goodness of fit
#' computed by `AssessGoF()` (and will automatically run `AssessGoF()` on the fly if
#' no GoF table is present yet):
#' \itemize{
#'   \item `"volcano"` (default): Feature-level triage plotting probability-scale error
#'   (`wrmse_fn`, or `mean_signed_resid_fn` if `signed = TRUE`) against `-log10(padj)`
#'   for structural lack of fit (or replicate overdispersion when single-timepoint
#'   designs have zero structural degrees of freedom).
#'   \item `"sample_qc"`: Per-sample boxplots of residuals (`obs_fn - pred_fn`,
#'   standardized residuals `std_resid`, or logit residuals), immediately flagging
#'   outlier samples, mis-specified label times, s4U dropout, or compartment contamination.
#'   \item `"obs_vs_pred"`: Global parity scatter/density plot of observed vs. predicted
#'   fractions (or logit fractions) against the `y = x` identity line.
#'   \item `"residual_trend"`: Median and interquartile ribbon of residuals across
#'   label time `tl` (or predicted fraction), revealing transcriptome-wide splicing/export
#'   lags or long-label toxicity.
#' }
#'
#' @param obj An `EZbakRData` object on which `AssessGoF()` has been run (or containing
#' `EstimateFractions()` + `AverageAndRegularize()` / `EZDynamics()` outputs).
#' @param plot_type Type of diagnostic plot: `"volcano"`, `"sample_qc"`, `"obs_vs_pred"`,
#' or `"residual_trend"`.
#' @param source_type Optional filter (`"averages"` or `"dynamics"`) selecting which GoF
#' table to plot when multiple exist.
#' @param source_table Optional exact upstream table name passed to `EZget()`.
#' @param features Optional feature vector passed to `EZget()`.
#' @param parameter Optional parameter name passed to `EZget()`.
#' @param kstrat Optional kinetic strategy passed to `EZget()`.
#' @param graph Optional `EZDynamics` graph matrix passed to `EZget()`.
#' @param grouping_features Optional `grouping_features` vector passed to `EZget()`.
#' @param repeatID Optional repeat ID passed to `EZget()`.
#' @param exactMatch Logical; passed to `EZget()` (default `FALSE` so partial queries match cleanly).
#' @param metric Which adjusted p-value to display on `"volcano"`: `"auto"` (default;
#' prioritizes `padj_lof`, then `padj_lof_joint`, then `padj_rep`), `"lof"`, `"joint"`, or `"rep"`.
#' @param FDR_cutoff Adjusted p-value threshold for highlighting poor fits on `"volcano"` (default `0.05`).
#' @param rmse_cutoff Fraction-scale weighted RMSE (or `|mean_signed_resid_fn|`) threshold
#' for highlighting poor fits on `"volcano"` (default `0.05`, i.e., 5 percentage points).
#' @param signed Logical (default `FALSE`); if `TRUE` on `"volcano"`, plots signed mean
#' residual (`obs - pred`) on X instead of `wrmse_fn`.
#' @param size Optional point size for scatter plots.
#' @param features_to_highlight Optional vector of feature IDs (or data frame of feature
#' key columns) to circle on `"volcano"`.
#' @param highlight_shape Shape of highlight overlay points (default `21`).
#' @param highlight_size_diff Size offset added to highlighted points (default `1`).
#' @param highlight_stroke Stroke width for highlighted points (default `0.7`).
#' @param highlight_fill Fill color for highlighted points (default `NA`).
#' @param highlight_color Border color for highlighted points (default `"black"`).
#' @param residual_scale Scale for residuals in `"sample_qc"`, `"obs_vs_pred"`, and
#' `"residual_trend"`: `"fraction"` (default), `"standardized"`, or `"logit"`.
#' @param facet_by Optional column name(s) in `gof_residuals` to facet `"obs_vs_pred"`
#' or `"residual_trend"` by.
#' @return A `ggplot2` object.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @examples
#' simdata <- EZSimulate(30, nreps = 2)
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#' ezbdo <- EstimateFractions(ezbdo)
#' ezbdo <- EstimateKinetics(ezbdo)
#' ezbdo <- AverageAndRegularize(ezbdo)
#' ezbdo <- AssessGoF(ezbdo)
#' EZGoFPlot(ezbdo, plot_type = "volcano")
#' EZGoFPlot(ezbdo, plot_type = "sample_qc")
#' @export
EZGoFPlot <- function(obj,
                      plot_type = c("volcano", "sample_qc", "obs_vs_pred", "residual_trend"),
                      source_type = NULL,
                      source_table = NULL,
                      features = NULL,
                      parameter = NULL,
                      kstrat = NULL,
                      graph = NULL,
                      grouping_features = NULL,
                      repeatID = NULL,
                      exactMatch = FALSE,
                      metric = c("auto", "lof", "joint", "rep"),
                      FDR_cutoff = 0.05,
                      rmse_cutoff = 0.05,
                      signed = FALSE,
                      size = NULL,
                      features_to_highlight = NULL,
                      highlight_shape = 21,
                      highlight_size_diff = 1,
                      highlight_stroke = 0.7,
                      highlight_fill = NA,
                      highlight_color = "black",
                      residual_scale = c("fraction", "standardized", "logit"),
                      facet_by = NULL){

  # Silence R CMD check NSE notes
  x_val <- y_val <- conclusion <- highlight <- sample_ord <- resid_val <- NULL
  group_fill <- pred_x <- obs_y <- x_bin <- med_r <- q25_r <- q75_r <- trend_group <- NULL

  plot_type <- match.arg(plot_type)
  metric <- match.arg(metric)
  residual_scale <- match.arg(residual_scale)

  if(length(obj[["metadata"]][["gof"]]) == 0L){
    message("No GoF tables found in `obj`; running `AssessGoF()` on the fly...")
    obj <- AssessGoF(
      obj,
      type = source_type %||% "auto",
      features = features,
      parameter = parameter,
      kstrat = kstrat,
      graph = graph,
      grouping_features = grouping_features,
      repeatID = repeatID
    )
  }

  always_chk <- !is.null(source_type) || !is.null(source_table) || !is.null(features) ||
    !is.null(parameter) || !is.null(kstrat) || !is.null(graph) ||
    !is.null(grouping_features) || !is.null(repeatID)

  gof_name <- EZget(
    obj,
    type = "gof",
    source_type = source_type,
    source_table = source_table,
    features = features,
    parameter = parameter,
    kstrat = kstrat,
    graph = graph,
    grouping_features = grouping_features,
    repeatID = repeatID,
    exactMatch = exactMatch,
    alwaysCheck = always_chk,
    returnNameOnly = TRUE
  )

  if(is.null(gof_name) || length(gof_name) == 0L){
    stop("No `gof` tables in `obj` matched your search criteria!")
  }
  if(length(gof_name) > 1L){
    stop("Multiple `gof` tables matched (", paste(gof_name, collapse = ", "),
         ")! Specify `source_type`, `source_table`, `features`, or `parameter`.")
  }

  gof_tbl <- obj[["gof"]][[gof_name]]
  gof_meta <- obj[["metadata"]][["gof"]][[gof_name]]

  ### 1. VOLCANO / FEATURE TRIAGE PLOT
  if(plot_type == "volcano"){
    if(metric == "auto"){
      if(any(gof_tbl$df_lof > 0L, na.rm = TRUE) && any(!is.na(gof_tbl$padj_lof))){
        eff_metric <- "lof"
      }else if("padj_lof_joint" %in% colnames(gof_tbl) &&
               any(gof_tbl$df_lof_joint > 0L, na.rm = TRUE) &&
               any(!is.na(gof_tbl$padj_lof_joint))){
        eff_metric <- "joint"
      }else if(any(gof_tbl$df_rep > 0L, na.rm = TRUE) && any(!is.na(gof_tbl$padj_rep))){
        eff_metric <- "rep"
      }else{
        stop("Selected GoF table has 0 degrees of freedom for both model lack-of-fit and replicate dispersion!")
      }
    }else{
      eff_metric <- metric
    }

    if(eff_metric == "lof"){
      padj_vec <- gof_tbl$padj_lof
      y_lab <- "-log10(padj) [Model lack of fit]"
    }else if(eff_metric == "joint"){
      if(!("padj_lof_joint" %in% colnames(gof_tbl))){
        stop("`metric = 'joint'` requires an `EZDynamics` GoF table with `include_coverage_gof = TRUE`!")
      }
      padj_vec <- gof_tbl$padj_lof_joint
      y_lab <- "-log10(padj) [Joint fraction + coverage GoF]"
    }else{
      padj_vec <- gof_tbl$padj_rep
      y_lab <- "-log10(padj) [Replicate dispersion]"
    }

    if(all(is.na(padj_vec))){
      stop("All adjusted p-values for `metric = '", eff_metric, "'` are NA (0 degrees of freedom)!")
    }

    df_plot <- gof_tbl
    df_plot$x_val <- if(signed) df_plot$mean_signed_resid_fn else df_plot$wrmse_fn
    df_plot$y_val <- -log10(pmax(padj_vec, 1e-300))
    df_plot <- df_plot[is.finite(df_plot$x_val) & is.finite(df_plot$y_val), , drop = FALSE]

    # Optional feature highlighting (supports vector on 1-feature or grouping_features)
    if(!is.null(features_to_highlight)){
      f_cols <- gof_meta[["features"]]
      if(!is.data.frame(features_to_highlight)){
        if(length(f_cols) == 1L){
          hl_df <- dplyr::tibble(!!f_cols := as.character(features_to_highlight))
          join_hl <- f_cols
        }else if(!is.null(gof_meta[["grouping_features"]]) &&
                 length(gof_meta[["grouping_features"]]) == 1L){
          gf1 <- gof_meta[["grouping_features"]]
          hl_df <- dplyr::tibble(!!gf1 := as.character(features_to_highlight))
          join_hl <- gf1
        }else{
          stop("When `gof` table has multiple feature columns (", paste(f_cols, collapse = ", "),
               "), `features_to_highlight` must be a data frame with those columns!")
        }
      }else{
        join_hl <- intersect(colnames(features_to_highlight), colnames(df_plot))
        if(length(join_hl) == 0L){
          stop("`features_to_highlight` shares no feature columns with the `gof` table!")
        }
        hl_df <- dplyr::distinct(features_to_highlight[, join_hl, drop = FALSE])
      }
      hl_df$highlight <- TRUE
      df_plot <- dplyr::left_join(df_plot, hl_df, by = join_hl)
      df_plot$highlight[is.na(df_plot$highlight)] <- FALSE
    }

    if(is.null(size)){
      nr <- nrow(df_plot)
      size <- if(nr < 1000L) 1.2 else if(nr > 10000L) 0.5 else (0.5 / 9000) * nr + (4 / 9)
    }

    if(signed){
      colors <- c("Overpredicted" = "deepskyblue4",
                  "Underpredicted" = "darkorange",
                  "Good Fit" = "darkgray")
      df_plot$conclusion <- factor(dplyr::case_when(
        df_plot$x_val < -abs(rmse_cutoff) & (10^(-df_plot$y_val)) < FDR_cutoff ~ "Overpredicted",
        df_plot$x_val > abs(rmse_cutoff) & (10^(-df_plot$y_val)) < FDR_cutoff ~ "Underpredicted",
        .default = "Good Fit"
      ), levels = c("Overpredicted", "Underpredicted", "Good Fit"))
      x_lab <- "Mean residual (observed - predicted fraction)"
    }else{
      colors <- c("Poor Fit" = "darkorange",
                  "Good Fit" = "darkgray")
      df_plot$conclusion <- factor(dplyr::case_when(
        df_plot$x_val >= abs(rmse_cutoff) & (10^(-df_plot$y_val)) < FDR_cutoff ~ "Poor Fit",
        .default = "Good Fit"
      ), levels = c("Poor Fit", "Good Fit"))
      x_lab <- "Weighted RMSE (fraction high)"
    }

    gg <- ggplot(df_plot, aes(x = x_val, y = y_val, color = conclusion)) +
      geom_point(size = size) +
      scale_color_manual(values = colors) +
      theme_classic() +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_hline(yintercept = -log10(FDR_cutoff),
                 color = "darkred", linetype = "dotted", linewidth = 0.75)

    if(abs(rmse_cutoff) > 0){
      gg <- gg +
        geom_vline(xintercept = abs(rmse_cutoff),
                   color = "darkred", linetype = "dotted", linewidth = 0.75)
      if(signed){
        gg <- gg +
          geom_vline(xintercept = -abs(rmse_cutoff),
                     color = "darkred", linetype = "dotted", linewidth = 0.75)
      }
    }

    if(!is.null(features_to_highlight)){
      gg <- gg +
        geom_point(
          data = df_plot[df_plot$highlight, , drop = FALSE],
          aes(x = x_val, y = y_val),
          inherit.aes = FALSE,
          shape = highlight_shape,
          size = size + highlight_size_diff,
          stroke = highlight_stroke,
          fill = highlight_fill,
          color = highlight_color
        )
    }

    return(gg)
  }

  ### Residual-table-backed plot modes ("sample_qc", "obs_vs_pred", "residual_trend")
  res_tbl <- obj[["gof_residuals"]][[gof_name]]
  if(is.null(res_tbl)){
    stop("`plot_type = '", plot_type, "'` requires sample-level residuals! Re-run `AssessGoF(..., store_residuals = TRUE)`.")
  }

  # Helper to choose residual column and y-axis label
  if(residual_scale == "fraction"){
    r_col <- "resid_fn"
    r_lab <- "Residual (observed - predicted fraction)"
  }else if(residual_scale == "standardized"){
    r_col <- "std_resid"
    r_lab <- "Standardized residual (profile z-score)"
  }else{
    r_col <- "resid_logit_fn"
    r_lab <- "Logit residual (obs - pred logit fraction)"
  }

  # Identify sample metadata grouping columns for coloring/ordering
  reserved_res_cols <- c("sample", gof_meta[["features"]], gof_meta[["raw_features"]],
                         "measured_specie", "tl", "obs_fn", "pred_fn", "resid_fn",
                         "obs_logit_fn", "pred_logit_fn", "resid_logit_fn", "std_resid",
                         "ci_lower_fn", "ci_upper_fn", "n", "obs_log10_cov",
                         "pred_log10_cov", "resid_log10_cov", "std_resid_cov")
  sf_real <- if(!identical(gof_meta[["sample_feature"]], "imputed_sample_feature")) gof_meta[["sample_feature"]] else NULL
  cond_cols <- unique(c(sf_real, setdiff(colnames(res_tbl), reserved_res_cols)))
  cond_cols <- cond_cols[cond_cols %in% colnames(res_tbl)]

  ### 2. SAMPLE-BY-SAMPLE QC PLOT
  if(plot_type == "sample_qc"){
    df_sq <- res_tbl
    df_sq$resid_val <- as.numeric(df_sq[[r_col]])
    df_sq <- df_sq[is.finite(df_sq$resid_val), , drop = FALSE]

    if(length(cond_cols) > 0L){
      df_sq$group_fill <- apply(df_sq[, cond_cols, drop = FALSE], 1L, paste, collapse = ":")
      fill_title <- paste(cond_cols, collapse = ":")
    }else{
      df_sq$group_fill <- paste0("tl=", df_sq$tl)
      fill_title <- "Label time"
    }

    samp_meta <- dplyr::distinct(df_sq[, c("sample", "group_fill", "tl"), drop = FALSE])
    samp_meta <- samp_meta[order(samp_meta$group_fill, samp_meta$tl, samp_meta$sample), , drop = FALSE]
    u_samps <- unique(samp_meta$sample)
    df_sq$sample_ord <- factor(df_sq$sample, levels = u_samps)

    gg <- ggplot(df_sq, aes(x = sample_ord, y = resid_val, fill = group_fill)) +
      geom_hline(yintercept = 0, color = "darkred", linetype = "dotted", linewidth = 0.8) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.35, alpha = 0.8) +
      scale_fill_viridis_d() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Sample") +
      ylab(r_lab) +
      guides(fill = guide_legend(title = fill_title))

    if(residual_scale == "standardized"){
      gg <- gg +
        geom_hline(yintercept = c(-1.96, 1.96), color = "gray40", linetype = "dashed", linewidth = 0.5)
    }

    return(gg)
  }

  ### 3. OBSERVED VS PREDICTED PARITY PLOT
  if(plot_type == "obs_vs_pred"){
    df_op <- res_tbl
    if(residual_scale == "logit"){
      df_op$pred_x <- as.numeric(df_op$pred_logit_fn)
      df_op$obs_y <- as.numeric(df_op$obs_logit_fn)
      x_lab <- "Model-predicted logit(fraction high)"
      y_lab <- "Observed logit(fraction high)"
    }else{
      df_op$pred_x <- as.numeric(df_op$pred_fn)
      df_op$obs_y <- as.numeric(df_op$obs_fn)
      x_lab <- "Model-predicted fraction high"
      y_lab <- "Observed fraction high"
    }
    df_op <- df_op[is.finite(df_op$pred_x) & is.finite(df_op$obs_y), , drop = FALSE]

    if(is.null(facet_by) && identical(gof_meta[["source_type"]], "dynamics")){
      auto_f <- c(sf_real, if(length(unique(df_op$measured_specie)) > 1L) "measured_specie" else NULL)
      if(length(auto_f) > 0L) facet_by <- auto_f
    }

    pt_sz <- size %||% (if(nrow(df_op) > 20000L) 0.4 else 0.8)
    pt_al <- if(nrow(df_op) > 20000L) 0.15 else 0.3

    gg <- ggplot(df_op, aes(x = pred_x, y = obs_y)) +
      geom_point(size = pt_sz, alpha = pt_al, color = "deepskyblue4") +
      geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed", linewidth = 0.8) +
      theme_classic() +
      xlab(x_lab) +
      ylab(y_lab)

    if(residual_scale != "logit"){
      gg <- gg + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    }

    if(!is.null(facet_by)){
      miss_f <- setdiff(facet_by, colnames(df_op))
      if(length(miss_f) > 0L){
        stop("`facet_by` column(s) not found in `gof_residuals`: ", paste(miss_f, collapse = ", "))
      }
      gg <- gg + facet_wrap(stats::as.formula(paste("~", paste(facet_by, collapse = " + "))),
                            labeller = label_both)
    }

    return(gg)
  }

  ### 4. RESIDUAL TREND PLOT (across tl or pred_fn)
  df_rt <- res_tbl
  df_rt$resid_val <- as.numeric(df_rt[[r_col]])
  df_rt <- df_rt[is.finite(df_rt$resid_val), , drop = FALSE]

  grp_cols_rt <- character(0)
  if(identical(gof_meta[["source_type"]], "dynamics")){
    if(!is.null(sf_real) && length(unique(df_rt[[sf_real]])) > 1L){
      grp_cols_rt <- c(grp_cols_rt, sf_real)
    }
    if("measured_specie" %in% colnames(df_rt) && length(unique(df_rt$measured_specie)) > 1L){
      grp_cols_rt <- c(grp_cols_rt, "measured_specie")
    }
  }else if(length(cond_cols) > 0L){
    grp_cols_rt <- cond_cols
  }

  if(length(grp_cols_rt) > 0L){
    df_rt$trend_group <- apply(df_rt[, grp_cols_rt, drop = FALSE], 1L, paste, collapse = ":")
    grp_title <- paste(grp_cols_rt, collapse = ":")
  }else{
    df_rt$trend_group <- "All observations"
    grp_title <- "Group"
  }

  multi_tl <- length(unique(df_rt$tl)) > 1L
  if(multi_tl){
    df_rt$x_bin <- as.numeric(df_rt$tl)
    x_lab <- "Label time (tl)"
  }else{
    # Bin pred_fn into deciles/bins when single label time
    breaks <- seq(0, 1, length.out = 11L)
    bin_id <- findInterval(pmin(pmax(df_rt$pred_fn, 0), 1), breaks, rightmost.closed = TRUE)
    mids <- 0.5 * (breaks[-length(breaks)] + breaks[-1L])
    df_rt$x_bin <- mids[bin_id]
    x_lab <- "Model-predicted fraction high (binned)"
  }

  group_vars_sum <- c("trend_group", "x_bin", facet_by)
  group_vars_sum <- unique(group_vars_sum[group_vars_sum %in% colnames(df_rt)])

  sum_rt <- df_rt %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_sum))) %>%
    dplyr::summarise(
      med_r = stats::median(resid_val, na.rm = TRUE),
      q25_r = stats::quantile(resid_val, probs = 0.25, na.rm = TRUE, names = FALSE),
      q75_r = stats::quantile(resid_val, probs = 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  gg <- ggplot(sum_rt, aes(x = x_bin, y = med_r, color = trend_group, fill = trend_group)) +
    geom_hline(yintercept = 0, color = "darkred", linetype = "dotted", linewidth = 0.8) +
    geom_ribbon(aes(ymin = q25_r, ymax = q75_r), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_classic() +
    xlab(x_lab) +
    ylab(paste0(r_lab, " (median & IQR)")) +
    guides(color = guide_legend(title = grp_title),
           fill = guide_legend(title = grp_title))

  if(!is.null(facet_by)){
    gg <- gg + facet_wrap(stats::as.formula(paste("~", paste(facet_by, collapse = " + "))),
                          labeller = label_both)
  }

  return(gg)
}


#' Plot model-predicted trajectories and sample-level fraction fits for selected features
#'
#' `EZFitPlot()` overlays smooth model-predicted fraction high `<muttype>` curves
#' (from either `AverageAndRegularize()` kinetic models or `EZDynamics()` ODE reaction
#' graphs) onto the individual sample estimates and 95\% profile-likelihood confidence
#' intervals from `EstimateFractions()`.
#'
#' Features can either be specified explicitly via `features_to_plot` or automatically
#' ranked via `rank_by` (`"worst_lof"`, `"worst_wrmse"`, `"worst_rep"`, or `"best_fit"`).
#' When `x_axis = "auto"`, time-course and `EZDynamics()` models plot continuous
#' trajectories over label time `tl`, whereas single-timepoint multi-condition
#' `AverageAndRegularize()` fits display replicate points alongside condition-level
#' model predictions on a categorical condition axis.
#'
#' @param obj An `EZbakRData` object containing `AssessGoF()` output (or upstream
#' `EstimateFractions()` + `AverageAndRegularize()` / `EZDynamics()` tables).
#' @param features_to_plot Optional character vector of feature names (when a single
#' feature or grouping feature identifies genes) or data frame of feature columns to plot.
#' If `NULL` (default), selects the top `top_n` features according to `rank_by`.
#' @param rank_by Automatic feature ranking criterion when `features_to_plot = NULL`:
#' `"worst_lof"` (default; smallest structural lack-of-fit p-value), `"worst_wrmse"`
#' (largest weighted fraction RMSE), `"worst_rep"` (smallest replicate dispersion p-value),
#' or `"best_fit"` (smallest weighted fraction RMSE among well-covered features).
#' @param top_n Number of features to display when `features_to_plot = NULL` (default `6`).
#' @param source_type Optional filter (`"averages"` or `"dynamics"`) passed to `EZget()`.
#' @param source_table Optional upstream table name passed to `EZget()`.
#' @param features Optional feature vector passed to `EZget()`.
#' @param parameter Optional parameter name passed to `EZget()`.
#' @param kstrat Optional kinetic strategy passed to `EZget()`.
#' @param graph Optional `EZDynamics` graph matrix passed to `EZget()`.
#' @param grouping_features Optional `grouping_features` vector passed to `EZget()`.
#' @param repeatID Optional repeat ID passed to `EZget()`.
#' @param exactMatch Logical; passed to `EZget()` (default `FALSE`).
#' @param x_axis Layout mode: `"auto"` (default), `"time"` (continuous label-time
#' trajectory curves starting at `t = 0`), or `"condition"` (categorical condition comparison).
#' @param n_grid Number of time-grid evaluation points for smooth trajectory curves (default `150`).
#' @param show_ci Logical (default `TRUE`); if `TRUE`, draws 95\% profile-likelihood
#' confidence intervals on observed sample points.
#' @param annotate_stats Logical (default `TRUE`); if `TRUE`, includes `wRMSE` in each
#' facet strip header.
#' @param point_size Point size for observed sample fractions (default `2.2`).
#' @param line_width Line width for model trajectory curves (default `0.9`).
#' @param ncol Optional number of facet columns passed to `facet_wrap()`.
#' @return A `ggplot2` object.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @examples
#' simdata <- EZSimulate(30, nreps = 2)
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#' ezbdo <- EstimateFractions(ezbdo)
#' ezbdo <- EstimateKinetics(ezbdo)
#' ezbdo <- AverageAndRegularize(ezbdo)
#' ezbdo <- AssessGoF(ezbdo)
#' EZFitPlot(ezbdo, top_n = 4)
#' @export
EZFitPlot <- function(obj,
                      features_to_plot = NULL,
                      rank_by = c("worst_lof", "worst_wrmse", "worst_rep", "best_fit"),
                      top_n = 6,
                      source_type = NULL,
                      source_table = NULL,
                      features = NULL,
                      parameter = NULL,
                      kstrat = NULL,
                      graph = NULL,
                      grouping_features = NULL,
                      repeatID = NULL,
                      exactMatch = FALSE,
                      x_axis = c("auto", "time", "condition"),
                      n_grid = 150,
                      show_ci = TRUE,
                      annotate_stats = TRUE,
                      point_size = 2.2,
                      line_width = 0.9,
                      ncol = NULL){

  # Silence R CMD check NSE notes
  tl <- pred_fn <- obs_fn <- ci_lower_fn <- ci_upper_fn <- curve_group <- NULL
  panel_label <- cond_x <- NULL

  rank_by <- match.arg(rank_by)
  x_axis <- match.arg(x_axis)

  if(length(obj[["metadata"]][["gof"]]) == 0L){
    message("No GoF tables found in `obj`; running `AssessGoF()` on the fly...")
    obj <- AssessGoF(
      obj,
      type = source_type %||% "auto",
      features = features,
      parameter = parameter,
      kstrat = kstrat,
      graph = graph,
      grouping_features = grouping_features,
      repeatID = repeatID,
      store_residuals = TRUE
    )
  }

  always_chk <- !is.null(source_type) || !is.null(source_table) || !is.null(features) ||
    !is.null(parameter) || !is.null(kstrat) || !is.null(graph) ||
    !is.null(grouping_features) || !is.null(repeatID)

  gof_name <- EZget(
    obj,
    type = "gof",
    source_type = source_type,
    source_table = source_table,
    features = features,
    parameter = parameter,
    kstrat = kstrat,
    graph = graph,
    grouping_features = grouping_features,
    repeatID = repeatID,
    exactMatch = exactMatch,
    alwaysCheck = always_chk,
    returnNameOnly = TRUE
  )

  if(is.null(gof_name) || length(gof_name) == 0L){
    stop("No `gof` tables in `obj` matched your search criteria!")
  }
  if(length(gof_name) > 1L){
    stop("Multiple `gof` tables matched (", paste(gof_name, collapse = ", "),
         ")! Specify `source_type`, `source_table`, `features`, or `parameter`.")
  }

  gof_tbl <- obj[["gof"]][[gof_name]]
  gof_meta <- obj[["metadata"]][["gof"]][[gof_name]]
  res_tbl <- obj[["gof_residuals"]][[gof_name]]
  if(is.null(res_tbl)){
    stop("`EZFitPlot()` requires sample-level residuals! Re-run `AssessGoF(..., store_residuals = TRUE)`.")
  }

  f_cols <- gof_meta[["features"]]

  ### Step 1: Select features to plot (preserving rank / user-specified order)
  if(is.null(features_to_plot)){
    n_take <- min(as.integer(top_n), nrow(gof_tbl))
    if(rank_by == "worst_lof"){
      if(any(!is.na(gof_tbl$pval_lof))){
        ord <- order(is.na(gof_tbl$pval_lof), gof_tbl$pval_lof, -gof_tbl$wrmse_fn)
      }else if("pval_lof_joint" %in% colnames(gof_tbl) && any(!is.na(gof_tbl$pval_lof_joint))){
        ord <- order(is.na(gof_tbl$pval_lof_joint), gof_tbl$pval_lof_joint, -gof_tbl$wrmse_fn)
      }else if(any(!is.na(gof_tbl$pval_rep))){
        ord <- order(is.na(gof_tbl$pval_rep), gof_tbl$pval_rep, -gof_tbl$wrmse_fn)
      }else{
        ord <- order(-gof_tbl$wrmse_fn)
      }
    }else if(rank_by == "worst_wrmse"){
      ord <- order(is.na(gof_tbl$wrmse_fn), -gof_tbl$wrmse_fn)
    }else if(rank_by == "worst_rep"){
      if(any(!is.na(gof_tbl$pval_rep))){
        ord <- order(is.na(gof_tbl$pval_rep), gof_tbl$pval_rep, -gof_tbl$wrmse_fn)
      }else{
        ord <- order(-gof_tbl$wrmse_fn)
      }
    }else{
      ord <- order(is.na(gof_tbl$wrmse_fn), gof_tbl$wrmse_fn)
    }
    sel_gof <- gof_tbl[ord[seq_len(n_take)], , drop = FALSE]
  }else{
    if(!is.data.frame(features_to_plot)){
      if(length(f_cols) == 1L){
        q_df <- dplyr::tibble(!!f_cols := as.character(features_to_plot))
        join_cols <- f_cols
      }else if(!is.null(gof_meta[["grouping_features"]]) &&
               length(gof_meta[["grouping_features"]]) == 1L){
        gf1 <- gof_meta[["grouping_features"]]
        q_df <- dplyr::tibble(!!gf1 := as.character(features_to_plot))
        join_cols <- gf1
      }else{
        stop("When `gof` table has multiple feature columns (", paste(f_cols, collapse = ", "),
             "), `features_to_plot` must be a data frame with those columns!")
      }
    }else{
      join_cols <- intersect(colnames(features_to_plot), f_cols)
      if(length(join_cols) == 0L){
        stop("`features_to_plot` does not contain matching feature columns (", paste(f_cols, collapse = ", "), ")!")
      }
      q_df <- dplyr::distinct(features_to_plot[, join_cols, drop = FALSE])
    }
    gof_chr <- gof_tbl
    for(jc in join_cols) gof_chr[[jc]] <- as.character(gof_chr[[jc]])
    for(jc in join_cols) q_df[[jc]] <- as.character(q_df[[jc]])
    sel_gof <- dplyr::inner_join(q_df, gof_chr, by = join_cols)
    if(nrow(sel_gof) == 0L){
      stop("None of the requested `features_to_plot` were found in the `gof` table!")
    }
  }

  sel_keys <- sel_gof[, f_cols, drop = FALSE]
  key_label_raw <- apply(sel_keys, 1L, paste, collapse = " | ")
  if(annotate_stats){
    panel_lvls <- sprintf("%s (wRMSE=%.3f)", key_label_raw, sel_gof$wrmse_fn)
  }else{
    panel_lvls <- key_label_raw
  }
  panel_map <- sel_keys
  panel_map$panel_label <- factor(panel_lvls, levels = unique(panel_lvls))

  res_chr <- res_tbl
  for(fc in f_cols){
    res_chr[[fc]] <- as.character(res_chr[[fc]])
    panel_map[[fc]] <- as.character(panel_map[[fc]])
  }
  obs_sub <- dplyr::inner_join(res_chr, panel_map, by = f_cols)

  # Identify sample metadata condition columns
  reserved_cols <- c("sample", f_cols, gof_meta[["raw_features"]], "measured_specie",
                     "tl", "obs_fn", "pred_fn", "resid_fn", "obs_logit_fn",
                     "pred_logit_fn", "resid_logit_fn", "std_resid",
                     "ci_lower_fn", "ci_upper_fn", "n", "obs_log10_cov",
                     "pred_log10_cov", "resid_log10_cov", "std_resid_cov", "panel_label")
  sf_real <- if(!identical(gof_meta[["sample_feature"]], "imputed_sample_feature")) gof_meta[["sample_feature"]] else NULL
  cond_cols <- setdiff(colnames(obs_sub), c(reserved_cols, sf_real))

  if(x_axis == "auto"){
    if(identical(gof_meta[["source_type"]], "dynamics") ||
       length(unique(obs_sub$tl)) > 1L ||
       identical(gof_meta[["kstrat"]], "pulse-chase") ||
       length(cond_cols) == 0L){
      eff_xaxis <- "time"
    }else{
      eff_xaxis <- "condition"
    }
  }else{
    eff_xaxis <- x_axis
  }

  ### Layout A: Categorical condition X-axis (single-tl multi-condition averages)
  if(eff_xaxis == "condition"){
    if(length(cond_cols) > 0L){
      obs_sub$cond_x <- apply(obs_sub[, cond_cols, drop = FALSE], 1L, paste, collapse = ":")
      x_title <- paste(cond_cols, collapse = ":")
    }else{
      obs_sub$cond_x <- "All samples"
      x_title <- "Condition"
    }

    pred_cond <- dplyr::distinct(obs_sub[, c("panel_label", "cond_x", "pred_fn"), drop = FALSE])
    pos_jd <- position_jitter(width = 0.12, height = 0, seed = 1L)

    gg <- ggplot(obs_sub, aes(x = cond_x, y = obs_fn, color = cond_x)) +
      geom_errorbar(
        data = if(show_ci) obs_sub else obs_sub[0, , drop = FALSE],
        aes(ymin = ci_lower_fn, ymax = ci_upper_fn),
        width = 0.15, alpha = 0.55, position = pos_jd
      ) +
      geom_point(size = point_size, alpha = 0.85, position = pos_jd) +
      geom_crossbar(
        data = pred_cond,
        aes(x = cond_x, y = pred_fn, ymin = pred_fn, ymax = pred_fn),
        inherit.aes = FALSE,
        color = "darkred",
        width = 0.48,
        linewidth = 0.75
      ) +
      scale_color_viridis_d() +
      coord_cartesian(ylim = c(0, 1)) +
      facet_wrap(~ panel_label, ncol = ncol) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            legend.position = "none") +
      xlab(x_title) +
      ylab("Fraction high (points = samples, red bar = model)")

    return(gg)
  }

  ### Layout B: Continuous time-course trajectories + sample points & 95% CIs
  curves_df <- predict_gof_trajectories(
    obj = obj,
    gof_name = gof_name,
    selected_keys_df = sel_gof[, f_cols, drop = FALSE],
    n_grid = as.integer(n_grid)
  )
  for(fc in f_cols) curves_df[[fc]] <- as.character(curves_df[[fc]])
  curves_df <- dplyr::inner_join(curves_df, panel_map, by = f_cols)

  if(identical(gof_meta[["source_type"]], "dynamics")){
    sf_names <- names(gof_meta[["modeled_to_measured"]])
    if(!is.null(sf_real) && length(sf_names) > 1L){
      obs_sub$curve_group <- paste0(obs_sub[[sf_real]], ":", obs_sub$measured_specie)
      leg_title <- paste0(sf_real, ":specie")
    }else{
      obs_sub$curve_group <- obs_sub$measured_specie
      leg_title <- "Measured specie"
    }
  }else{
    if(length(cond_cols) > 0L){
      obs_sub$curve_group <- apply(obs_sub[, cond_cols, drop = FALSE], 1L, paste, collapse = ":")
      leg_title <- paste(cond_cols, collapse = ":")
    }else{
      obs_sub$curve_group <- "All samples"
      leg_title <- "Group"
    }
  }

  gg <- ggplot() +
    geom_line(
      data = curves_df,
      aes(x = tl, y = pred_fn, color = curve_group, group = curve_group),
      linewidth = line_width
    )

  if(show_ci){
    gg <- gg +
      geom_errorbar(
        data = obs_sub,
        aes(x = tl, ymin = ci_lower_fn, ymax = ci_upper_fn, color = curve_group),
        width = 0.02 * max(obs_sub$tl, na.rm = TRUE),
        alpha = 0.45
      )
  }

  gg <- gg +
    geom_point(
      data = obs_sub,
      aes(x = tl, y = obs_fn, color = curve_group),
      size = point_size,
      alpha = 0.85
    ) +
    scale_color_viridis_d() +
    coord_cartesian(ylim = c(0, 1)) +
    facet_wrap(~ panel_label, ncol = ncol) +
    theme_classic() +
    xlab("Label time (tl)") +
    ylab("Fraction high") +
    guides(color = guide_legend(title = leg_title))

  return(gg)
}

