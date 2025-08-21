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

  if(is.null(comparison_name)){
    stop("No tables from your EZbakR analysis match your search criteria!")
  }

  comparison <- obj[['comparisons']][[comparison_name]]
  metadata <- obj[['metadata']][['comparisons']][[comparison_name]]


  # Track which features user wants highlighted
  if(!is.null(features_to_highlight)){

    features <- obj[['metadata']][['comparisons']][['features']]

    if(is.null(features)){
      # Mainly for backwards compatibility where a bug prevented the features
      # entry of the comparisons metadata from getting updated. Just try to
      # join by the features that exist and hope for the best.

      comparison <- comparison %>%
        dplyr::left_join(
          features_to_highlight %>%
            dplyr::mutate(
              highlight = TRUE
            )
        ) %>%
        tidyr::replace_na(
          list(
            highlight = FALSE
          )
        )

    }else{

      if(is.vector(features_to_highlight)){
        # Make into one column data frame is user provided a single vector of
        # feature names

        if(length(features) > 1){

          stop("You provided a vector to features_to_highlight but your comparison table has more than one feature!
              features_to_highlight must be a data frame with columns for each of the feature
              types tracked in the relevant comparison table.")

        }

        features_to_highlight <- dplyr::tibble(
          !!features := features_to_highlight
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

    features <- obj[['metadata']][['comparisons']][['features']]

    if(is.null(features)){
      # Mainly for backwards compatibility where a bug prevented the features
      # entry of the comparisons metadata from getting updated. Just try to
      # join by the features that exist and hope for the best.

      comparison <- comparison %>%
        dplyr::left_join(
          features_to_highlight %>%
            dplyr::mutate(
              highlight = TRUE
            )
        ) %>%
        tidyr::replace_na(
          list(
            highlight = FALSE
          )
        )

    }else{

      if(is.vector(features_to_highlight)){
        # Make into one column data frame is user provided a single vector of
        # feature names

        if(length(features) > 1){

          stop("You provided a vector to features_to_highlight but your comparison table has more than one feature!
              features_to_highlight must be a data frame with columns for each of the feature
              types tracked in the relevant comparison table.")

        }

        features_to_highlight <- dplyr::tibble(
          !!features := features_to_highlight
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


  if(metadf_cols_to_use == "all"){
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
      variance = var(!!dplyr::sym(data_col))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      variance > quantile(variance, probs = seq(0, 1, 0.1))[variance_decile]
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
  fn_pca_df <- dplyr::tibble(PC1 = fn_PC1,
                      PC2 = fn_PC2) %>%
    dplyr::bind_cols(obj$metadf %>%
                       dplyr::filter(
                         sample %in% samps_to_keep
                       ) %>%
                       dplyr::select(
                         !!metadf_cols_to_use
                       )
                ) %>%
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
