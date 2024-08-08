#' Average parameter estimates across replicates, and regularize variance estimates
#'
#' `AverageAndRegularize` uses a linear model to average estimates of a parameter
#' of interest over replicates, and to get averages for all of a set of conditions
#' specified by the user. This specification is done through formulas that will
#' be used to create a design matrix for parameters to estimate. Currently,
#' some design matrices will yield much longer runtimes than others. This is due
#' to some design matrices necessitating numerically estimating a set of maximum
#' likelihood parameters, while others yield simple, analytic approximations that
#' are quick to calculate. For example, formulas such as `~ treatment` or
#' `~ treatment:duration` will be quick to fit, and formulas such as
#' `~ treatment + batch` are currently much slower to fit. This will eventually
#' be remedied through approximation of the solution to all full models.
#'
#' @param obj An `EZbakRFractions` object, which is an `EZbakRData` object on
#' which `EstimateFractions()` has been run.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param parameter Parameter to average across replicates of a given condition.
#' @param type What type of table is the parameter found in? Default is "kinetics",
#' but can also set to "fractions".
#' @param kstrat If `type == "kinetics"`, then `kstrat` specifies the kinetic parameter
#' inference strategy.
#' @param repeatID If multiple `kinetics` or `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param formula_mean An R formula object specifying how the `parameter` of interest
#' depends on the sample characteristics specified in `obj`'s metadf. The most common formula
#' will be `~ treatment` or `~ treatment:duration`, where `treatment` and `duration` would
#' be replaced with whatever you called the relevant sample characteristics in your metadf.
#' `~ treatment` means that an average value of `parameter` should be estimated for each
#' set of samples with the same value for `treatment` in the metadf. `~ treatment:duration` specifies
#' that an average value of `parameter` should be estimated for each set of samples with the same
#' combination of `treatment` and `duration` values in the metadf. An example of the latter
#' case is a situation where you have two or more treatments (e.g., drug treated and untreated control)
#' which were applied for different durations of time (e.g., 4 and 8 hours).
#' @param formula_sd Same as `formula_mean`, but this time specifying how the variance
#' in the replicate estimates of `parameter` depends on the sample characteristics specified
#' in `obj`'s metadf. Unlike standard linear modeling, this allows you to specify
#' a heteroskedastic model. I suggest allowing a `parameter`'s variance to depend on
#' the "treatment" condition, as changes in relative RNA abundance can impact `parameter` variance,
#' so differential expression caused by your "treatment" could impact `parameter` variance.
#' @param include_all_parameters If TRUE, an additional table will be saved with the prefix `fullfit_`,
#' which includes all of the parameters estimated throughout the course of linear modeling and
#' regularization. This can be nice for visualizing the regularized mean-variance trend.
#' @param sd_reg_factor Determines how strongly variance estimates are shrunk towards trend.
#' Higher numbers lead to more regularization. Eventually, this will be replaced with estimation
#' of how much variance there seems to be in the population of variances.
#' @param error_if_singular If TRUE, linear model will throw an error if parameters
#' cannot be uniquely identified. This is most often caused by parameters that cannot
#' be estimated from the data, e.g., due to limited replicate numbers or correlated
#' sample characteristics (i.e., all treatment As also correspond to batch As, and
#' all treatment Bs correspond to batch Bs).
#' @param min_reads Minimum number of reads in all samples for a feature to be kept.
#' @param force_optim Perform numerical likelihood estimation, even if a more efficient,
#' approximate, analytical strategy is possible given `formula_mean` and `formula_sd`.
#' @param conservative If TRUE, conservative variance regularation will be performed.
#' In this case, variances below the trend will be regularized up to the trend, and
#' variances above the trend will be left unregularized. This avoids underestimation
#' of variances.
#' @param character_limit Limit on the number of characters of the name given to the
#' output table. Will attempt to concatenate the parameter name with the names of all
#' of the features. If this is too long, only the parameter name will be used.
#' @param feature_lengths Table of effective lengths for each feature combination in your
#' data. For example, if your analysis includes features named GF and XF, this
#' should be a data frame with columns GF, XF, and length.
#' @param overwrite If TRUE, identical, existing output will be overwritten.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
AverageAndRegularize <- function(obj, features = NULL, parameter = "log_kdeg",
                                 type = "kinetics",
                                 kstrat = NULL,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 exactMatch = TRUE,
                                 repeatID = NULL,
                                 formula_mean = NULL, formula_sd = NULL,
                                 include_all_parameters = TRUE,
                                 sd_reg_factor = 10,
                                 error_if_singular = TRUE,
                                 min_reads = 10,
                                 force_optim = FALSE,
                                 conservative = FALSE,
                                 character_limit = 20,
                                 feature_lengths = NULL,
                                 overwrite = TRUE){



  if(parameter == "log_TILAC_ratio"){

    TILAC <- TRUE

  }

  obj <- general_avg_and_reg(obj = obj, features = features, parameter = parameter,
                             type = type, kstrat = kstrat, populations = populations,
                             exactMatch = exactMatch,
                             repeatID = repeatID,
                             fraction_design = fraction_design,
                             formula_mean = formula_mean, formula_sd = formula_sd,
                             include_all_parameters = include_all_parameters,
                             sd_reg_factor = sd_reg_factor,
                             error_if_singular = error_if_singular, min_reads = min_reads,
                             force_optim = force_optim,
                             conservative = conservative,
                             TILAC = TILAC,
                             character_limit = character_limit,
                             feature_lengths = feature_lengths,
                             overwrite = overwrite)

  return(obj)

}


# See if a factor referenced in a formula object has only a single factor.
# This will break the linear modeling strategies I employ in this function
checkSingleLevelFactors <- function(formula, data) {

  variables <- all.vars(formula)

  # Initialize a vector to keep track of single-level factors
  singleLevelFactors <- logical(length = 0)

  # Loop through each variable to check if it's a factor present
  # in the data with only one level
  for (var in variables) {

    if (var %in% names(data) && is.factor(data[[var]])) {

      if (length(levels(data[[var]])) == 1) {

        singleLevelFactors <- c(singleLevelFactors, TRUE)

      } else {

        singleLevelFactors <- c(singleLevelFactors, FALSE)

      }
    }
  }

  # Return TRUE if any single-level factors are found, FALSE otherwise
  any(singleLevelFactors)
}


# Function to get normal distribution posterior mean
# This is log-scale standard deviation regularization
# Currently have implemented a conservative strategy where
# standard deviations less than the trend are regularized
# all the way up to the trend, and standard deviations above the
# trend are not regularized at all. Similar to the old
# DESeq model (pre-DESeq2).
get_sd_posterior <- function(n = 1, sd_est, sd_var,
                             fit_var, fit_mean,
                             conservative = FALSE){


  denom <- (n/sd_var + 1/fit_var)
  num <- sd_est/sd_var + fit_mean/fit_var

  if(conservative){

    output <- ifelse(sd_est > fit_mean,
                     sd_est,
                     fit_mean)

  }else{

    output <- ifelse(sd_est > fit_mean,
                     num/denom,
                     fit_mean)
  }

  return(output)

}



general_avg_and_reg <- function(obj, features, parameter,
                                type = "kinetics",
                                kstrat = NULL,
                                populations = NULL,
                                exactMatch = NULL,
                                repeatID = NULL,
                                fraction_design = NULL,
                                formula_mean, formula_sd,
                                include_all_parameters,
                                sd_reg_factor,
                                error_if_singular,
                                TILAC = FALSE, min_reads = 10,
                                force_optim = FALSE,
                                conservative = FALSE,
                                character_limit = 20,
                                feature_lengths = NULL,
                                overwrite = TRUE){


  ### Hack to deal with devtools::check() NOTEs
  n <- log_normalized_reads <- logsd <- logsd_from_uncert <- replicates <- NULL
  coverage <- se_mean <- se_logsd <- coverages <- parameters <- normalized_reads <- NULL

  `.` <- list


  ### Get name of standard error column

  parameter_se <- paste0("se_", parameter)


  ### Extract kinetic parameters of interest

  # metadf for covariates
  metadf <- obj$metadf

  ### Figure out which table to use


  # Function is in Helpers.R
  param_name <- EZget(obj,
                      type = type,
                      features = features,
                      kstrat = kstrat,
                      exactMatch = exactMatch,
                      populations = populations,
                      fraction_design = fraction_design,
                      repeatID = repeatID,
                      returnNameOnly = TRUE)


  # Get fractions
  kinetics <- obj[[type]][[param_name]]

  features_to_analyze <- obj[["metadata"]][[type]][[param_name]][["features"]]


  if(type == "fractions"){

    normalized_reads <- get_normalized_read_counts(obj = obj,
                                                   features_to_analyze = features_to_analyze,
                                                   fractions_name = param_name,
                                                   feature_lengths = feature_lengths) %>%
      dplyr::as_tibble()

    # Get the kinetic parameter data frame
    kinetics <- kinetics %>%
      dplyr::inner_join(normalized_reads[,c("sample",
                                            features_to_analyze,
                                            "normalized_reads")],
                        by = c("sample", features_to_analyze)) %>%
      dplyr::mutate(log_normalized_reads = log10(normalized_reads))


    # Which samples need to get filtered out
    mcols <- colnames(metadf)
    tl_cols <- mcols[grepl("^tl", mcols) | grepl("^tpulse", mcols)]

    samples_with_no_label <- metadf %>%
      dplyr::rowwise() %>%
      dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) == 0)) %>%
      dplyr::ungroup() %>%
      dplyr::select(sample) %>%
      unlist() %>%
      unname()

    # Filter out label-less samples
    kinetics <- kinetics %>%
      dplyr::filter(!(sample %in% samples_with_no_label))

  }else{

    # Get the kinetic parameter data frame
    kinetics <- kinetics %>%
      dplyr::mutate(log_normalized_reads = log10(normalized_reads))


  }



  if(is.null(formula_mean)){

    # Add kinetic parameter column to formula\
    condition_vars <- colnames(metadf)[!grepl("tl", colnames(metadf)) &
                                         (colnames(metadf) != "sample")]


    formula_mean <- stats::as.formula(paste0("~", paste(condition_vars, collapse = "+")))

  }


  formula_mean <- stats::as.formula(paste0(paste(c(parameter, formula_mean), collapse = ""), "-1"))


  if(is.null(formula_sd)){

    # Add kinetic parameter column to formula\
    condition_vars <- colnames(metadf)[!grepl("tl", colnames(metadf)) &
                                         (colnames(metadf) != "sample")]


    formula_sd <- stats::as.formula(paste0("~", paste(condition_vars, collapse = "+")))

  }

  formula_reads <- stats::as.formula(paste0(paste(c("log_normalized_reads", formula_sd), collapse = ""), "-1"))


  formula_sd <- stats::as.formula(paste0(paste(c(parameter, formula_sd), collapse = ""), "-1"))



  ### Fit linear, potentially heteroskedastic model
  meta_cols <- colnames(metadf)
  tl_cols <- meta_cols[grepl("^tl", meta_cols)]

  # Add covariates to kinetics
  kinetics <- kinetics %>%
    dplyr::inner_join(metadf, by = "sample")



  ### Filter out features that are not present in all samples

  num_samps <- length(unique(kinetics$sample))

  features_to_keep <- kinetics %>%
    dplyr::filter(n > min_reads) %>%
    dplyr::select(-n) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(features_to_analyze))) %>%
    dplyr::count() %>%
    dplyr::filter(n == num_samps) %>%
    dplyr::select(!!features_to_analyze)


  kinetics <- kinetics %>%
    dplyr::inner_join(features_to_keep,
                      by = features_to_analyze)


  ### Need to cover for edge case where there are single level factors
  single_level_mean <- checkSingleLevelFactors(kinetics,
                                               formula_mean)
  single_level_sd <- checkSingleLevelFactors(kinetics,
                                             formula_sd)



  message("Fitting linear model")
  if(single_level_mean | single_level_sd){

    if(length(all.vars(formula_mean)) != 2 | length(all.vars(formula_sd)) != 2 ){

      stop("You have specified a multi-covariate model where one or more of the covariates
           has a single level. The mean and standard deviations of all of the groups
           in such a model cannot be estimated. For example, if you have 4 samples,
           and two covariates, call them treatment and batch, and treatment has the values
           of 'A', 'A', 'B', 'B' for the 4 samples, and batch has a value of 'Z' for all 4 samples,
           there is no way to estimate the average value of a parameter in batch 'Z' because
           the average value in the first two samples is the mu_A + mu_Z, and in the
           second two samples is the sum of mu_B and mu_Z, where mu_x is the average
           value of the parameter of interest in group x. In this case, there are
           three parameters to estimate (mu_A, mu_B, and mu_Z) but only two
           unique groups of data points (A+Z and B+Z).
           SOLUTION: remove the single-level covariate(s) or specify a model
           in which there is only a single covariate, a single-level one.")


    }

    ### TO-DO; DEAL WITH THIS EDGE-CASE EFFECTIVELY; IMPORTANT FOR TILAC
    model_fit <- kinetics %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(features_to_analyze))) %>%
      dplyr::summarise()

  }else{

    ### If there is only one element of formula, then just do simple averaging.
    ### Else, will have to estimate maximum likelihood parameters
    mean_vars <- all.vars(formula_mean)
    sd_vars <- all.vars(formula_sd)

    if(length(mean_vars) == 2 & length(sd_vars) == 2 & !force_optim){

      # It's much faster this way
      model_fit <- kinetics %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(mean_vars[2], features_to_analyze)))) %>%
        dplyr::summarise(mean = mean(!!dplyr::sym(parameter)),
                         logsd = log(stats::sd(!!dplyr::sym(parameter)) ),
                         logsd_from_uncert = log(mean(!!dplyr::sym(parameter_se))),
                         coverage = mean(log_normalized_reads),
                         replicates = dplyr::n()) %>%
        dplyr::mutate(logsd = log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2) + 0.025/sqrt(replicates))) %>%
        # dplyr::mutate(logsd = dplyr::case_when(
        #   logsd < logsd_from_uncert ~ log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2)),
        #   #logsd < logsd_from_uncert ~ logsd_from_uncert,
        #   .default = logsd
        # )) %>%
        dplyr::select(-logsd_from_uncert) %>%
        dplyr::mutate(se_mean = exp(logsd)/sqrt(replicates),
                      se_logsd = 1/sqrt(2*(replicates - 1)),
                      logsd = log(exp(logsd)/sqrt(replicates))) %>%
        dplyr::select(-replicates) %>%
        tidyr::pivot_wider(names_from = !!mean_vars[2],
                           values_from = c(mean, logsd, coverage, se_mean, se_logsd),
                           names_sep = paste0("_", mean_vars[2]))


    }else if(length(mean_vars) == 2 & length(sd_vars) == 1 & !force_optim){

      # It's faster this way
      # Bit tricky as you have to:
      # 1) Estimate mean of sd in each group
      # 2) Average sd across groups
      # Can't just do sd across all groups, as this will be a function of both
      # replicate variability and difference across groups. So will overestimate sd.
      # 3) Calculate se of sd based on number of replicates across all groups
      # 4) Calculate se of mean based on number of replicates within a group
      model_fit <- kinetics %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(mean_vars[2], features_to_analyze)))) %>%
        dplyr::summarise(mean = mean(!!dplyr::sym(parameter)),
                         logsd = log(stats::sd(!!dplyr::sym(parameter)) ), # Lower MSE estimator
                         coverage = mean(coverage),
                         se_logsd = mean(se_logsd)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(features_to_analyze)))) %>%
        dplyr::mutate(logsd = mean(logsd),
                      logsd_from_uncert = log(mean(!!dplyr::sym(parameter_se))),
                      coverage = mean(log_normalized_reads)) %>%
        dplyr::mutate(logsd = log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2) + 0.025/sqrt(dplyr::n()))) %>%
        # dplyr::mutate(logsd = dplyr::case_when(
        #   logsd < logsd_from_uncert ~ log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2)),
        #   #logsd < logsd_from_uncert ~ logsd_from_uncert,
        #   .default = logsd
        # )) %>%
        dplyr::select(-logsd_from_uncert) %>%
        dplyr::mutate(se_logsd = 1/sqrt(2*(dplyr::n() - 1))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(mean_vars[2], features_to_analyze)))) %>%
        dplyr::mutate(se_mean = exp(logsd)/sqrt(dplyr::n()),
                      logsd = log(exp(logsd)/sqrt(dplyr::n()))) %>%
        tidyr::pivot_wider(names_from = !!mean_vars[2],
                           values_from = c(mean, logsd, coverage, se_mean, se_logsd),
                           names_sep = paste0("_", mean_vars[2]))

    }else if(interaction_only(formula_mean) & interaction_only(formula_sd) & length(mean_vars) == length(sd_vars) & !force_optim){

      ### CURRENTLY A SLIGHT OVERSIMPLIFICATION THAT ASSUMES A USER WILL HAVE THE
      ### SAME MEAN AND SD FORMULAS

      # It's much faster this way
      model_fit <- kinetics %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(mean_vars[2:length(mean_vars)], features_to_analyze)))) %>%
        dplyr::summarise(mean = mean(!!dplyr::sym(parameter)),
                         logsd = log(stats::sd(!!dplyr::sym(parameter)) ),
                         logsd_from_uncert = log(mean(!!dplyr::sym(parameter_se))),
                         coverage = mean(log_normalized_reads),
                         replicates = dplyr::n()) %>%
        dplyr::mutate(logsd = log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2) + 0.025/sqrt(replicates))) %>%
        # dplyr::mutate(logsd = dplyr::case_when(
        #   logsd < logsd_from_uncert ~ log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2)),
        #   #logsd < logsd_from_uncert ~ logsd_from_uncert,
        #   .default = logsd
        # )) %>%
        dplyr::select(-logsd_from_uncert) %>%
        dplyr::mutate(logsd = log(exp(logsd)/sqrt(replicates)),
                      se_mean = exp(logsd)/sqrt(replicates),
                      se_logsd = 1/sqrt(2*(replicates - 1)) ) %>%
        dplyr::select(-replicates) %>%
        tidyr::pivot_wider(names_from = !!mean_vars[2:length(mean_vars)],
                           values_from = c(mean, logsd, coverage, se_mean, se_logsd),
                           names_glue = paste0("{.value}_", paste(paste0(mean_vars[2:length(mean_vars)],
                                                                         "{", mean_vars[2:length(mean_vars)],"}"),
                                                                  collapse = ":") ))

    }
    else{

      model_fit <- kinetics %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(features_to_analyze))) %>%
        dplyr::do(dplyr::tibble(parameters = list(fit_heteroskedastic_linear_model(.,
                                                                                   formula_mean = formula_mean,
                                                                                   formula_sd = formula_sd,
                                                                                   dependent_var = parameter,
                                                                                   error_if_singular = error_if_singular)),
                                coverages = list(calc_avg_coverage(.,
                                                                   formula = formula_reads))))

      model_fit <- model_fit %>%
        tidyr::unnest_wider(coverages) %>%
        tidyr::unnest_wider(parameters)

      # Add carriage return so next message is separate from progress bar
      message("")

    }






  }



  ### Estimate coverage vs. variance trends for each standard deviation estimate

  message("Estimating coverage vs. variance trend")


  # Step 1: Filter column names for relevant patterns
  sd_columns <- names(model_fit)[grepl("^logsd_", names(model_fit))]
  covariate_names <- substring(sd_columns, 7)
  coverage_columns <- names(model_fit)[grepl("^coverage_", names(model_fit))]
  relevant_columns <- union(sd_columns, coverage_columns)

  # Step 2: Iterate and perform regression
  regression_results <- purrr::map(covariate_names, ~ {

    # If modeling fraction news, then mean absolute value of fraction new should
    # be included in the regression. In general, this is something I would like
    # to include, but that is difficult to include if the fraction new information
    # is not present in the object being modeled. I could go back and get it, but
    # that is likely a large refactor
    if(type == "fractions"){

      formula_str <- paste("`logsd_", .x, "`", " ~ `coverage_", .x, "`",
                           " + abs(`mean_", .x, "`)", sep = "")

    }else{

      formula_str <- paste("`logsd_", .x, "`", " ~ `coverage_", .x, "`", sep = "")

    }

    formula <- stats::as.formula(formula_str)

    # Perform linear regression
    lm_result <- stats::lm(formula, data = model_fit)

    # Return result
    return(list(covariate = .x, lm_result = lm_result))
  })

  names(regression_results) <- covariate_names


  # Iterate over the regression results to add predicted values to the dataframe
  for(result in regression_results) {
    # Extract the covariate name
    covariate <- result$covariate

    # Generate the fitted values using the predict function
    fitted_values <- stats::predict(result$lm_result)

    # Create the new column name
    new_column_name <- paste("logsd_", covariate, "_fit", sep = "")

    # Add the fitted values to the dataframe as a new column
    model_fit[[new_column_name]] <- fitted_values
  }


  ### Regularize variance estimates

  message("Regularizing variance estimates")


  # Loop over each covariate's sd to be regularized
  cols_to_keep <- c()
  coverage_cols <- c()
  for(c in covariate_names){

    # Names of columns to pass to regularization function
    col_name <- paste0("logsd_", c, "_posterior")
    natural_col_name <- paste0("sd_", c, "_posterior")
    sd_est_name <- paste0("logsd_", c)
    sd_var_name <- paste0("se_logsd_", c)
    fit_mean_name <- paste0("logsd_", c, "_fit")
    sd_mean_name <- paste0("se_mean_", c)


    # Regularize
    model_fit <- model_fit %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!col_name := get_sd_posterior(sd_est = !!dplyr::sym(sd_est_name),
                                                   sd_var = (!!dplyr::sym(sd_var_name)) ^ 2,
                                                   fit_var = stats::var(regression_results[[c]]$lm_result$residuals) / sd_reg_factor,
                                                   fit_mean = !!dplyr::sym(fit_mean_name),
                                                   conservative = conservative)) %>%
      dplyr::mutate(!!natural_col_name := exp((!!dplyr::sym(col_name))))

    mean_est <- paste0("mean_", c)
    coverage_est <- paste0("coverage_", c)

    cols_to_keep <- c(cols_to_keep, mean_est, natural_col_name)
    coverage_cols <- c(coverage_cols, coverage_est)
  }


  # Want to stick the coverage columns at the end
  final_output <- model_fit %>%
    dplyr::select(!!features_to_analyze, !!cols_to_keep, !!coverage_cols)



  # What should output be named?
  avg_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")
  avg_vect <- paste0(gsub("_","", parameter), "_", avg_vect)

  if(nchar(avg_vect) > character_limit){

    psearch <- paste0("^", gsub("_", "", parameter), "[1-9]")

    num_avgs <- sum(grepl(psearch, names(obj[['averages']])))
    avg_vect <- paste0("averages", num_avgs + 1)

  }

  # Are there any metadata or fractions objects at this point?
  if(length(obj[['metadata']][['averages']]) > 0){

    avg_vect <- decide_output(obj,
                              proposed_name = avg_vect,
                              type = "averages",
                              features = features_to_analyze,
                              parameter = parameter,
                              mean_vars = mean_vars[2:length(mean_vars)], # 1st element == parameter
                              sd_vars = sd_vars[2:length(sd_vars)], # 1st element == parameter
                              overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'averages',
                               features = features_to_analyze,
                               parameter = parameter,
                               returnNameOnly = TRUE,
                               mean_vars = mean_vars[2:length(mean_vars)], # 1st element == parameter
                               sd_vars = sd_vars[2:length(sd_vars)], # 1st element == parameter
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }


  # Save
  obj[['averages']][[avg_vect]] <- dplyr::as_tibble(final_output)

  # Save metadata
  obj[['metadata']][['averages']][[avg_vect]] <- list(features = features_to_analyze,
                                                      parameter = parameter,
                                                      mean_vars = mean_vars[2:length(mean_vars)], # 1st element == parameter
                                                      sd_vars = sd_vars[2:length(sd_vars)], # 1st element == parameter
                                                      repeatID = repeatID)


  if(include_all_parameters){

    output_name <- paste0("fullfit_", avg_vect)
    obj[['averages']][[output_name]] <- dplyr::as_tibble(model_fit)

  }

  if(!methods::is(obj, "EZbakRFit")){

    class(obj) <- c( "EZbakRFit", class(obj))

  }

  return(obj)


}


#' Get contrasts of estimated parameters
#'
#' @param obj An `EZbakRFit` object, which is an `EZbakRFractions` object on
#' which `AverageAndRegularize()` has been run.
#' @param condition Name of factor from `metadf` whose parameter estimates at
#' different factor values you would like to compare.
#' @param reference Name of reference `condition` factor level value. Difference
#' will be calculated as `experimental` - `reference`.
#' @param experimental Name of `condition` factor level value to compare to reference.
#' Difference will be calculated as `experimental` - `reference`.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param mean_vars Sample features from metadf that were used in formula for
#' parameter average linear model.
#' @param sd_vars Sample features from metadf that were used in formula for
#' parameter standard deviation linear model.
#' @param repeatID If multiple `averages` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param parameter Parameter to average across replicates of a given condition.
#' @param overwrite If TRUE, then identical output will be overwritten if it exists.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
CompareParameters <- function(obj, condition, reference, experimental,
                              features = NULL, exactMatch = TRUE, repeatID = NULL,
                              mean_vars = NULL, sd_vars = NULL,
                              overwrite = TRUE, parameter = "log_kdeg"){


  ### Hack to deal with annoying devtools::check() NOTE

  difference <- uncertainty <- pval <- padj <- avg_coverage <- NULL


  ### Extract kinetic parameters of interest

  # metadf for covariates
  metadf <- obj$metadf


  # Function is in Helpers.R
  averages_name <- EZget(obj,
                         type = "averages",
                         features = features,
                         parameter = parameter,
                         mean_vars = mean_vars,
                         sd_vars = sd_vars,
                         exactMatch = exactMatch,
                         repeatID = repeatID,
                         returnNameOnly = TRUE)


  # Get fractions
  parameter_est <- obj[["averages"]][[averages_name]]

  features_to_analyze <- obj[["metadata"]][["averages"]][[averages_name]][["features"]]


  ### Perform comparative analysis of interest

  ref_mean <- paste0("mean_", condition, reference)
  ref_sd <- paste0("sd_", condition, reference, "_posterior")
  ref_cov <- paste0("coverage_", condition, reference)

  exp_mean <- paste0("mean_", condition, experimental)
  exp_sd <- paste0("sd_", condition, experimental, "_posterior")
  exp_cov <- paste0("coverage_", condition, experimental)

  comparison <- parameter_est %>%
    dplyr::mutate(difference = !!dplyr::sym(exp_mean) - !!dplyr::sym(ref_mean),
                  uncertainty = sqrt( (!!dplyr::sym(ref_sd))^2 + (!!dplyr::sym(exp_sd))^2 ),
                  stat = difference/uncertainty,
                  pval = 2*stats::pnorm(-abs(stat)),
                  avg_coverage = ((!!dplyr::sym(ref_cov)) + (!!dplyr::sym(exp_cov))) / 2) %>%
    dplyr::mutate(padj = stats::p.adjust(pval, method = "BH")) %>%
    dplyr::select(!!features_to_analyze, difference, uncertainty, stat, pval, padj, avg_coverage)



  ### Add output to object

  num_comparisons <- length(obj[['comparisons']])

  output_name <- paste0("comparison", num_comparisons + 1)


  if(length(obj[['comparisons']]) > 0){
    output_name <- decide_output(obj, output_name, type = "comparisons",
                                 features = features, parameter = parameter,
                                 reference = reference, experimental = experimental,
                                 condition = condition,
                                 overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'comparisons',
                               features = features,
                               parameter = parameter,
                               condition = condition,
                               reference = reference,
                               experimental = experimental,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }

  obj[["comparisons"]][[output_name]] <- dplyr::as_tibble(comparison)

  obj[["metadata"]][["comparisons"]][[output_name]] <- list(condition = condition,
                                                            reference = reference,
                                                            experimental = experimental,
                                                            features = features_to_analyze,
                                                            parameter = parameter,
                                                            repeatID = repeatID)


  if(!methods::is(obj, "EZbakRCompare")){

    class(obj) <- c( "EZbakRCompare", class(obj))

  }

  return(obj)


}


# Define the likelihood function
heteroskedastic_likelihood <- function(params, y, X_mean, X_sd, debug = FALSE) {

  if(debug){
    browser()
  }

  n <- length(y)
  beta <- params[1:ncol(X_mean)]
  log_sigma <- params[(ncol(X_mean) + 1):length(params)]

  mu <- X_mean %*% beta
  sigma <- exp(X_sd %*% log_sigma)

  # Gaussian log-likelihood
  logL <- -sum(stats::dnorm(y, mean=mu, sd=sigma, log=TRUE))
  return(logL)
}


fit_heteroskedastic_linear_model <- function(formula_mean, formula_sd, data,
                                             dependent_var,
                                             error_if_singular = TRUE) {


  # Parse formula objects into model matrices
  designMatrix_mean <- stats::model.matrix(formula_mean, data)
  designMatrix_sd <- stats::model.matrix(formula_sd, data)


  if(error_if_singular){

    # Check for singularity in the design matrices
    if(qr(designMatrix_mean)$rank < ncol(designMatrix_mean)) {
      stop("Design matrix for mean is singular")
    }
    if(qr(designMatrix_sd)$rank < ncol(designMatrix_sd)) {
      stop("Design matrix for standard deviation is singular")
    }

  }


  # Initial parameter guesses; keeping it realistic to hopefully increase odds of convergence
  startParams <- c(rep(mean(data[[dependent_var]]), times = ncol(designMatrix_mean)),
                   rep(log(stats::sd(data[[dependent_var]]) + 0.001), times = ncol(designMatrix_sd)))

  # Try BFGS
  opt <- stats::optim(startParams, heteroskedastic_likelihood,
                      #gr = heteroskedastic_gradient,
                      y = data[[dependent_var]],
                      X_mean = designMatrix_mean,
                      X_sd = designMatrix_sd,
                      method = "BFGS",
                      hessian = TRUE)



  if(opt$convergence != 0) {



    # Sometimes L-BFGS-B converges when BFGS doesn't, just out of pure luck
    # of it being a rougher approximation I guess.
    # The reason I try BFGS first is because I found that it was far more
    # likely to converge on any given data subset. Kinda surprised by how often
    # L-BFGS-B doesn't converge, but these are some quasi-non-trivial models I am
    # fitting.
    # I will note that "failed convergence" almost always meant that the
    # max number of iterations was hit (code 1).
    opt <- stats::optim(startParams, heteroskedastic_likelihood,
                        #\gr = heteroskedastic_gradient,
                        y = data[[dependent_var]],
                        X_mean = designMatrix_mean,
                        X_sd = designMatrix_sd,
                        method = "L-BFGS-B",
                        upper = rep(7, times = length(startParams)),
                        lower = rep(-7, times = length(startParams)),
                        hessian = TRUE)


    if(opt$convergence != 0){

      warning("Model did not converge for one feature set!")


    }

  }

  if(error_if_singular){

    # Singularity check using the determinant of the Hessian matrix
    svd_hessian <- svd(opt$hessian)
    if(min(svd_hessian$d) < 1e-10){
      stop("Model parameters are unidentifiable! Specify a simpler
         design matrix.")
    }

  }

  tryCatch({
    solve(opt$hessian)
  }, error = function(e) {
    print(e)
    print("Looks like you don't have enough data to estimate all of your parameters!
          Specify a simpler model and try again.")
  })

  # Estimate uncertainty from Hessian
  ses <- sqrt(diag(solve(opt$hessian)))

  estimates <- c(opt$par, ses)
  names(estimates) <- c(paste0("mean_", colnames(designMatrix_mean)),
                        paste0("logsd_", colnames(designMatrix_sd)),
                        paste0("se_mean_", colnames(designMatrix_mean)),
                        paste0("se_logsd_", colnames(designMatrix_sd)))

  return(as.list(estimates))
}

calc_avg_coverage <- function(data, formula){

  # Calculate average read counts in each group
  fit <- stats::lm(formula = formula,
                   data = data)


  estimates <- fit$coefficient
  names(estimates) <- c(paste0("coverage_", names(estimates)))

  return(as.list(estimates))

}


# Are there only interaction terms in a formula?
interaction_only <- function(formula) {

  # Extract terms information
  terms_info <- stats::terms(formula)

  # Convert to character for easier handling
  terms_char <- attr(terms_info, "term.labels")

  # Check if all terms are interaction terms or '-1'
  all_interactions_or_no_intercept <- all(sapply(terms_char, function(term) {
    grepl("[:*]", term) || term == "-1"
  }))

  return(all_interactions_or_no_intercept)
}


# Heteroskedastic regression log-likelihood gradient
heteroskedastic_gradient <- function(params, y, X_mean, X_sd){


  n <- length(y)
  beta <- params[1:ncol(X_mean)]
  log_sigma <- params[(ncol(X_mean) + 1):length(params)]

  mu <- X_mean %*% beta
  sigma <- exp(X_sd %*% log_sigma)

  # Vectors to be Hadamard multiplied with desigm matrices
  b_vects <- ((y - mu)/(sigma^2))
  s_vects <- ((y - mu)^2)/(sigma^3)

  # Hadamard product for betas
  # Because derivative is 0 if design matrix is 0 and derivative of
  # normal log-likelihood otherwise
  nr <- nrow(X_mean)
  ncb <- ncol(X_mean)
  ncs <- ncol(X_sd)

  dprod_beta <- matrix(0,ncol=ncb, nrow=nr)
  dprod_sigma <- matrix(0,ncol=ncs, nrow=nr)
  for(i in 1:nr){

    dprod_beta[i,] <- X_mean[i,1:ncb]*b_vects[i]
    dprod_sigma[i,] <- X_sd[i,1:ncs]*s_vects[i]*sigma[i]

  }


  # Sum derivatives of each data point to get total log-likelihood derivatives
  dlldlb <- colSums(dprod_beta)
  dlldls <- colSums(dprod_sigma)


  return(c(-dlldlb, -dlldls))



}


