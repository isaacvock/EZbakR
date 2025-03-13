#' Average parameter estimates across replicates, and regularize variance estimates
#'
#' `AverageAndRegularize` fits a generalized linear model to your data to effectively
#' average parameter estimates across replicates and get overall uncertainty estimates
#' for those parameters. The linear model to which your data is fit is specified via
#' an R formula object supplied to the `formula_mean` parameter. Uncertainty estimates
#' are regularized via a hierarchical modeling strategy originally introduced with
#' bakR, though slightly improved upon since then.
#'
#' The EZbakR website has an extensive vignette walking through various use cases
#' and model types that you can fit with `AverageAndRegularize()`: [vignette link](https://isaacvock.github.io/EZbakR/articles/Linear-modeling.html).
#' EZbakR improves upon bakR by balancing extra conservativeness in several steps
#' with a more highly powered statistical testing scheme in its `CompareParameters()`
#' function. In particular, the following changes to the variance regularization
#' scheme were made:
#' \itemize{
#'  \item Sample-specific parameter uncertainties are used to generate conservative estimates
#'  of feature-specific replicate variabilties. In addition, a small floor is set to ensure
#'  that replicate variance estimates are never below a certain level, for the same reason.
#'  \item Condition-wide replicate variabilities are regressed against both read coverage and
#'  either a) |logit(estimate)| when modeling average fraction labeled. This captures
#'  the fact thta estimates are best around a logit(fraction labeled) of 0 and get
#'  worse for more extreme fraction labeled's.; b) log(kdeg) when modeling log degradation
#'  rate constants. At first, I considered a strategy similar to the fraction labeled
#'  modeling, but found that agreement between a fully rigorous MCMC sampling approach
#'  and EZbakR was significantly improved by just regressing hee value of the log kientic parameter,
#'  likely due to the non-linear transformation of fraction labeled to log(kdeg);
#'  and c) only coverage in all other cases.
#'  \item Features with replicate variabilities below the inferred trend have their replicate
#'  variabilites set equal to that predicted by the trend. This helps limit underestimation
#'  of parameter variance. Features with above-trend replicate variabilties have their
#'  replicate variabilities regularized with a Normal prior Normal likelihood Bayesian
#'  model, as in bakR (so the log(variance) is the inferred mean of this distribution, and
#'  the known variance is inferred from the amount of variance about the linear dataset-wide
#'  trend).
#' }
#' All of this allows `CompareParameters()` to use a less conservative statistical test
#' when calculating p-values, while still controlling false discovery rates.
#'
#' @param obj An `EZbakRFractions` or `EZbakRKinetics` object, which is an `EZbakRData` object on
#' which `EstimateFractions()` or `EstimateKinetics()` has been run.
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
#' @param populations Character vector of the set of mutational populations
#' that you want to infer the fractions of. Only relevant if type == "fractions".
#' @param fraction_design "Design matrix" specifying which RNA populations exist
#' in your samples. Only relevant if type == "fractions".
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
#'
#' NOTE: EZbakR automatically removes any intercept terms from the model. That way,
#' there is no ambiguity about what parameter is defined as the reference.
#' @param sd_grouping_factors What metadf columns should data be grouped by when estimating
#' standard deviations across replicates? If this is NULL, then EZbakR will check to see
#' if the `formula_mean` specifies a formula that cleanly stratifies samples into disjoint
#' groups. For example, the formula `~ treatment` will assign each sample to a single factor
#' (its value for the metadf's `treatment` column). In this case, standard deviations can
#' be calculated for sets of replicates in each `treatment` group. If such a stratification
#' does not exist, a single standard deviation will be estimated for each feature (i.e.,
#' homoskedasticity will be assumed as in standard linear modeling).
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
#' @param convert_tl_to_factor If a label time variable is included in the `formula_mean`,
#' convert its values to factors so as to avoid performing continuous regression on label
#' times. Defaults to TRUE as including label time in the regression is often meant to
#' stratify samples by their label time if, for example, you are averaging logit(fractions).
#' @param regress_se_with_abs If TRUE, and if `type == "fractions"`, then standard error
#' will be regressed against logit fraction rather than magnitude of logit fraction.
#' Makes sense to set this to FALSE if analyzing certain site-specific mutational probing
#' methods when high mutation content things are likely low variance SNPs.
#' @param force_lm Certain formula lend them selves to efficient approximations of the
#' full call to `lm()`. Namely, formulas that stratify samples into disjoint groups where
#' a single parameter of the model is effectively estimated from each group can be tackled
#' via simple averaging of data from each from group. If you would like to force EZbakR
#' to fit the fully rigorous linear model though, set `force_lm` to `TRUE`.
#' @param force_optim Old parameter that is now passed the value `force_lm` and will be
#' deprecated in later releases
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
#' @return `EZbakRData` object with an additional "averages" table, as well as a
#' fullfit table under the same list heading, which includes extra information about
#' the priors used for regularization purposes.
#' @export
AverageAndRegularize <- function(obj, features = NULL, parameter = "log_kdeg",
                                 type = "kinetics",
                                 kstrat = NULL,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 exactMatch = TRUE,
                                 repeatID = NULL,
                                 formula_mean = NULL,
                                 sd_grouping_factors = NULL,
                                 include_all_parameters = TRUE,
                                 sd_reg_factor = 10,
                                 error_if_singular = TRUE,
                                 min_reads = 10,
                                 convert_tl_to_factor = TRUE,
                                 regress_se_with_abs = TRUE,
                                 force_lm = FALSE,
                                 force_optim = force_lm,
                                 conservative = FALSE,
                                 character_limit = 20,
                                 feature_lengths = NULL,
                                 overwrite = TRUE){



  if(parameter == "log_TILAC_ratio"){

    TILAC <- TRUE

  }


  ### Hack to deal with devtools::check() NOTEs
  n <- log_normalized_reads <- logse <- logse_from_uncert <- replicates <- NULL
  coverage <- se_mean <- se_logse <- coverages <- parameters <- normalized_reads <- NULL
  logsd_from_uncert <- logsd <- NULL

  `.` <- list


  ### Get name of standard error column

  parameter_se <- paste0("se_", parameter)


  ### Extract kinetic parameters of interest

  # metadf for covariates
  metadf <- obj$metadf


  # Which samples need to get filtered out
  mcols <- colnames(metadf)
  tl_cols <- mcols[grepl("^tl", mcols) | grepl("^tpulse", mcols)]


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


  ### Check to see if simple averaging is compatible with specified model
  variables <- all.vars(formula_mean)
  variables <- variables[2:length(variables)]

  # Filter out -label data
  metadf <- metadf  %>%
    dplyr::rowwise() %>%
    dplyr::filter(all(dplyr::c_across(dplyr::all_of(tl_cols)) != 0))

  # Convert tl to factor if necessary
  if(convert_tl_to_factor & any(tl_cols %in% variables)){

    cols_to_convert <- tl_cols[tl_cols %in% variables]
    metadf[cols_to_convert] <- lapply(metadf[cols_to_convert], as.factor)

  }

  can_simply_average <- TRUE
  X <- stats::model.matrix(formula_mean,
                    metadf %>%
                      dplyr::mutate(!!parameter := 1))


  if(is.null(sd_grouping_factors)){


    # If there is a clean parameter to sample mapping, we can infer
    # sd_grouping_factors
    if(all(rowSums(X!= 0) == 1)){

      sd_grouping_factors <- variables
      can_simply_average <- TRUE

    }else{

      can_simply_average <- FALSE

    }

  }else{

    # Easiest case to accommodate is if sd_grouping factors was set to all of
    # the factors in formula_mean and formula_mean permits simple averaging.
    # Could also accommodate case where sd_grouping_factors are distinct from
    # formula_mean factors, but won't worry about that for now as it would be
    # a weird decision for a user to seek out anyway.
    can_simply_average <- identical(sd_grouping_factors, variables) &
      all(rowSums(X) == 1)

  }



  ### Fit linear, potentially heteroskedastic model

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


  message("Fitting linear model")

  if(single_level_mean){

    if(length(all.vars(formula_mean)) != 2 ){

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
    ### Else, will have to run lm()
    mean_vars <- all.vars(formula_mean)
    sd_vars <- sd_grouping_factors

    # median_uncert <- median(kinetics[[parameter_se]])

    if(can_simply_average & !force_optim){


      # It's much faster this way
      model_fit <- kinetics %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(mean_vars[2:length(mean_vars)], features_to_analyze)))) %>%
        dplyr::summarise(mean = mean(!!dplyr::sym(parameter)),
                         logsd = log(stats::sd(!!dplyr::sym(parameter)) ),
                         logsd_from_uncert = log(mean(!!dplyr::sym(parameter_se))),
                         coverage = mean(log_normalized_reads),
                         replicates = dplyr::n()) %>%
        # Helps to be conservative set a floor on what you could expect the sd to be
        # That's what the `0.025/sqrt(replicates)` does
        dplyr::mutate(logsd = log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2) + 0.025/sqrt(replicates))) %>%
        # dplyr::mutate(logsd = dplyr::case_when(
        #   logsd < logsd_from_uncert & logsd_from_uncert < median_uncert ~ log(sqrt(exp(logsd_from_uncert)^2 + exp(logsd)^2) + 0.025/sqrt(replicates)),
        #   .default = log(exp(logsd) + 0.025/sqrt(replicates))
        # )) %>%
        dplyr::select(-logsd_from_uncert) %>%
        dplyr::mutate(se_mean = exp(logsd)/sqrt(replicates),
                      se_logse = 1/sqrt(2*(replicates - 1)),
                      logse = log(exp(logsd)/sqrt(replicates)))

      # Only two possibilities: you have exclusively interaction terms
      # or you have a single independent variable
      if(length(mean_vars) == 2){

        model_fit <- model_fit %>%
          dplyr::select(-replicates, -logsd) %>%
          tidyr::pivot_wider(names_from = !!mean_vars[2],
                             values_from = c(mean, logse, coverage, se_mean, se_logse),
                             names_sep = paste0("_", mean_vars[2]))


      }else{


        model_fit <- model_fit %>%
          dplyr::select(-logsd, -replicates) %>%
          tidyr::pivot_wider(names_from = !!mean_vars[2:length(mean_vars)],
                             values_from = c(mean, logse, coverage, se_mean, se_logse),
                             names_glue = paste0("{.value}_", paste(paste0(mean_vars[2:length(mean_vars)],
                                                                           "{", mean_vars[2:length(mean_vars)],"}"),
                                                                    collapse = ":") ))

      }



    }else{


      parameter_se_col <- paste0("se_", parameter)

      model_fit <- kinetics %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(features_to_analyze))) %>%
        dplyr::do(dplyr::tibble(parameters = list(
          fit_ezbakR_linear_model(.,
                                  formula_mean = formula_mean,
                                  sd_groups = sd_grouping_factors,
                                  coverage_col = "normalized_reads",
                                  uncertainties_col = parameter_se_col,
                                  error_if_singular = TRUE)
        ))) %>%
        tidyr::unnest_wider(parameters)

      # Add carriage return so next message is separate from progress bar
      message("")

    }



  }



  ### Estimate coverage vs. variance trends for each standard deviation estimate

  message("Estimating coverage vs. variance trend")


  # Step 1: Filter column names for relevant patterns
  sd_columns <- names(model_fit)[grepl("^logse_", names(model_fit))]
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

      if(regress_se_with_abs){

        formula_str <- paste("`logse_", .x, "`", " ~ `coverage_", .x, "`",
                             " + abs(`mean_", .x, "`)", sep = "")

      }else{

        formula_str <- paste("`logse_", .x, "`", " ~ `coverage_", .x, "`",
                             " + `mean_", .x, "`", sep = "")

      }


    }else if(parameter == "log_kdeg"){

      formula_str <- paste("`logse_", .x, "`", " ~ `coverage_", .x, "`",
                           " + `mean_", .x, "`", sep = "")

    }else{

      formula_str <- paste("`logse_", .x, "`", " ~ `coverage_", .x, "`", sep = "")

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
    new_column_name <- paste("logse_", covariate, "_fit", sep = "")

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
    col_name <- paste0("logse_", c, "_posterior")
    natural_col_name <- paste0("sd_", c, "_posterior")
    sd_est_name <- paste0("logse_", c)
    sd_var_name <- paste0("se_logse_", c)
    fit_mean_name <- paste0("logse_", c, "_fit")
    sd_mean_name <- paste0("se_mean_", c)


    # Regularize
    model_fit <- model_fit %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!col_name := get_sd_posterior(sd_est = !!dplyr::sym(sd_est_name),
                                                   sd_var = (!!dplyr::sym(sd_var_name)) ^ 2,
                                                   fit_var = stats::var(regression_results[[c]]$lm_result$residuals),
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

  # Strip formula of environment to ensure easy comparisons via EZget()
  environment(formula_mean) <- NULL

  # Are there any metadata or fractions objects at this point?
  if(length(obj[['metadata']][['averages']]) > 0){

    avg_vect <- decide_output(obj,
                              proposed_name = avg_vect,
                              type = "averages",
                              features = features_to_analyze,
                              parameter = parameter,
                              fit_params = covariate_names,
                              formula_mean = formula_mean,
                              kstrat = kstrat,
                              populations = populations,
                              fraction_design = fraction_design,
                              sd_grouping_factors = sd_grouping_factors,
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
                               fit_params = covariate_names,
                               formula_mean = formula_mean,
                               sd_grouping_factors = sd_grouping_factors,
                               kstrat = kstrat,
                               populations = populations,
                               fraction_design = fraction_design,
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
                                                      fit_params = covariate_names,
                                                      formula_mean = formula_mean,
                                                      sd_grouping_factors = sd_grouping_factors,
                                                      kstrat = kstrat,
                                                      populations = populations,
                                                      fraction_design = fraction_design,
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
                             fit_var_min = 0.01,
                             conservative = FALSE){


  if(fit_var <= 0){
    fit_var <- fit_var_min
  }


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



#' Get contrasts of estimated parameters
#'
#' `CompareParameters()` calculates differences in parameters estimated by
#' `AverageAndRegularize()` or `EZDynamics()` and performs null hypothesis
#' statistical testing, comparing their values to a null hypothesis of 0.
#'
#' The EZbakR website has an extensive vignette walking through various use cases
#' and parameters you can compare with `CompareParameters()`: [vignette link](https://isaacvock.github.io/EZbakR/articles/Linear-modeling.html).
#'
#' There are essentially 3 scenarios that `CompareParameters()` can handle:
#' \itemize{
#'  \item Pairwise comparisons: compare `reference` to `experimental` parameter
#'  estimates of a specified `design_factor` from `AverageAndRegularize()`. log(`experimental` / `reference`) is
#'  the computed "difference" in this case.
#'  \item Assess the value of a single parameter, which itself should represent a
#'  difference between other parameters. The name of this parameter can be specified
#'  via the `param_name` argument. This is useful for various interaction models
#'  where some of the parameters of these models may represent things like "effect of
#'  A on condition X".
#'  \item Pairwise comparison of dynamical systems model parameter estimate: similar
#'  to the first case listed above, but now when `type == "dynamics"`. `design_factor` can
#'  now be a vector of all the `metadf` columns you stratified parameter estimates by.
#' }
#' Eventually, a 4th option via the currently non-functional `param_function` argument
#' will be implemented, which will allow you to specify functions of parameters to be assessed,
#' which can be useful for certain interaction models.
#'
#' `CompareParameters()` calculates p-values using what is essentially an asymptotic Wald test,
#' meaning that a standard normal distribution is integrated. P-values are then multiple-test
#' adjusted using the method of Benjamini and Hochberg, implemented in R's `p.adjust()`
#' function.
#'
#' @param obj An `EZbakRFit` object, which is an `EZbakRFractions` object on
#' which `AverageAndRegularize()` has been run.
#' @param design_factor Name of factor from `metadf` whose parameter estimates at
#' different factor values you would like to compare. If you specify this, you need
#' to also specify `reference` and `experimental`. If `type` == "dynamics", this
#' can have multiple values, being the names of all of the factors you would like
#' to stratify a group by.
#' @param reference Name of reference `design_factor` factor level value. Difference
#' will be calculated as `experimental` - `reference`. If type == "dynamics", then this should specify the levels
#' of all of the `design_factor`(s) reference group. For example, if you have
#' multiple `design_factor`'s, then `reference` must be a named character vector with one element
#' per `design_factor`, with elements named the corresponding `design_factor`. For
#' example, if `design_factor` is c("genotype", "treatment"), and you would like
#' to compare genotype = "WT" and treatment = "untreated" (reference) to genotype = "KO" and
#' treatment = "treated", then `reference` would need to be a vector with
#' one element named "genotype", equal to "WT" and one element named "treatment"
#' equal to "untreated" (this example could be created with, `c(genotype = "WT", treatment = "untreated")`).
#' @param experimental Name of `design_factor` factor level value to compare to reference.
#' Difference will be calculated as `experimental` - `reference`. If type == "dynamics", then this should specify the levels
#' of all of the `design_factor`(s) reference group. For example, if you have
#' multiple `design_factor`'s, then `experimental` must be a named character vector with one element
#' per `design_factor`, with elements named the corresponding `design_factor`. For
#' example, if `design_factor` is c("genotype", "treatment"), and you would like
#' to compare genotype = "WT" and treatment = "untreated" (reference) to genotype = "KO" and
#' treatment = "treated", then `experimental` would need to be a vector with
#' one element named "genotype", equal to "KO" and one element named "treatment"
#' equal to "treated" (this example could be created with, `c(genotype = "KO", treatment = "treated")`).
#' @param param_name If you want to assess the significance of a single parameter,
#' rather than the comparison of two parameters, specify that one parameter's name
#' here.
#' @param parameter Parameter to average across replicates of a given condition.
#' @param type Type of table to use. Can either be "averages" or "dynamics".
#' @param param_function NOT YET IMPLEMENTED. Will allow you to specify more complicated
#' functions of parameters when hypotheses you need to test are combinations of parameters
#' rather than individual parameters or simple differences in two parameters.
#' @param condition Same as `design_factor`, will be deprecated in favor of the
#' former in later release.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param repeatID If multiple `averages` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param formula_mean An R formula object specifying how the `parameter` of interest
#' depends on the sample characteristics for the averages object you want to use.
#' @param sd_grouping_factors Metadf columns should data was grouped by when estimating
#' standard deviations across replicates for the averages object you want to use.
#' @param fit_params Character vector of parameter names in the averages object
#' you would like to use.
#' @param normalize_by_median If TRUE, then median difference will be set equal to 0.
#' This can account for global biases in parameter estimates due to things like differences
#' in effective label times. Does risk eliminating real signal though, so user discretion
#' is advised.
#' @param reference_levels Same as `reference`, but exclusively parsed in case of
#' `type` == "dynamics, included for backwards compatibility.
#' @param experimental_levels Same as `experimental`, but exclusively parsed in case of
#' `type` == "dynamics, included for backwards compatibility.
#' @param overwrite If TRUE, then identical output will be overwritten if it exists.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
CompareParameters <- function(obj, design_factor, reference, experimental,
                              param_name,
                              parameter = "log_kdeg",
                              type = "averages",
                              param_function,
                              condition = NULL,
                              features = NULL, exactMatch = TRUE, repeatID = NULL,
                              formula_mean = NULL, sd_grouping_factors = NULL,
                              fit_params = NULL, normalize_by_median = FALSE,
                              reference_levels = NULL,
                              experimental_levels = NULL,
                              overwrite = TRUE){

  ### Hack to deal with annoying devtools::check() NOTE

  difference <- uncertainty <- pval <- padj <- avg_coverage <- NULL
  exp_par <- ref_par <- exp_se <- ref_se <- NULL

  ### Determine what strategy to use for "comparisons"

  # Deal with parameter name change to make life easier for older users
  if(missing(design_factor) & !is.null(condition)){
    design_factor <- condition
  }

  if(type == "dynamics"){

    strategy <- "dynamics"

    # Catch missing arg for sake backwards compatibility
    if(missing(reference) & !is.null(reference_levels)){
      reference <- reference_levels
    }

    if(missing(experimental) & !is.null(experimental_levels)){
      experimental <- experimental_levels
    }


  }else{

    if(missing(design_factor) | missing(reference) | missing(experimental)){

      if(missing(param_name)){

        if(missing(param_function)){

          stop("You either need to specify: 1) condition, reference, and experimental,
             2) param_name or 3) param_function!")

        }else{

          strategy = "function"

        }

      }else{

        strategy = "single_param"

      }

    }else{

      strategy = "contrast"

    }

  }



  ### Extract kinetic parameters of interest

  # metadf for covariates
  metadf <- obj$metadf


  # Function is in Helpers.R
  if(!is.null(formula_mean)){
    environment(formula_mean) <- NULL
  }
  averages_name <- EZget(obj,
                         type = type,
                         features = features,
                         parameter = parameter,
                         formula_mean = formula_mean,
                         fit_params = fit_params,
                         sd_grouping_factors = sd_grouping_factors,
                         exactMatch = exactMatch,
                         repeatID = repeatID,
                         returnNameOnly = TRUE)


  # Get fractions
  parameter_est <- obj[[type]][[averages_name]]


  # Get features
  if(type == "dynamics"){

    features_to_analyze <- obj[["metadata"]][[type]][[averages_name]][["grouping_features"]]

  }else{

    features_to_analyze <- obj[["metadata"]][[type]][[averages_name]][["features"]]

  }



  ### Perform comparative analysis of interest

  if(strategy == "contrast"){

    ref_mean <- paste0("mean_", design_factor, reference)
    ref_sd <- paste0("sd_", design_factor, reference, "_posterior")
    ref_cov <- paste0("coverage_", design_factor, reference)

    exp_mean <- paste0("mean_", design_factor, experimental)
    exp_sd <- paste0("sd_", design_factor, experimental, "_posterior")
    exp_cov <- paste0("coverage_", design_factor, experimental)

    comparison <- parameter_est %>%
      dplyr::mutate(difference = !!dplyr::sym(exp_mean) - !!dplyr::sym(ref_mean),
                    uncertainty = sqrt( (!!dplyr::sym(ref_sd))^2 + (!!dplyr::sym(exp_sd))^2 ),
                    stat = difference/uncertainty,
                    pval = 2*stats::pnorm(-abs(stat)),
                    avg_coverage = ((!!dplyr::sym(ref_cov)) + (!!dplyr::sym(exp_cov))) / 2) %>%
      dplyr::mutate(padj = stats::p.adjust(pval, method = "BH")) %>%
      dplyr::select(!!features_to_analyze, difference, uncertainty, stat, pval, padj, avg_coverage)


  }

  if(strategy == "single_param"){

    param_mean <- paste0("mean_", param_name)
    param_sd <- paste0("sd_", param_name, "_posterior")
    param_cov <- paste0("coverage_", param_name)

    comparison <- parameter_est %>%
      dplyr::mutate(difference = !!dplyr::sym(param_mean),
                    uncertainty = !!dplyr::sym(param_sd),
                    stat = difference/uncertainty,
                    pval = 2*stats::pnorm(-abs(stat)),
                    avg_coverage = !!dplyr::sym(param_cov)) %>%
      dplyr::mutate(padj = stats::p.adjust(pval, method = "BH")) %>%
      dplyr::select(!!features_to_analyze, difference, uncertainty, stat, pval, padj, avg_coverage)


  }

  if(strategy == "dynamics"){

    # See if design_factor needs to be imputed
    if(missing(design_factor)){

      design_factors <- obj[["metadata"]][[type]][[averages_name]][["dynamics_design_factors"]]

      if(length(design_factors) > 1){
        stop("Need to specify design_factor if you have more than one such
             factor in your 'dynamics' table!")
      }else{

        design_factor <- design_factors

      }


    }

    # Makes more sense for it to be plural in this context
    design_factors <- design_factor



    # Helper function to create a filter expression for the given levels
    create_filter_expr <- function(cols, levels) {
      exprs <- purrr::map2(cols, cols, ~rlang::expr(!!sym(.x) == !!levels[[.y]]))
      purrr::reduce(exprs, `&`)
    }

    # Convert reference_levels and experimental_levels to named lists
    if(is.null(names(reference))){
      reference_list <- stats::setNames(as.list(reference), design_factors)
    }else{
      reference_list <- as.list(reference)
    }

    if(is.null(names(experimental))){
      experimental_list <- stats::setNames(as.list(experimental), design_factors)
    }else{
      experimental_list <- as.list(experimental)
    }


    # Names of new columns
    par_se <- paste0("se_", parameter)


    # Calculate the differences for each unique set of feature values
    comparison <- parameter_est %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(features_to_analyze))) %>%
      dplyr::summarise(
        ref_par = stats::weighted.mean( (!!dplyr::sym(parameter))[!!create_filter_expr(design_factors, reference_list)],
                                w = (!!dplyr::sym(par_se))[!!create_filter_expr(design_factors, reference_list)],
                                na.rm = TRUE),
        exp_par = stats::weighted.mean( (!!dplyr::sym(parameter))[!!create_filter_expr(design_factors, experimental_list)],
                                w = (!!dplyr::sym(par_se))[!!create_filter_expr(design_factors, experimental_list)],
                                na.rm = TRUE),
        ref_se = sqrt(sum( ( (!!dplyr::sym(par_se))[!!create_filter_expr(design_factors, reference_list)])^2 )),
        exp_se = sqrt(sum( ( (!!dplyr::sym(par_se))[!!create_filter_expr(design_factors, experimental_list)])^2 )),
      ) %>%
      dplyr::mutate(
        difference = exp_par - ref_par,
        uncertainty = sqrt(exp_se^2 + ref_se^2),
        stat = difference/uncertainty,
        pval = 2*stats::pnorm(-abs(stat))
      ) %>%
      dplyr::mutate(
        padj = stats::p.adjust(pval, method = "BH")
      )

  }


  if(strategy == "function"){


    stop("Not implemented yet!!!")


  }


  if(normalize_by_median){

    # subtract median difference to account for potential global biases
    comparison <- comparison %>%
      dplyr::ungroup() %>%
      dplyr::mutate(difference = difference - stats::median(difference)) %>%
      dplyr::mutate(
        stat = difference/uncertainty,
        pval = 2*stats::pnorm(-abs(stat))
      ) %>%
      dplyr::mutate(
        padj = stats::p.adjust(pval, method = "BH")
      )

  }



  ### Add output to object

  num_comparisons <- length(obj[['comparisons']])

  output_name <- paste0("comparison", num_comparisons + 1)


  # Impute missing values
  if(missing(design_factor)){
    design_factor <- NULL
  }
  if(missing(reference)){
    reference <- NULL
  }
  if(missing(experimental)){
    experimental <- NULL
  }
  if(missing(param_name)){
    param_name <- NULL
  }
  if(missing(param_function)){
    param_function <- NULL
  }




  if(length(obj[['comparisons']]) > 0){
    output_name <- decide_output(obj, output_name, type = "comparisons",
                                 features = features, parameter = parameter,
                                 design_factor = design_factor,
                                 reference = reference,
                                 experimental = experimental,
                                 param_name = param_name,
                                 param_function = param_function,
                                 cstrat = strategy,
                                 normalize_by_median = normalize_by_median,
                                 overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'comparisons',
                               features = features,
                               parameter = parameter,
                               design_factor = design_factor,
                               reference = reference,
                               experimental = experimental,
                               param_name = param_name,
                               param_function = param_function,
                               cstrat = strategy,
                               normalize_by_median = normalize_by_median,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }

  obj[["comparisons"]][[output_name]] <- dplyr::as_tibble(comparison)

  obj[["metadata"]][["comparisons"]][[output_name]] <- list(features = features,
                                                            parameter = parameter,
                                                            design_factor = design_factor,
                                                            reference = reference,
                                                            experimental = experimental,
                                                            param_name = param_name,
                                                            param_function = param_function,
                                                            features = features_to_analyze,
                                                            parameter = parameter,
                                                            cstrat = strategy,
                                                            normalize_by_median = normalize_by_median,
                                                            repeatID = repeatID)


  if(!methods::is(obj, "EZbakRCompare")){

    class(obj) <- c( "EZbakRCompare", class(obj))

  }

  return(obj)


}



fit_ezbakR_linear_model <- function(data, formula_mean, sd_groups,
                                    coverage_col,
                                    uncertainties_col,
                                    error_if_singular = TRUE){

  parameter <- fsvector <- NULL


  fit <- stats::lm(formula_mean, data,
            singular.ok = !error_if_singular)

  means <- summary(fit)$coef[,"Estimate"]
  names(means) <- paste0("mean_", names(means))


  if(is.null(sd_groups)){

    # My version of robust se estimation: add parameter uncertainty to avoid underestimation
    # I like to be conservative in bioinformatics, because there are always unaccounted
    # for sources of variance.
    X <- stats::model.matrix(fit)
    s2 <- stats::sigma(fit)^2 + data[[uncertainties_col]]^2 + (0.05^2)/(nrow(data)/ncol(X))
    vce <- solve(t(X) %*% X) %*% (t(X) %*% (s2*diag(nrow(data))) %*% X) %*% solve(t(X) %*% X)
    ses <- log(sqrt(diag(vce)))


    # Standard errors in the standard error. Accuracy of this is not super important
    # as it only slighly influences the strength of regularization and nothing else.
    # Result is delta approximation for log(standard devitation)
    replicates <- colSums(X)
    se_logses <- 1/sqrt(2*replicates - 1)
    names(se_logses) <- paste0("se_logse_", names(se_logses))

    # Get average coverage for regression model
    coverage <- rep(log10(mean(data[[coverage_col]])),
                    times = length(ses))

    # Give informative names to output
    names(coverage) <- paste0("coverage_", names(ses))
    names(ses) <- paste0("logse_", names(ses))

    estimates <- c(means,
                   ses, log10(coverage), se_logses)

  }else{

    # Estimate standard deviations in each group
    sds <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(sd_groups))) %>%
      dplyr::mutate(group_sd = stats::sd(!!dplyr::sym(parameter)),
             group_coverage = mean(!!dplyr::sym(coverage_col)))

    # Robust standard errors
    X <- stats::model.matrix(fit)
    s2 <- sds$group_sd^2 + data[[uncertainties_col]]^2 #+ (0.025^2)/(nrow(data)/ncol(X))
    vce <- solve(t(X) %*% X) %*% (t(X) %*% (s2*diag(nrow(data))) %*% X) %*% solve(t(X) %*% X)
    ses <- log(sqrt(diag(vce)))
    names(ses) <- paste0("logse_", names(ses))

    # Standard errors in the standard error. Accuracy of this is not super important
    # as it only slighly influences the strength of regularization and nothing else.
    # Result is delta approximation for log(standard devitation)
    replicates <- colSums(X)
    se_logses <- 1/sqrt(2*replicates - 1)
    names(se_logses) <- paste0("se_logse_", names(se_logses))


    # Idea: Average coverages over group presence in each's parameter column of design
    # matrix.
    coverages <- t(X) %*% matrix(sds$group_coverage, nrow = nrow(data))
    coverages <- coverages / replicates
    rownames(coverages) <- paste0("coverage_", rownames(coverages))

    estimates <- c(means, ses,
                   log10(coverages[,1]),
                   se_logses)

  }

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

