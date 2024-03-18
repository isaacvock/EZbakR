#' Average parameter estimates across replicates, and regularize variance estimates
#'
#' `AverageAndRegularize` uses a linear model to average estimates of a parameter
#' of interest over replicates, and to get averages for all of a set of conditions
#' specified by the user. This specification is done through formulas that will
#' be used to create a design matrix for parameters to estimate.
#'
#' @param obj An `EZbakRFractions` object, which is an `EZbakRData` object on
#' which `EstimateFractions()` has been run.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param parameter Parameter to average across replicates of a given condition.
#' @param formula_mean An R formula object specifying how the `parameter` of interest
#' depends on the sample characteristics specified in `obj`'s metadf. The most formula
#' will be `~ treatment` or `~ treatment + batch`, where `treatment` and `batch` would
#' be replaced with whatever you called the main treatment and batch identifiers in
#' your metadf.
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
#' @import data.table
#' @importFrom magrittr %>%
AverageAndRegularize <- function(obj, features = "all", parameter = "log_kdeg",
                            formula_mean = NULL, formula_sd = NULL,
                            include_all_parameters = TRUE,
                            sd_reg_factor = 10,
                            error_if_singular = TRUE){



  if(parameter == "log_TILAC_ratio"){

    obj <- general_avg_and_reg(obj = obj, features = features, parameter = parameter,
                             formula_mean = formula_mean, formula_sd = formula_sd,
                             include_all_parameters = include_all_parameters,
                             sd_reg_factor = sd_reg_factor,
                             error_if_singular = error_if_singular,
                             TILAC = TRUE)

  }else{

    obj <- general_avg_and_reg(obj = obj, features = features, parameter = parameter,
                        formula_mean = formula_mean, formula_sd = formula_sd,
                        include_all_parameters = include_all_parameters,
                        sd_reg_factor = sd_reg_factor,
                        error_if_singular = error_if_singular)

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
# This is on log-scale, so have to return on natural
# scale to make it usable in downstream test statistic
# calculations.
get_sd_posterior <- function(n = 1, sd_est, sd_var,
                             fit_var, fit_mean){

  denom <- (n/sd_var + 1/fit_var)
  num <- sd_est/sd_var + fit_mean/fit_var

  return(num/denom)

}



general_avg_and_reg <- function(obj, features, parameter,
                                formula_mean, formula_sd,
                                include_all_parameters,
                                sd_reg_factor,
                                error_if_singular,
                                TILAC = FALSE){


  ### Extract kinetic parameters of interest

  # metadf for covariates
  metadf <- obj$metadf

  # cB columns
  cB <- obj$cB
  cB_cols <- colnames(cB)

  # Mutation columns in cB
  mutcounts_in_cB <- find_mutcounts(obj)

  # Base count columns in cB
  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # feature columns
  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB,
                                             basecounts_in_cB,
                                             "sample", "n"))]

  # Need to determine which columns of the cB to group reads by
  if(features == "all"){

    features_to_analyze <- features_in_cB

  }else{

    if(!all(features %in% features_in_cB)){

      stop("features includes columns that do not exist in your cB!")

    }else{

      features_to_analyze <- features

    }

  }

  # Get the kinetic parameter data frame
  kinetics_name <- paste(c("kinetics", features_to_analyze), collapse = "_")
  kinetics <- obj[[kinetics_name]]
  kinetics <- kinetics %>%
    dplyr::mutate(log_normalized_reads = log10(normalized_reads))



  # Add kinetic parameter column to formula
  formula_mean <- as.formula(paste(c(parameter, formula_mean), collapse = ""))
  formula_reads <- as.formula(paste(c("log_normalized_reads", formula_sd), collapse = ""))
  formula_sd <- as.formula(paste(c(parameter, formula_sd), collapse = ""))

  ### Fit linear, potentially heteroskedastic model
  meta_cols <- colnames(metadf)
  tl_cols <- meta_cols[grepl("^tl", meta_cols)]

  # Add covariates to kinetics
  kinetics <- kinetics %>%
    dplyr::inner_join(metadf %>%
                 dplyr::select(-!!tl_cols), by = "sample")


  ### Need to coverage for edge case where there are single level factors
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

    model_fit <- kinetics %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(features_to_analyze))) %>%
      dplyr::summarise()

  }else{

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

  }



  ### Estimate coverage vs. variance trends for each standard deviation estimate

  message("")
  message("Estimating coverage vs. variance trend")

  # Step 1: Filter column names for relevant patterns
  sd_columns <- names(model_fit)[grepl("^logsd_", names(model_fit))]
  covariate_names <- substring(sd_columns, 7)
  coverage_columns <- names(df)[grepl("^coverage_", names(model_fit))]
  relevant_columns <- union(sd_columns, coverage_columns)

  # Step 2: Iterate and perform regression
  regression_results <- purrr::map(covariate_names, ~ {

    # Dynamically create formula
    formula_str <- paste("logsd_", .x, " ~ coverage_", .x, sep = "")
    formula <- as.formula(formula_str)

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
      dplyr::mutate(!!col_name := get_sd_posterior(sd_est = !!dplyr::sym(sd_est_name),
                                                   sd_var = (!!dplyr::sym(sd_var_name)) ^ 2,
                                                   fit_var = var(regression_results[[c]]$lm_result$residuals) / sd_reg_factor,
                                                   fit_mean = !!dplyr::sym(fit_mean_name))) %>%
      dplyr::mutate(!!natural_col_name := exp((!!dplyr::sym(col_name))))

    mean_est <- paste0("mean_", c)
    coverage_est <- paste0("coverage_", c)

    cols_to_keep <- c(cols_to_keep, mean_est, natural_col_name)
    coverage_cols <- c(coverage_cols, coverage_est)
  }


  # Want to stick the coverage columns at the end
  final_output <- model_fit %>%
    dplyr::select(!!features_to_analyze, !!cols_to_keep, !!coverage_cols)



  # Prep output
  output_name <- paste0("average_", parameter, "_", features_to_analyze)

  obj[[output_name]] <- final_output

  if(include_all_parameters){

    output_name <- paste0("fullfit_average_", parameter, "_", features_to_analyze)
    obj[[output_name]] <- model_fit

  }

  if(!is(obj, "EZbakRFit")){

    class(obj) <- c( "EZbakRFit", class(obj))

  }

  return(obj)


}


#' Get contrasts of estimated parameters
#'
#' @param obj An `EZbakRFit` object, which is an `EZbakRFractions` object on
#' which `AverageAndRegularize()` has been run.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param parameter Parameter to average across replicates of a given condition.
#' @import data.table
#' @importFrom magrittr %>%
CompareParameters <- function(obj, features = "all", parameter = "log_kdeg",
                              condition, reference, experimental){

  ### Extract kinetic parameters of interest

  # metadf for covariates
  metadf <- obj$metadf

  # cB columns
  cB <- obj$cB
  cB_cols <- colnames(cB)

  # Mutation columns in cB
  mutcounts_in_cB <- find_mutcounts(obj)

  # Base count columns in cB
  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  # feature columns
  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB,
                                             basecounts_in_cB,
                                             "sample", "n"))]

  # Need to determine which columns of the cB to group reads by
  if(features == "all"){

    features_to_analyze <- features_in_cB

  }else{

    if(!all(features %in% features_in_cB)){

      stop("features includes columns that do not exist in your cB!")

    }else{

      features_to_analyze <- features

    }

  }

  # Get the kinetic parameter data frame
  parameter_name <- paste(c("average", parameter, features_to_analyze), collapse = "_")
  parameter_est <- obj[[parameter_name]]


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
           pval = 2*pnorm(-abs(stat)),
           avg_coverage = ((!!dplyr::sym(ref_cov)) + (!!dplyr::sym(exp_cov))) / 2) %>%
    dplyr::mutate(padj = stats::p.adjust(pval, method = "BH")) %>%
    dplyr::select(!!features_to_analyze, difference, uncertainty, stat, pval, padj, avg_coverage)



  ### Add output to object

  output_name <- paste0(paste(features_to_analyze, collapse = "_"), "_",
                        parameter, "_", condition, "_", experimental, "_vs_", reference)

  obj[["comparisons"]][[output_name]] <- comparison


  if(!is(obj, "EZbakRCompare")){

    class(obj) <- c( "EZbakRCompare", class(obj))

  }

  return(obj)


}


# Define the likelihood function
heteroskedastic_likelihood <- function(params, y, X_mean, X_sd) {
  n <- length(y)
  beta <- params[1:ncol(X_mean)]
  log_sigma <- params[(ncol(X_mean) + 1):length(params)]

  mu <- X_mean %*% beta
  sigma <- exp(X_sd %*% log_sigma)

  # Gaussian log-likelihood
  logL <- -sum(dnorm(y, mean=mu, sd=sigma, log=TRUE))
  return(logL)
}


fit_heteroskedastic_linear_model <- function(formula_mean, formula_sd, data,
                                             dependent_var,
                                             error_if_singular = TRUE) {


  # Parse formula objects into model matrices
  designMatrix_mean <- model.matrix(formula_mean, data)
  designMatrix_sd <- model.matrix(formula_sd, data)


  if(error_if_singular){

    # Check for singularity in the design matrices
    if(qr(designMatrix_mean)$rank < ncol(designMatrix_mean)) {
      stop("Design matrix for mean is singular")
    }
    if(qr(designMatrix_sd)$rank < ncol(designMatrix_sd)) {
      stop("Design matrix for standard deviation is singular")
    }

  }


  # Initial parameter guesses
  startParams <- rep(0, ncol(designMatrix_mean) + ncol(designMatrix_sd))

  # Optimization
  opt <- optim(startParams, heteroskedastic_likelihood,
               y = data[[dependent_var]],
               X_mean = designMatrix_mean,
               X_sd = designMatrix_sd,
               method = "BFGS",
               hessian = TRUE)



  if(opt$convergence != 0) {
    stop("Model did not converge!")
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
    # Call your function inside tryCatch
    solve(opt$hessian)
  }, error = function(e) {
    # If an error occurs, print the error message
    print(e)
    # Then enter the RStudio debugger
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
  fit <- lm(formula = formula,
            data = data)


  estimates <- fit$coefficient
  names(estimates) <- c(paste0("coverage_", names(estimates)))

  return(as.list(estimates))

}
