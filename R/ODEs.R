#' Generalized dynamical systems modeling
#'
#' `EZDynamics()` estimates parameters of a user-specified dynamical systems
#' model. The dynamical system model is specified through an adjacency matrix,
#' which is an NxN matrix described below (see `graph` documentation). Modeling
#' can either be done for species all assayed in each sample, or species that
#' are assayed across a set of independent samples (e.g., subcellular fractionation
#' involves assaying different species in different samples).
#'
#' When running `AverageAndRegularize()` to produce input for `EZDynamics()`, you
#' must set `parameter` to "logit_fraction_high<muttype>" (<muttype> = type of mutation
#' modeled by `EstimateFractions()`, e.g., TC), and you must include the label time
#' in your regression formula. `EZDynamics()` models the logit(fraction high <muttype),
#' and this will depend on the label time (longer label time = higher fraction), which
#' is why these two conditions must be met.
#'
#' @param obj Currently must be an EZbakRData object on which `AverageAndRegularize`
#' has been run. In the future, will also support (in case where all species are
#' assayed in every sample) providing output of just `EstimateFractions()` as input,
#' acting as a generalization of `EstimateKinetics()` in that case.
#' @param graph An NxN adjacency matrix, where N represents the number of species
#' being modeled. One of these species must be called "0" and represent the "no RNA"
#' species. This is the species from which some species are synthesized (e.g., 0 -> P,
#' means premature RNA is synthesized from no RNA), and the species to which
#' some species are degraded (e.g., M -> 0 means mature RNA is converted to "no RNA"
#' via degradation). The rows and columns of this matrix must be the names of all
#' modled species, and rownames(graph) == colnames(graph). Entry i,j of the matrix
#' is either 0 if species i cannot be converted into species j under your model,
#' and an integer from 1:npars (where npars = total number of parameters to be
#' estimated) if it can.
#'
#' For example, the model 0 -> P -> M -> 0 would have the `graph`:
#' `matrix(c(0, 1, 0, 0, 0, 2, 3, 0, 0), nrow = 3, ncol = 3, byrow = TRUE`.
#' @param sub_features Which feature columns distinguish between the different
#' measured species? Note, the measured species need not have the same name,
#' and may not be directly equivalent to, the modeled species. The relationship
#' between the modeled species in `graph` and `sub_features` needs to be specified
#' in `modeled_to_measured` if the names are not equivalent though.
#' @param grouping_features Which features are the overarching feature assignments
#' by which `sub_features` should be grouped? This will usually be the feature
#' columns corresponding to full-gene assignments, as well as any higher order
#' assignments (e.g., chromosome). A `sub_feature` can be included in `grouping_features`
#' if it never has the value of `unassigned_name` ("__no_feature" by default). Only
#' one `sub_feature` should ever fulfill this criterion though.
#' @param sample_feature If different samples involve assaying different species,
#' this must be the name of the metadf column that distinguishes the different
#' classes of samples. For example, if analyzing a subcellular fractionation
#' dataset, you likely included a column in your metadf called "compartment".
#' This would then be your `sample_feature`, assuming you ran `AverageAndRegularize()`,
#' with a `mean_formula` that included compartment as a term.
#' @param modeled_to_measured If `sub_features` is not identical to the non "0"
#' species names in `graph`, then you must specify the relationship between
#' `sub_features`, `sample_feature` (if specified), and the species in `graph`.
#' This is done through a list of formulas, whose left side is a `sub_feature`
#' and whose right side is a function of species in `graph`. If `sample_feature`
#' is not specified, then `modeled_to_measured` should be a nameless list of
#' formulae. If `sample_feature` is specified, then `modeled_to_measured` must
#' be a named list of formulas where the names correspond to unique values of
#' `sample_feature`. In this latter case, the elements, should be the mapping
#' of measured features (`sub_features`) to modeled species (names in `graph`
#' that aren't "0").
#'
#' For example, if your model is 0 -> P -> M -> 0, where P stands for premature RNA
#' and M stands for mature RNA, and you have a column called
#' GF that corresponds to assignment anywhere in the gene, and XF that corresponds
#' to assignment of exclusively exon mapping reads, then your `modeled_to_measured`
#' should be `list(GF ~ P + M, XF ~ M).`.
#'
#' As another example, if your model is 0 -> N -> C -> 0, where N stands for
#' "nuclear RNA" and C stands for "cytoplasmic RNA", and your `sample_feature`
#' takes on values of "nuclear", "cytoplasmic", and "total", and you have a single
#' `sub_feature` called XF, then your `modeled_to_measured` shouild be
#' `list(nuclear = GF ~ N, cytoplasmic = GF ~ C, total = GF ~ C + N)`. This is
#' interpreted as meaning in groups of samples for which `sample_feature` == "nuclear",
#' reads assigned to a GF are from the N species (nuclear RNA). When
#' `sample_feature` == "cytoplasmic", reads assigned to GF correspond to the C
#' species (cytoplasmic RNA). When `sample_feature` == "total", reads assigned
#' to GF correspond to a combination of N and C (nuclear and cytoplasmic RNA).
#' @param parameter_names Vector of names you would like to give to the estimated
#' parameters. ith element should correspond to name of parameter given the ID
#' i in `graph`. By default, this is just ki, where i is this numerical index.
#' @param unassigned_name What value will a `sub_feature` column have if a read
#' was not assigned to said feature? "__no_feature" by default.
#' @param type What type of table would you like to use? Currently only supports
#' "averages", but will support "fractions" in the near future.
#' @param features Character vector of the set of features you want to stratify
#' reads by and estimate proportions of each RNA population. The default of "all"
#' will use all feature columns in the `obj`'s cB.
#' @param populations Mutational populations that were analyzed to generate the
#' fractions table to use. For example, this would be "TC" for a standard
#' s4U-based nucleotide recoding experiment.
#' @param fraction_design "Design matrix" specifying which RNA populations exist
#' in your samples. By default, this will be created automatically and will assume
#' that all combinations of the `mutrate_populations` you have requested to analyze are
#' present in your data. If this is not the case for your data, then you will have
#' to create one manually. See docs for `EstimateFractions` (run ?EstimateFractions()) for more details.
#' @param parameter Parameter to average across replicates of a given condition.
#' Has to be "logit_fraction_high<muttype>", where <muttype> is the type of
#' mutation modeled in `EstimateFractions()` (e.g, TC) in this case.
#' @param mean_vars Sample features from metadf that were used in formula for
#' parameter average linear model.
#' @param sd_vars Sample features from metadf that were used in formula for
#' parameter standard deviation linear model.
#' @param repeatID If multiple `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @importFrom magrittr %>%
#' @export
EZDynamics <- function(obj,
                     graph,
                     sub_features,
                     grouping_features,
                     sample_feature = NULL,
                     modeled_to_measured = NULL,
                     parameter_names = paste0("k", 1:max(graph)),
                     unassigned_name = "__no_feature",
                     type = "averages",
                     features = NULL,
                     populations = NULL,
                     fraction_design = NULL,
                     parameter = NULL,
                     mean_vars = NULL,
                     sd_vars = NULL,
                     repeatID = NULL,
                     exactMatch = TRUE,
                     overwrite = TRUE){

  ##### ORDER OF OPERATIONS
  # 1) Infer homogeneous ODE system matrix representation (A)
  # 2) Fit model, lol.

  ##### Challenges to overcome
  # 1) Need to generalize steady-state estimation which involves figuring
  # out what the vector of zeroth-order parameters looks like
  # 2) Going to use averages as input, so need to figure out how to best
  # model average of normalized read counts.
  # 3) Still struggling with a complete, efficient generalization of inference
  # of the matrix A. NOT ANY MORE: Create a diagonal matrix from the row
  # sums of the reduced adjacency matrix and add this to the transpose of the reduced
  # adjacency matrix.

  ##### DESIRED FEATURES
  # 1) Should be able to take averages or fractions as input. Latter only applicable
  # if no cross-sample integration needs to be done (e.g., modeling P -> M).
  # 2) Would love to assess dynamic range of parameter estimates using log-likelihood
  # profile to let users know what they can say about certain parameters

  ### Input

  # Can try to infer replicate numbers to get more informative coverage likelihood
  metadf <- obj$metadf
  meta_groups <- c("tl", sample_feature)
  reps <- metadf %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(meta_groups))) %>%
    dplyr::summarise(nreps = dplyr::n())


  # Tidy averages table (currently makes hard assumption about interaction terms being only terms)
  table_name <- EZget(obj = obj, type = type, features = features,
                      populations = populations, fraction_design = fraction_design,
                      parameter = parameter, mean_vars = mean_vars,
                      sd_vars = sd_vars, repeatID = repeatID,
                      exactMatch = exactMatch, returnNameOnly = TRUE)

  if(length(table_name) > 1){
    stop("More than one table matches your search criterion!")
  }

  table <- obj[[type]][[table_name]]

  if(type == "averages"){

    pivot_columns <- obj[['metadata']][['averages']][[table_name]][['mean_vars']]
    features <- obj[['metadata']][['averages']][[table_name]][['features']]

    # Currently making hard assumption of interactions only...
    # TO-DO: Loosen this restriction
    # Also assuming that mean_formula == sd_formula
    pattern <- paste0("(.*)_", paste(paste0(pivot_columns, "(.+)"), collapse = ":"))

    # Need to remove _posterior suffix (maybe shouldn't be there at all?)
    # to make general pivoting strategy work
    colnames(table) <- gsub("_posterior$", "", colnames(table))

    # Tidy data!
    tidy_avgs <- table %>%
      tidyr::pivot_longer(
        cols = -!!features,
        names_to = c(".value", pivot_columns),
        names_pattern = pattern
      )


    # formula list
    if(is.null(modeled_to_measured)){

      species <- rownames(graph)[rownames(graph) != "0"]
      modeled_to_measured <- vector(mode = "list", length = length(species))
      for(s in seq_along(species)){
        modeled_to_measured[[s]] <- as.formula(paste0(species[s], "~", species[s]))
      }

    }

    if(is.null(names(modeled_to_measured))){

      modeled_to_measured <- list(total = modeled_to_measured)
      tidy_avgs <- tidy_avgs %>%
        dplyr::mutate(imputed_sample_feature = 'total')

      sample_feature <- "imputed_sample_feature"

    }


    ##### STEPS
    # 1) Infer general solution as I do in simulation
    # 2) Probably need to tidy up the averages table, which may be easier said then done
    # There should be a mean, sd, and coverage column, as well as the feature columns,
    # and finally any group detail columns (tl, compartment, etc.)



    ### FIT MODEL

    # Determine initial vector of parameters and their upper and lower bounds
    npars <- max(graph)

    # Have to add a bit of noise because equal parameters means general solution
    # breaks down as their aren't N eigenvalues. Could generalize to case
    # of equal parameters, but thinking of parameter estimates as continuous
    # random variables that thus have probability of 0 of being equal.
    lower_bounds <- rnorm(npars,
                          -10,
                          0.01)
    upper_bounds <- rnorm(npars,
                          10,
                          0.01)
    starting_values <- rnorm(npars,
                             0, 0.01)


    ### Infer which measured species each row of tidy_avgs belongs to
    feature_mat <- tidy_avgs %>%
      dplyr::select(!!sub_features) %>%
      as.matrix()

    feature_mat <- feature_mat == unassigned_name

    # Are any species always assigned?
    fm_cols <- colnames(feature_mat)
    always_assigned <- fm_cols[which(colSums(feature_mat) == nrow(feature_mat))]

    if(length(always_assigned) > 1){

      stop("There are multiple sub_features that never have a value of `unassigned_name`!
           There is not way for EZbakR to determine which measured species a row
           corresponds to in this case.")

    }else if(length(always_assigned) == 0){

      sometimes_assigned <- fm_cols

    }else{

      sometimes_assigned <- fm_cols[fm_cols != always_assigned]

    }


    measured_species <- apply(feature_mat, 1, function(row) {


      valid_cols <- sometimes_assigned[!row[!(fm_cols %in% always_assigned)]]

      if(length(valid_cols) > 0){

        valid_cols[1]

      }else{

        always_assigned

      }

    })

    tidy_avgs$measured_specie <- measured_species


    # Fit model
    cols_to_group_by <- c(grouping_features,
                       pivot_columns[!(pivot_columns %in% c("tl", sample_feature))])


    dynfit <- tidy_avgs  %>%
      dplyr::mutate(tl = as.numeric(tl)) %>%
      dplyr::inner_join(reps,
                        by = meta_groups) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group_by))) %>%
      dplyr::summarise(fit = list(I(stats::optim(starting_values,
                                          fn = dynamics_likelihood,
                                          graph = graph,
                                          formula_list = modeled_to_measured,
                                          logit_fn = mean,
                                          logit_fn_sd = sd,
                                          coverage = coverage,
                                          nreps = nreps,
                                          tls = as.numeric(tl),
                                          sample_features = !!dplyr::sym(sample_feature),
                                          feature_types = measured_specie,
                                          lower = lower_bounds,
                                          upper = upper_bounds,
                                          method = "L-BFGS-B",
                                          hessian = TRUE,
                                          use_coverage = TRUE))))


    # Get parameter estimates and uncertainties
    for(n in 1:npars){

      par_name <- parameter_names[n]
      par_se_name <- paste0(par_name, "_se")

      dynfit <- dynfit %>% dplyr::mutate(
        !!par_name := purrr::map_dbl(fit, ~ .x$par[n]),
        !!par_se_name := tryCatch(
          {
          purrr::map_dbl(fit, ~ sqrt(diag(solve(.x$hessian)))[n])
          },
          error = function(e) {
            Inf
          }
          ),
      )
    }


  }else{


    features_to_analyze <- obj[["metadata"]][["fractions"]][[table_name]][["features"]]


    # Determine which column to use for kinetic parameter estimation
    fraction_cols <- colnames(table)

    fraction_of_interest <- fraction_cols[grepl("^logit_fraction_high", fraction_cols)]
    fractionse <- fraction_cols[grepl("^se_logit_fraction_high", fraction_cols)]


    ### Normalize read counts


    # TO-DO: ALLOW USERS TO JUST USE THE TPM FROM ISOFORM QUANTIFICATION
    reads_norm <- get_normalized_read_counts(obj = obj,
                                             features_to_analyze = features_to_analyze,
                                             fractions_name = table_name)


    ### Prep tables for kinetic parameter estimation

    kinetics <- setDT(data.table::copy(table))

    # Add label time info
    metadf <- obj$metadf

    metadf <- setDT(data.table::copy(metadf))
    setkey(metadf, sample)
    setkey(kinetics, sample)

    kinetics <- kinetics[metadf[,c("sample", "tl")], nomatch = NULL]

    # Make sure no -s4U controls made it through
    kinetics <- kinetics[tl > 0]


    ### Infer which measured species each row belongs to
    feature_mat <- kinetics %>%
      dplyr::select(!!sub_features) %>%
      as.matrix()

    feature_mat <- feature_mat == unassigned_name

    # Are any species always assigned?
    fm_cols <- colnames(feature_mat)
    always_assigned <- fm_cols[which(colSums(feature_mat) == nrow(feature_mat))]

    if(length(always_assigned) > 1){

      stop("There are multiple sub_features that never have a value of `unassigned_name`!
           There is not way for EZbakR to determine which measured species a row
           corresponds to in this case.")

    }else if(length(always_assigned) == 0){

      sometimes_assigned <- fm_cols

    }else{

      sometimes_assigned <- fm_cols[fm_cols != always_assigned]

    }


    measured_species <- apply(feature_mat, 1, function(row) {


      valid_cols <- sometimes_assigned[!row[!(fm_cols %in% always_assigned)]]

      if(length(valid_cols) > 0){

        valid_cols[1]

      }else{

        always_assigned

      }

    })

    kinetics$measured_specie <- measured_species


    ### Impute dummy value for modeled_to_measured
    if(is.null(names(modeled_to_measured))){

      modeled_to_measured <- list(total = modeled_to_measured)
      kinetics <- kinetics %>%
        dplyr::mutate(imputed_sample_feature = 'total')

      sample_feature <- "imputed_sample_feature"

    }

    ### Fit model

    # Determine initial vector of parameters and their upper and lower bounds
    npars <- max(graph)

    # Have to add a bit of noise because equal parameters means general solution
    # breaks down as their aren't N eigenvalues. Could generalize to case
    # of equal parameters, but thinking of parameter estimates as continuous
    # random variables that thus have probability of 0 of being equal.
    lower_bounds <- rnorm(npars,
                          -10,
                          0.01)
    upper_bounds <- rnorm(npars,
                          10,
                          0.01)
    starting_values <- rnorm(npars,
                             0, 0.01)

    cols_to_group_by = c("sample", grouping_features)

    dynfit <- kinetics %>%
      dplyr::mutate(tl = as.numeric(tl)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group_by))) %>%
      dplyr::summarise(fit = list(I(stats::optim(starting_values,
                                                 fn = dynamics_likelihood,
                                                 graph = graph,
                                                 formula_list = modeled_to_measured,
                                                 logit_fn = !!dplyr::sym(fraction_of_interest),
                                                 logit_fn_sd = !!dplyr::sym(fraction_se),
                                                 coverage = n,
                                                 nreps = 1,
                                                 tls = as.numeric(tl),
                                                 sample_features = !!dplyr::sym(sample_feature),
                                                 feature_types = measured_specie,
                                                 lower = lower_bounds,
                                                 upper = upper_bounds,
                                                 method = "L-BFGS-B",
                                                 hessian = TRUE,
                                                 use_coverage = TRUE))),
                       n = mean(n) # I'm not sure how to calculate read counts to be passed to AverageAndRegularize(). I think average should be fine though
                       )


    # Get parameter estimates and uncertainties
    for(n in 1:npars){

      par_name <- parameter_names[n]
      par_se_name <- paste0(par_name, "_se")

      dynfit <- dynfit %>% dplyr::mutate(
        !!par_name := purrr::map_dbl(fit, ~ .x$par[n]),
        !!par_se_name := tryCatch(
          {
            purrr::map_dbl(fit, ~ sqrt(diag(solve(.x$hessian)))[n])
          },
          error = function(e) {
            Inf
          }
        ),
      )

    }


  }


  ### Add output to object

  num_dynamics <- length(obj[['dynamics']])

  output_name <- paste0("dynamics", num_dynamics + 1)


  if(num_dynamics > 0){
    output_name <- decide_output(obj, output_name, type = "dynamics",
                                 features = features,
                                 sub_features = sub_features,
                                 grouping_features = grouping_features,
                                 sample_feature = sample_feature,
                                 graph = graph, parameter = parameter,
                                 mean_vars = mean_vars, sd_vars = sd_vars)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'dynamics',
                               features = features,
                               sub_features = sub_features,
                               grouping_features = grouping_features,
                               sample_feature = sample_feature,
                               graph = graph, parameter = parameter,
                               mean_vars = mean_vars, sd_vars = sd_vars,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }

  obj[["dynamics"]][[output_name]] <- dplyr::as_tibble(dynfit)

  obj[["metadata"]][["dynamics"]][[output_name]] <- list(features = features,
                                                         sub_features = sub_features,
                                                         grouping_features = grouping_features,
                                                         sample_feature = sample_feature,
                                                         graph = graph, parameter = parameter,
                                                         mean_vars = mean_vars, sd_vars = sd_vars,
                                                         repeatID = repeatID)


  if(!methods::is(obj, "EZbakRDynamics")){

    class(obj) <- c( "EZbakRDynamics", class(obj))

  }


  return(obj)

}



# Function to relate measured species to actual species by evaluating the
# user provided formulas
evaluate_formulas2 <- function(original_vector, formula) {

  new_vector <- numeric()
  response <- as.character(formula[[2]])
  terms <- all.vars(formula[[3]])
  expr <- formula[[3]]
  env <- list2env(as.list(original_vector))
  value <- eval(expr, envir = env)

  # Will return vector with names of measured species
  new_vector[response] <- value
  return(new_vector)

}

# Likelihood function for generalized dynamical systems modeling
# with averages tables
dynamics_likelihood <- function(parameter_ests, graph, formula_list = NULL,
                                logit_fn, logit_fn_sd, coverage, nreps = 2,
                                tls, sample_features, feature_types,
                                use_coverage = TRUE, alt_coverage = FALSE,
                                coverage_sd = NULL, scale_factor = NULL){

  ### Step 0, check to see if single replicate of data is being passed
  if(nreps == 1){
    single_replicate <- TRUE
  }


  ### Step 1, construct A

  # Parameters are on log-scale for ease of optimization
  param_extend <- c(0, exp(parameter_ests))
  param_graph <- matrix(param_extend[graph + 1],
                        nrow = nrow(graph),
                        ncol = ncol(graph),
                        byrow = FALSE)

  # A is same size as graph minus the "0" row and column
  A <- matrix(0,
              nrow = nrow(graph) - 1,
              ncol = ncol(graph) - 1)

  rownames(A) <- rownames(graph[-1,])
  colnames(A) <- rownames(graph[-1,])



  zero_index <- which(colnames(graph) == "0")

  # Left as an exercise to the reader
  # Just kidding, I'll document this somewhere.
  # TO-DO: Should really compile documentation for all
  # non-obvious mathematical results
  diag(A) <- -rowSums(param_graph[-zero_index,])
  A <- A + t(param_graph[-zero_index,-zero_index])


  ### Step 2: infer general solution

  Rss <- solve(a = A,
               b = -param_graph[zero_index,-zero_index])


  ev <- eigen(A)

  lambda <- ev$values
  V<- ev$vectors
  cs <- solve(V, -Rss)


  ### Step 3: Infer data for actual measured species
  all_ss <- c()
  all_fns <- c()

  for(n in seq_along(tls)){

    tl <- tls[n]
    sample_feature <- sample_features[n]
    feature_type <- feature_types[n]


    exp_lambda <- exp(lambda*tl)

    scaled_eigenvectors <- V %*% diag(exp_lambda*cs)

    result_vector <- rowSums(scaled_eigenvectors) + Rss

    names(result_vector) <- rownames(A)

    # Evaluate the formulas
    sample_formula <- formula_list[[sample_feature]]

    sample_formula <- sample_formula[[which(sapply(1:length(sample_formula),
                                                   function(x) all.vars(sample_formula[[x]])[1] == feature_type))]]

    measured_levels <- evaluate_formulas2(result_vector, sample_formula)

    names(Rss) <- rownames(A)
    measured_ss <- evaluate_formulas2(Rss, sample_formula)

    all_fns <- c(all_fns, measured_levels/measured_ss)
    all_ss <- c(all_ss, measured_ss)

  }

  all_fns <- dplyr::case_when(
    all_fns > 0.9999 ~ 0.9999,
    all_fns < 0.0001 ~ 0.0001,
    .default = all_fns
  )


  ### Step 4 calculate likelihood
  ll <- stats::dnorm(logit_fn,
              logit(all_fns),
              logit_fn_sd,
              log = TRUE)

  if(use_coverage){

    if(alt_coverage){

      all_reads <- all_ss*scale_factor

      ll <- ll +
        stats::dnorm(coverage,
                     log10(all_ss),
                     sd = coverage_sd,
                     log = TRUE)
      # ll <- ll +
      #   stats::dnorm(coverage,
      #                log10(all_reads),
      #                (1/((10^coverage)*(log(10)^2)))/sqrt(nreps),
      #                log = TRUE)


    }else if(single_replicate){

      ll <- ll +
        stats::dpois(coverage,
                     log = TRUE)


    } else{

      ll <- ll +
        stats::dnorm(coverage,
                     log10(all_ss),
                     sqrt((1/((10^coverage)*(log(10)^2))))/sqrt(nreps),
                     log = TRUE)

    }

  }

  if(!is.finite(sum(ll))){
    browser()
  }

  return(-sum(ll))



}
