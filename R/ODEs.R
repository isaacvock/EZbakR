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
#' modeled by `EstimateFractions()`, e.g., TC). If you have multiple distinct label times,
#' you must also include the label time (`tl` of your `metadf`)
#' in your regression formula. `EZDynamics()` models the logit(fraction high <muttype),
#' and this will depend on the label time (longer label time = higher fraction), which
#' is why these two conditions must be met. If you only have a single label time though,
#' `EZDynamics` will be able to impute this one value for all samples from your `metadf`.
#' You can also include additional interaction terms in your `AverageAndRegularize()`
#' model for different experimental conditions in which experiments were conducted,
#' so that inferred kinetic parameters can be compared across these conditions. Currently,
#' more complex modeling beyond simple stratification of samples by one or more condition
#' is not possible with `EZDynamics()`.
#'
#' For normalization purposes, especially if analyzing pre-mRNA processing dynamics,
#' you will need to provide `AverageAndRegularize()` with a table of feature lengths
#' via the `feature_lengths` parameter. This will be used in all cases to length
#' normalize read counts. Even in the case when you are just modeling mature mRNA
#' dynamics, this is technically necessary for accurate estimation of scale factors.
#'
#' The first step of `EZDynamics()` is attempted inference of normalization scale
#' factors for read counts. If you have scale factors you calculated yourself,
#' e.g. via specialized spike-in protocols, you can provide these via the `scale_factors`
#' parameter. If not, `EZDynamics()` will try to infer these from the fraction labeled's
#' in each `sample_feature` (e.g., in different subcellular compartments). This relies on
#' having some `sample_feature`'s that are a combination of other `sample_feature`'s. For
#' example, if analyzing subcellular fractionation data, you may have 1) total RNA; 2)
#' cytoplasmic RNA; and 3) nuclear RNA. Total RNA = cytoplasmic + nuclear RNA, and thus
#' the fraction of reads from labeled RNA is a function of the total cytoplasmic and
#' nuclear fraction labeled's, as well as the relative molecular abundances of cytoplasmic
#' and nuclear RNA. The latter is precisely the scale factors we need to estimate.
#' If you do not have sufficient combinations of data to perform this scale factor estimation,
#' `EZDynamics()` will only use the fraction labeled's for modeling kinetic parameters.
#' It can then perform post-hoc normalization to estimate a single synthesis rate constant,
#' using the downstream rate constants to infer the unknown normalization scale factor necessary
#' to combine kinetic parameter estimates and read counts to infer this rate constant.
#'
#' For estimating kinetic parameters, `EZDynamics()` infers the solution of the linear
#' system of ODEs specified in your `graph` matrix input. This is done by representing
#' the system of equations as a matrix, and deriving the general solution of the levels
#' of each modeled species of RNA from the eigenvalues and eigenvectors of this matrix.
#' While this makes `EZDynamics()` orders of magnitude more efficient than if it had to
#' numerically infer the solution for each round of optimization, needing to compute
#' eigenvalues and eigenvectors in each optimization iteration is still non-trivial,
#' meaning that `EZDynamics()` may take anywhere from 10s of minutes to a couple hours
#' to run, depending on how complex your model is and how many distinct set of samples
#' and experimental conditions you have.
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
#' @param scale_factors Data frame mapping samples to factors by which to multiply
#' read counts so as ensure proper normalization between different RNA populations.
#' Only relevant if you are modeling relationships between distinct RNA populations,
#' for example RNA from nuclear and cytoplasmic fractions. Will eventually be inferred
#' automatically.
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
#' @param prior_means Mean of log-Normal prior for kinetic parameters. Should
#' be vector where ith value is mean for ith parameter, i = index in `graph`
#' @param prior_sds Std. dev. of log-Normal prior for kinetic parameter.
#' Should be vector where ith value is mean for ith parameter, i = index in `graph`.
#' @param avg_params_tokeep Names of parameters in averages table that you would
#' like to keep. Other parameters will be discarded. Don't include the prefixes "mean_",
#' "sd_", or "coverage_"; these will be imputed automatically. In other words,
#' this should be the base parameter names.
#' @param avg_params_todrop Names of parameters in averages table that you would
#' like to drop. Other parameters will be kept. Don't include the prefixes "mean_",
#' "sd_", or "coverage_"; these will be imputed automatically. In other words,
#' this should be the base parameter names.
#' @param label_time_name Name of relevant label time column that will be found
#' in the parameter names. Defaults to the standard "tl".
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
#' @param repeatID If multiple `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param feature_lengths Table of effective lengths for each feature combination in your
#' data. For example, if your analysis includes features named GF and XF, this
#' should be a data frame with columns GF, XF, and length.
#' @param use_coverage If TRUE, normalized read counts will be used to inform
#' kinetic parameter estimates. If FALSE, only fraction news will be used, which
#' will leave some parameters (e.g., synthesis rate) unidentifiable, though has
#' the advantage of avoiding the potential biases induced by normalization problems.
#' @param normalization_repeatID For extracting the `fractions` table needed for
#' normalization of multi-sample data. If multiple `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param normalization_exactMatch For extracting the `fractions` table needed for
#' normalization of multi-sample data. If TRUE, then `features` and `populations`
#' have to exactly match those for a given fractions table for that table to be used.
#' Means that you can't specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param species_to_sf List mapping individual RNA species in `graph` to
#' different sample_feature values (sf). This is relevant if you are modeling both
#' pre- and mature RNA dynamics in subcellular fractionation data. EZbakR can usually
#' automatically infer this, but if not, then you can manually specify this mapping.
#' Should be a list with one element per unique sample_feature type. Each element should
#' be a vector of the RNA species in `graph` that belong to that sample_feature type. For
#' example, if you have whole cell, cytoplasmic, and nuclear fraction data, this should have
#' one element called "cytoplasmic" and one element called "nuclear". The whole cell
#' sample_feature is a sum of cytoplasmic and nuclear and thus does not apply. The "cytoplasmic"
#' element of this list should be the set of RNA species that are present in the cytoplasm,
#' e.g. cytoplasmic pre-RNA (maybe referred to in `graph` as CP) and cytoplasmic mature
#' RNA (maybe referred to in `graph` as CM).
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @importFrom magrittr %>%
#' @import data.table
#' @return `EZbakRData` object with an additional "dynamics" table.
#' @export
EZDynamics <- function(obj,
                     graph,
                     sub_features,
                     grouping_features,
                     scale_factors = NULL,
                     sample_feature = NULL,
                     modeled_to_measured = NULL,
                     parameter_names = paste0("logk", 1:max(graph)),
                     unassigned_name = "__no_feature",
                     type = "averages",
                     prior_means = rep(-3, times = max(graph)),
                     prior_sds = c(3, rep(1, times = max(graph)-1)),
                     avg_params_tokeep = NULL,
                     avg_params_todrop = NULL,
                     label_time_name = "tl",
                     features = NULL,
                     populations = NULL,
                     fraction_design = NULL,
                     parameter = NULL,
                     repeatID = NULL,
                     exactMatch = TRUE,
                     feature_lengths = NULL,
                     use_coverage = TRUE,
                     normalization_repeatID = NULL,
                     normalization_exactMatch = TRUE,
                     species_to_sf = NULL,
                     overwrite = TRUE){


  # Hack to deal with devtools::check() NOTEs
  param_type <- coverage <- sd <- nreps <- tl <- measured_specie <- NULL
  fit <- ode_models <- GF <- NULL

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
  meta_groups <- c(label_time_name, sample_feature)



  # Tidy averages table (currently makes hard assumption about interaction terms being only terms)
  table_name <- EZget(obj = obj, type = type, features = features,
                      populations = populations, fraction_design = fraction_design,
                      parameter = parameter, repeatID = repeatID,
                      exactMatch = exactMatch, returnNameOnly = TRUE)

  if(length(table_name) > 1){
    stop("More than one table matches your search criterion!")
  }else if(is.null(table_name)){
    if(type == "fractions"){

      stop("No table matches your search criterion! Have you run EstimateFraction() yet?")

    }else if (type == "averages"){

      stop("No table matches your search criterion! Have you run AverageAndRegularize() yet?")

    }else{

      stop("type must be 'averages' or 'fractions'!")

    }
  }

  table <- obj[[type]][[table_name]]

  if(type == "averages"){


    ### Process input data for model

    # Is there only one label time? will have to impute later
    label_times <- metadf %>%
      dplyr::filter(!!dplyr::sym(label_time_name) > 0) %>%
      dplyr::select(!!label_time_name) %>%
      unlist() %>%
      unname()

    only_one_tl <- length(unique(label_times)) == 1


    # Figure out what sample details were included and thus need to be pivoted on

    formula_mean <- obj[['metadata']][['averages']][[table_name]][['formula_mean']]
    pivot_columns <- all.vars(formula_mean)
    pivot_columns <- pivot_columns[2:length(pivot_columns)]

    features <- obj[['metadata']][['averages']][[table_name]][['features']]

    # Get rid of columns for parameters users don't want to model
    if(!is.null(avg_params_tokeep) | !is.null(avg_params_todrop)){

      all_cols <- colnames(table)

      # All columns user wants to keep
      cols_keep <- c(features,
                     paste0("mean_", avg_params_tokeep),
                     paste0("sd_", avg_params_tokeep, "_posterior"),
                     paste0("coverage_", avg_params_tokeep))

      # All columns user wants to drop
      cols_drop <- c(paste0("mean_", avg_params_todrop),
                     paste0("sd_", avg_params_todrop, "_posterior"),
                     paste0("coverage_", avg_params_todrop))


      cols_keep <- cols_keep[cols_keep %in% all_cols]
      cols_drop <- cols_drop[cols_drop %in% all_cols]

      # Don't perform filtering if relevant param was NULL
      if(!is.null(avg_params_tokeep)){

        table <- table %>%
          dplyr::select(!!cols_keep)

      }

      if(!is.null(avg_params_todrop)){

        table <- table %>%
          dplyr::select(-!!cols_drop)

      }

      # Need to filter out no longer extant factors
      all_cols <- colnames(table)
      count <- 1
      final_pivot <- c()

      for(p in seq_along(pivot_columns)){

        if(any(grepl(pivot_columns[p], all_cols))){

          final_pivot[count] <- pivot_columns[p]
          count <- count + 1

        }

      }

      pivot_columns <- final_pivot

    }

    # Currently making hard assumption of interactions only...
    # TO-DO: Loosen this restriction
    pattern <- paste0("(.*)_", paste(paste0(pivot_columns, "(.+)"), collapse = ":"))

    # Need to remove _posterior suffix (maybe shouldn't be there at all?)
    # to make general pivoting strategy work
    colnames(table) <- gsub("_posterior$", "", colnames(table))


    # Tidy data!
    tidy_avgs <- table %>%
      tidyr::pivot_longer(
        cols = -!!features,
        names_to = "parameter"
      ) %>%
      dplyr::mutate(
        param_type = dplyr::case_when(
          grepl("^mean", parameter) ~ "mean",
          grepl("^sd", parameter) ~ "sd",
          grepl("^coverage", parameter) ~ "coverage"
        )
      ) %>%
      dplyr::mutate(
        parameter = sub("^(mean_|sd_|coverage_)", "", parameter)
      ) %>%
      tidyr::pivot_wider(names_from = param_type,
                         values_from = "value")

    # Look for substrings corresponding to modeled factors
    substrings <- pivot_columns
    query_strings <- tidy_avgs$parameter

    # Need to extract values of all factors (including label time)
    # associated with a parameter. This assumes that a given parameter
    # will have a set of factors separated by ":". So parameter names
    # are of the form "substr1<value>:substr2<value>:...".
    extract_values <- function(query, substrings) {

      parts <- strsplit(query, ":")[[1]]

      results <- vector("list", length(substrings))
      names(results) <- substrings

      # Look for each possible factor
      for (substr in substrings) {

        matched <- grep(paste0("^", substr), parts, value = TRUE)

        # Return NA if substring is not found
        if (length(matched) > 0) {
          results[[substr]] <- sub(paste0("^", substr), "", matched)
        } else {
          results[[substr]] <- NA
        }
      }

      return(results)
    }

    extracted_values <- lapply(query_strings, extract_values, substrings = substrings)

    # Annotate parameters with factor values
    tidy_avgs <- dplyr::bind_cols(tidy_avgs, dplyr::bind_rows(extracted_values))




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
    lower_bounds <- stats::rnorm(npars,
                          -10,
                          0.01)
    upper_bounds <- stats::rnorm(npars,
                          10,
                          0.01)
    starting_values <- seq(from = -5, to = -2,
                           length.out = npars)


    ### Infer which measured species each row of tidy_avgs belongs to
    feature_mat <- tidy_avgs %>%
      dplyr::select(!!sub_features) %>%
      as.matrix()

    feature_mat <- feature_mat != unassigned_name

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


      valid_cols <- sometimes_assigned[row[!(fm_cols %in% always_assigned)]]

      if(length(valid_cols) > 0){

        valid_cols[1]

      }else{

        always_assigned

      }

    })

    tidy_avgs$measured_specie <- measured_species



    # Metadata columns retained in the final output
    cols_to_group_by <- c(grouping_features,
                       pivot_columns[!(pivot_columns %in% c(label_time_name, sample_feature))])

    # For metadata
    design_factors <- pivot_columns[!(pivot_columns %in% c(label_time_name, sample_feature))]

    # Filter out features that don't have all measured species
    nspecies <- length(unique(measured_species))



    # Impute label time if there is only one
    if(only_one_tl){

      tidy_avgs <- tidy_avgs %>%
        dplyr::mutate(!!label_time_name := unique(label_times))

    }


    # Remove underspecified features
      # Maybe not the move as multi-label times
      # could render some systems identifiable;
      # need to test if this is the case though...
    tidy_avgs <- tidy_avgs %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(cols_to_group_by, label_time_name)))) %>%
      dplyr::filter(dplyr::n() >= (nrow(graph) - 1)) %>%
      dplyr::ungroup()


    # formula list
    if(is.null(modeled_to_measured)){

      species <- rownames(graph)[rownames(graph) != "0"]
      modeled_to_measured <- vector(mode = "list", length = length(species))
      for(s in seq_along(species)){
        modeled_to_measured[[s]] <- stats::as.formula(paste0(species[s], "~", species[s]))
      }

    }


    ### Get scale factors if not provided,
    ### and if sample_feature is specified (as this indicates normalization
    ### across distinct RNA populations)
    if(!is.null(sample_feature) & is.null(scale_factors)){

      species <- colnames(graph)
      species <- species[species != "0"]

      ### Summarise data for groups
      norm_df <- tidy_avgs %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(design_factors, sample_feature, label_time_name)))) %>%
        dplyr::summarise(
          global_logit_fraction = logit(sum(inv_logit(mean)*(10^coverage))/sum(10^coverage)),
          global_lf_sigma = mean(sd)/sqrt(dplyr::n())
        )

      ### Try to normalize
      scale_factors <- tryCatch(
        {
          normalize_EZDynamics(norm_df,
                               modeled_to_measured,
                               sample_feature = sample_feature,
                               species = species,
                               label_time_name = label_time_name,
                               species_to_sf = species_to_sf,
                               design_factors = design_factors[design_factors != sample_feature])
        },
        error = function(e){
          message(paste0("Estimation of normalization scale factors failed!
                         Will only estimate parameters that are identifiable without
                         read counts.
                         Error message thrown during normalization: ", e))
          NULL
        }
      )

    }


    # Don't use reads if can't normalize
    if(is.null(scale_factors) & !is.null(sample_feature)){
      use_coverage <- FALSE
    }


    # Add nuisance name if not multi-compartment modeling
    if(is.null(names(modeled_to_measured))){

      modeled_to_measured <- list(total = modeled_to_measured)
      tidy_avgs <- tidy_avgs %>%
        dplyr::mutate(imputed_sample_feature = 'total')

      sample_feature <- "imputed_sample_feature"

    }

    # Add names to modeled_to_measured to avoid shockingly intensive computation
    # in each optim() function call to find the relevant formula
    modeled_to_measured <- renameFormulas(modeled_to_measured)



    ksyn_index <- unique(graph["0",])
    ksyn_index <- ksyn_index[ksyn_index != 0]

    if(!use_coverage){

      # Don't estimate ksyn if coverage is not being modeled
      npars <- npars - 1


      prior_means <- prior_means[-ksyn_index]
      prior_sds <- prior_sds[-ksyn_index]

      starting_values <- starting_values[-ksyn_index]
      upper_bounds <- upper_bounds[-ksyn_index]
      lower_bounds <- lower_bounds[-ksyn_index]
    }else{
      starting_values[ksyn_index] <- mean(tidy_avgs$coverage, na.rm = TRUE)
      upper_bounds[ksyn_index] <- 15
      lower_bounds[ksyn_index] <- 0
    }


    ### Infer replicate numbers
    meta_groups <- c(label_time_name, unique(c(sample_feature, design_factors)))
    meta_groups <- meta_groups[meta_groups != "imputed_sample_feature"]
    reps <- metadf %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(meta_groups))) %>%
      dplyr::summarise(nreps = dplyr::n())


    # Fit model
    if(is.null(scale_factors)){

      dynfit <- tidy_avgs  %>%
        dplyr::mutate(tl = as.numeric(!!dplyr::sym(label_time_name))) %>%
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
                                                   scale_factor = 1,
                                                   method = "L-BFGS-B",
                                                   hessian = TRUE,
                                                   use_coverage = use_coverage))))

    }else{

      scale_column <- colnames(scale_factors)
      scale_column <- "scale"

      dynfit <- tidy_avgs  %>%
        dplyr::inner_join(scale_factors, by = unique(c(sample_feature, design_factors))) %>%
        dplyr::mutate(tl = as.numeric(!!dplyr::sym(label_time_name))) %>%
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
                                                   scale_factor = !!dplyr::sym(scale_column),
                                                   sample_features = !!dplyr::sym(sample_feature),
                                                   feature_types = measured_specie,
                                                   prior_means = prior_means,
                                                   prior_sds = prior_sds,
                                                   lower = lower_bounds,
                                                   upper = upper_bounds,
                                                   method = "L-BFGS-B",
                                                   hessian = TRUE,
                                                   use_coverage = use_coverage))))

    }



    # Get parameter estimates and uncertainties
    for(n in 1:npars){

      par_name <- parameter_names[n]
      par_se_name <- paste0("se_", par_name)

      dynfit <- dynfit %>% dplyr::mutate(
        !!par_name := purrr::map_dbl(fit, ~ .x$par[n]),
        !!par_se_name := tryCatch(
          {
          purrr::map_dbl(fit, ~ sqrt(abs(diag(solve(.x$hessian))))[n])
          },
          error = function(e) {
            Inf
          }
          ),
      )
    }


  }else{

    # For metadata later
    design_factors <- NULL

    features_to_analyze <- obj[["metadata"]][["fractions"]][[table_name]][["features"]]


    # Determine which column to use for kinetic parameter estimation
    fraction_cols <- colnames(table)

    fraction_of_interest <- fraction_cols[grepl("^logit_fraction_high", fraction_cols)]
    fraction_se <- fraction_cols[grepl("^se_logit_fraction_high", fraction_cols)]


    ### Normalize read counts


    # TO-DO: ALLOW USERS TO JUST USE THE TPM FROM ISOFORM QUANTIFICATION
    reads_norm <- get_normalized_read_counts(obj = obj,
                                             features_to_analyze = features_to_analyze,
                                             fractions_name = table_name,
                                             feature_lengths = feature_lengths)


    ### Prep tables for kinetic parameter estimation

    kinetics <- data.table::setDT(data.table::copy(table))

    # Add label time info
    metadf <- obj$metadf

    metadf <- data.table::setDT(data.table::copy(metadf))
    setkey(metadf, sample)
    setkey(kinetics, sample)

    kinetics <- kinetics[metadf[,c("sample", label_time_name)], nomatch = NULL]

    # Make sure no -s4U controls made it through
    kinetics <- kinetics[tl > 0]


    ### Infer which measured species each row belongs to
    feature_mat <- kinetics %>%
      dplyr::select(!!sub_features) %>%
      as.matrix()

    feature_mat <- feature_mat != unassigned_name

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


      valid_cols <- sometimes_assigned[row[!(fm_cols %in% always_assigned)]]

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


    # Add names to modeled_to_measured to avoid shockingly intensive computation
    # in each optim() function call to find the relevant formula
    modeled_to_measured <- renameFormulas(modeled_to_measured)

    ### Fit model

    # Determine initial vector of parameters and their upper and lower bounds
    npars <- max(graph)

    # Have to add a bit of noise because equal parameters means general solution
    # breaks down as their aren't N eigenvalues. Could generalize to case
    # of equal parameters, but thinking of parameter estimates as continuous
    # random variables that thus have probability of 0 of being equal.
    lower_bounds <- stats::rnorm(npars,
                          -10,
                          0.01)
    upper_bounds <- stats::rnorm(npars,
                          10,
                          0.01)
    starting_values <- stats::rnorm(npars,
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
                                                 prior_means = prior_means,
                                                 prior_sds = prior_sds,
                                                 lower = lower_bounds,
                                                 upper = upper_bounds,
                                                 method = "L-BFGS-B",
                                                 hessian = TRUE,
                                                 use_coverage = TRUE))),
                       n = mean(n) # I'm not sure how to calculate read counts to be passed to AverageAndRegularize(). I think average should be fine though
                       )



    # Get parameter estimates and uncertainties
    fits <- dynfit$fit
    for(n in 1:npars){

      par_name <- parameter_names[n]
      par_se_name <- paste0("se_", par_name)

      dynfit[[par_name]] <- purrr::map_dbl(fits, ~ .x$par[n])
      dynfit[[par_se_name]] <- tryCatch(
        {
          purrr::map_dbl(fits, ~ sqrt(diag(solve(.x$hessian)))[n])
        },
        error = function(e) {
          Inf
        }
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
                                 modeled_to_measured = modeled_to_measured,
                                 graph = graph,
                                 dynamics_design_factors = design_factors,
                                 scale_factors = scale_factors,
                                 feature_lengths = feature_lengths,
                                 parameter = parameter)

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
                               modeled_to_measured = modeled_to_measured,
                               graph = graph, parameter = parameter,
                               dynamics_design_factors = design_factors,
                               scale_factors = scale_factors,
                               feature_lengths = feature_lengths,
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
                                                         modeled_to_measured = modeled_to_measured,
                                                         graph = graph, parameter = parameter,
                                                         dynamics_design_factors = design_factors,
                                                         scale_factors = scale_factors,
                                                         feature_lengths = feature_lengths,
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


# If parameter estimates are nearly identical in generalized ODE
# likelihood, it will crash due to a singularity in the matrix solution.
# This function deals with these edge cases by spacing out nearly duplicate
# values in a vector a small bit
spread_duplicates <- function(vec, digits = 4, spread = 0.01) {

  rounded_vec <- round(vec, digits)
  dup_vals <- unique(rounded_vec[duplicated(rounded_vec)])

  result <- vec

  for (val in dup_vals) {

    indices <- which(rounded_vec == val)

    new_vals <- seq(val - spread, val + spread, length.out = length(indices))

    result[indices] <- new_vals

  }

  return(result)

}

# Likelihood function for generalized dynamical systems modeling
# with averages tables
dynamics_likelihood <- function(parameter_ests, graph, formula_list = NULL,
                                logit_fn, logit_fn_sd, coverage, nreps = 2,
                                tls, sample_features, feature_types,
                                design_factor = NULL,
                                prior_means = rep(-3, times = max(graph)),
                                prior_sds = c(3, rep(1, times = max(graph)-1)),
                                use_coverage = TRUE, alt_coverage = FALSE,
                                coverage_sd = NULL, scale_factor = NULL){



  ### Step 0, check to see if single replicate of data is being passed
  if(all(nreps == 1)){
    single_replicate <- TRUE
  }else{
    single_replicate <- FALSE
  }


  ### Step 1, construct A

  # Make sure parameters are all different values
  count <- 1
  s <- 0.01
  while(length(unique(round(parameter_ests, digits = 4))) != length(parameter_ests)){
    # parameter_ests <- stats::rnorm(length(parameter_ests),
    #                         mean = parameter_ests,
    #                         sd = 0.01)

    parameter_ests <- spread_duplicates(parameter_ests,
                                        spread = s)



    count <- count + 1
    s <- s + 0.01
    if(count > 5){
      stop("Infinite while loop!!")
    }
  }




  # Parameters are on log-scale for ease of optimization
  param_extend <- c(0, exp(parameter_ests))

  # Add back ksyn as dummy value if not modeling it
  if(!use_coverage){
    ksyn_index <- unique(graph["0",])
    ksyn_index <- ksyn_index[ksyn_index != 0]

    param_extend <- append(param_extend, 1, after = ksyn_index)
  }

  param_graph <- matrix(param_extend[graph + 1],
                        nrow = nrow(graph),
                        ncol = ncol(graph),
                        byrow = FALSE)

  # A is same size as graph minus the "0" row and column
  A <- matrix(0,
              nrow = nrow(graph) - 1,
              ncol = ncol(graph) - 1)

  rownames(A) <- rownames(graph[-1,,drop=FALSE])
  colnames(A) <- rownames(graph[-1,,drop=FALSE])



  zero_index <- which(colnames(graph) == "0")

  # Left as an exercise to the reader
  # Just kidding, I'll document this somewhere.
  # TO-DO: Should really compile documentation for all
  # non-obvious mathematical results
  diag(A) <- -rowSums(param_graph[-zero_index,,drop=FALSE])
  A <- A + t(param_graph[-zero_index,-zero_index,drop=FALSE])


  ### Step 2: infer general solution

  tryCatch(
    {
      Rss <- solve(a = A,
                   b = -param_graph[zero_index,-zero_index])

      names(Rss) <- rownames(A)


      ev <- eigen(A)

      lambda <- ev$values
      V<- ev$vectors


      cs <- solve(V, -Rss)

    },
    error = function(x){
      browser()
    }

  )




  ### Step 3: Infer data for actual measured species
  all_ss <- c()
  all_fns <- c()
  indices_to_remove <- c()

  for(n in seq_along(tls)){

    tl <- tls[n]
    sample_feature <- sample_features[n]
    feature_type <- feature_types[n]


    exp_lambda <- exp(lambda*tl)

    scaled_eigenvectors <- V %*% diag(exp_lambda*cs)

    result_vector <- rowSums(scaled_eigenvectors) + Rss

    names(result_vector) <- rownames(A)

    # Evaluate the formulas
    sample_formula <- formula_list[[sample_feature]][[feature_type]]

    # If feature is not in model, throw it out and data for it
    if(is.null(sample_formula)){

      indices_to_remove <- c(indices_to_remove, n)

    }else{



      ### Generalized even to modeling contamination
      # measured_levels <- evaluate_formulas2(result_vector, sample_formula)
      #
      # measured_ss <- evaluate_formulas2(Rss, sample_formula)
      species_to_sum <- all.vars(sample_formula[[3]])

      measured_levels <- sum(result_vector[species_to_sum])
      measured_ss <- sum(Rss[species_to_sum])

      all_fns <- c(all_fns, measured_levels/measured_ss)
      all_ss <- c(all_ss, measured_ss)

    }


  }


  all_fns <- pmax(0.0001, pmin(0.9999, all_fns))


  # Remove data for unmodeled features
  if(length(indices_to_remove) > 0){
    logit_fn <- logit_fn[-indices_to_remove]
    logit_fn_sd <- logit_fn_sd[-indices_to_remove]
  }



  ### Step 4 calculate likelihood
  ll <- stats::dnorm(logit_fn,
              logit(all_fns),
              logit_fn_sd,
              log = TRUE)

  if(use_coverage){

    if(length(indices_to_remove) > 0){
      coverage <- coverage[-indices_to_remove]
      scale_factor <- scale_factor[-indices_to_remove]
      nreps <- nreps[-indices_to_remove]
    }


    if(alt_coverage){

      if(length(indices_to_remove) > 0){
        coverage_sd <- coverage_sd[-indices_to_remove]
      }

      # all_reads <- all_ss*scale_factor

      ll <- ll +
        stats::dnorm(coverage + log10(scale_factor),
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
        stats::dpois(coverage * scale_factor,
                     all_ss,
                     log = TRUE)


    } else{


      # Simple, conservative dispersion parameter estimate of 50
      ll <- ll +
        stats::dnorm(coverage + log10(scale_factor),
                     log10(all_ss),
                     (log(10) * sqrt((1/10^coverage) + 1/100))/sqrt(nreps),
                     log = TRUE)

    }

  }

  # Add prior
  ll <- sum(ll) + sum(stats::dnorm(parameter_ests,
                          mean = prior_means,
                          sd = prior_sds,
                          log = TRUE))

  return(-ll)



}


# TO-DO: Generalize for multi-factor designs
normalize_EZDynamics <- function(norm_df,
                                 modeled_to_measured,
                                 sample_feature,
                                 label_time_name,
                                 species,
                                 design_factors,
                                 species_to_sf = NULL){


  # Hack to deal with devtools::check() NOTEs
  n <- valid <- NULL

  ### What label times are present across all compartments
  valid_label_times <- norm_df %>%
    dplyr::ungroup() %>%
    dplyr::select(!!label_time_name, !!sample_feature, !!design_factors) %>%
    dplyr::distinct() %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(label_time_name, design_factors)))) %>%
    dplyr::count() %>%
    dplyr::mutate(valid = n >= length(modeled_to_measured)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(label_time_name)))) %>%
    dplyr::summarise(valid = all(valid)) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!label_time_name) %>%
    unlist() %>%
    unname() %>%
    unique()

  if(length(valid_label_times) == 0){
    stop("You don't have a set of samples from all sample_features with
         the same label time! Can't estimate scale factors")
  }

  ### Get scale factors for each label time
  scale_list <- list()
  uncertainty_vect <- c()
  count <- 1


  for(t in valid_label_times){

    if(length(design_factors) > 0){

      df_lookup <- norm_df %>%
        dplyr::ungroup() %>%
        dplyr::select(!!design_factors) %>%
        dplyr::distinct()

      scale_df_final <- dplyr::tibble()
      uncertainty <- 0
      for(d in 1:nrow(df_lookup)){

        # Data for particular label time
        norm_df_t <- norm_df %>%
          dplyr::filter(!!dplyr::sym(label_time_name) == t) %>%
          dplyr::inner_join(df_lookup[d,],
                            by = design_factors)

        M <- infer_eqns(norm_df_t = norm_df_t,
                        modeled_to_measured = modeled_to_measured,
                        species = species,
                        species_to_sf = species_to_sf,
                        sample_feature = sample_feature)

        ### Estimate scale factors (and some nuisance parameters)
        # 2*N - 1 parameters; N = number of species:
        # N global length-normalized fractions (nuisance)
        # N - 1 scale factors representing the relative molecular abundances of each specie
        fit <- stats::optim(
          rep(0, times = 2*ncol(M) - 1),
          fn = scale_likelihood,
          M = M,
          y = norm_df_t[['global_logit_fraction']],
          sig = norm_df_t[['global_lf_sigma']],
          method = "L-BFGS-B",
          upper = c(rep(7, times = ncol(M)),
                    rep(5, times = ncol(M) - 1)),
          lower = c(rep(-7, times = ncol(M)),
                    rep(-5, times = ncol(M) - 1)),
          hessian = TRUE
        )


        # Check if singular; if so, go to next label time
        uncertainties <- tryCatch(
          {
            sqrt(diag(solve(fit$hessian)))
          },
          error = function(x){
            Inf
          }

        )

        if(any(is.nan(uncertainties))){
          uncertainties <- Inf
        }

        # Check if any of the parameters are at bounds


        if(!any(is.infinite(uncertainties))){

          uncertainty <- uncertainty + sqrt(sum(uncertainties[(ncol(M)+1):(length(fit$par))]^2))


          # Return estimated scale factors
          scale_factors <- c(1, exp(fit$par[(ncol(M)+1):(length(fit$par))]))
          scale_factors <- as.vector(M %*% scale_factors)
          scale_df <- dplyr::tibble(scale = scale_factors) %>%
            dplyr::bind_cols(df_lookup[d,])
          scale_df[[sample_feature]] <- norm_df_t[[sample_feature]]
          scale_df_final <- scale_df %>%
            dplyr::bind_rows(scale_df_final)



        }else{
          break
        }


      }

      if(d == nrow(df_lookup)){
        # Only increment in this case so that it doesn't increment if normalization
        # fails for one design_factors combination
        scale_list[[count]] <- scale_df_final
        uncertainty_vect[count] <- uncertainty
        count <- count + 1
      }


    }else{

      # Data for particular label time
      norm_df_t <- norm_df %>%
        dplyr::filter(!!dplyr::sym(label_time_name) == t)



      M <- infer_eqns(norm_df_t = norm_df_t,
                      modeled_to_measured = modeled_to_measured,
                      species = species,
                      species_to_sf = species_to_sf,
                      sample_feature = sample_feature)

      ### Estimate scale factors (and some nuisance parameters)
      # 2*N - 1 parameters; N = number of species:
      # N global length-normalized fractions (nuisance)
      # N - 1 scale factors representing the relative molecular abundances of each specie
      fit <- stats::optim(
        rep(0, times = 2*ncol(M) - 1),
        fn = scale_likelihood,
        M = M,
        y = norm_df_t[['global_logit_fraction']],
        sig = norm_df_t[['global_lf_sigma']],
        method = "L-BFGS-B",
        upper = c(rep(7, times = ncol(M)),
                  rep(5, times = ncol(M) - 1)),
        lower = c(rep(-7, times = ncol(M)),
                  rep(-5, times = ncol(M) - 1)),
        hessian = TRUE
      )


      # Check if singular; if so, go to next label time
      uncertainties <- tryCatch(
        {
          sqrt(diag(solve(fit$hessian)))
        },
        error = function(x){
          Inf
        }

      )

      if(any(is.nan(uncertainties))){
        uncertainties <- Inf
      }

      # Check if any of the parameters are at bounds


      if(!any(is.infinite(uncertainties))){

        uncertainty_vect[count] <- sqrt(sum(uncertainties[(ncol(M)+1):(length(fit$par))]^2))


        # Return estimated scale factors
        scale_factors <- c(1, exp(fit$par[(ncol(M)+1):(length(fit$par))]))
        scale_factors <- as.vector(M %*% scale_factors)
        scale_df <- dplyr::tibble(scale = scale_factors)
        scale_df[[sample_feature]] <- norm_df_t[[sample_feature]]
        scale_list[[count]] <- scale_df
        count <- count + 1

      }

    }










  }

  if(count == 1){
    stop("Set of species sampled is insufficient to estimate
         scale factors with!")
  }



  ### Return scale factors with lowest uncertainty
  # Could "average" these, but a bit unsure as to how to
  # do this rigorously
  element_to_keep <- which(uncertainty_vect == min(uncertainty_vect))[1]

  scale_df <- scale_list[[element_to_keep]]

  ### Currently, scale_factors is in terms of modeled species
  ### Need to convert it to a scale factor for each compartment

  return(scale_df)


}


# Infer matrix representation of system of equations for normalization
infer_eqns <- function(norm_df_t,
                       modeled_to_measured,
                       species,
                       species_to_sf,
                       sample_feature){

  # Repeat equations if there are repeat compartments
  actual_mtom <- sapply(norm_df_t[[sample_feature]],
                        function(x) modeled_to_measured[[x]],
                        simplify = FALSE)

  ### Infer system of equations

  # Number of unique sample types
  ns <- length(actual_mtom)


  M <- matrix(0, nrow = length(actual_mtom),
              ncol = length(species))



  # Do we need to map multiple species to one sample_feature?
  multi_mapping <- FALSE
  for(f in seq_along(actual_mtom)){

    if(length(actual_mtom[[f]]) > 1){
      multi_mapping <- TRUE
      break
    }

  }

  # If so, resolve multi-species to one sample_feature
  if(multi_mapping){

    if(is.null(species_to_sf)){
      # List where each element is a vector of species names that
      # belong to a given sample_feature (names of the elements).
      # For example, this might denote that nuclear pre-RNA and nuclear mature
      # RNA come from the "nuclear" sample_feature
      species_to_sf <- list()

      infer_stosf <- TRUE

    }else{
      infer_stosf <- FALSE
    }

    # First pass, figure out which species map to which sample_features
    for(f in seq_along(actual_mtom)){

      newrow <- rep(0, times = length(species))
      names(newrow) <- species
      for(x in seq_along(actual_mtom[[f]])){
        newrow <- newrow + get_coefficients(stats::as.formula(actual_mtom[[f]][[x]]), species)
      }

      M[f,] <- newrow


      if(infer_stosf){

        if(!all(newrow > 0)){

          # Sanity check; make sure if this sample_feature has been seen before,
          # that it provided the same result
          if(length(species_to_sf[[names(actual_mtom)[f]]]) > 0){

            if(!identical(newrow[newrow != 0], species_to_sf[[names(actual_mtom)[f]]])){
              stop("Inferred species groups contradicted one another! Try specifying
                   species_to_sf parameter manually.")
            }

          }else{

            species_to_sf[[names(actual_mtom)[f]]] <- names(newrow[newrow != 0])

          }

        }
      }

    }


    # Make sure that all species have been assigned to some group
    if(!(sum(sapply(species_to_sf, length)) == length(species))){
      stop("Was not able to automatically infer species_to_sf! Try specifying
             it manually.")
    }


    colnames(M) <- species

    M <- sapply(species_to_sf, function(cols){
      rowSums(M[, cols])
    })

    # Normalize (technically maybe unnecessary... Shouldn't hurt though)
    M <- M / rowSums(M)

  }else{

    ### SIMPLE MAPPING OF SPECIES TO SAMPLE_FEATURES

    # f is for formula
    for(f in seq_along(actual_mtom)){

      M[f,] <- get_coefficients(stats::as.formula(actual_mtom[[f]][[1]]), species)

    }

  }


  return(M)

}


# Bit of recursive (dark) magic
parse_expr <- function(expr) {
  if (is.call(expr)) {
    if (expr[[1]] == as.name("+")) {

      # Sum: parse each operand and combine
      lhs <- parse_expr(expr[[2]])
      rhs <- parse_expr(expr[[3]])
      return(c(lhs, rhs))

    } else if (expr[[1]] == as.name("*")) {

      # Product: extract coefficient and variable
      if (is.numeric(expr[[2]])) {

        coeff <- expr[[2]]
        var <- parse_expr(expr[[3]])

      } else if (is.numeric(expr[[3]])) {

        coeff <- expr[[3]]
        var <- parse_expr(expr[[2]])

      } else {

        # Multiplication of variables, not supported
        stop("Unsupported expression format: variables multiplied together")

      }

      var_name <- names(var)
      if (!is.null(var_name)) {

        return(stats::setNames(coeff, var_name))

      } else {

        stop("Unsupported expression format in multiplication")

      }
    } else {

      stop("Unsupported operator: ", expr[[1]])

    }
  } else if (is.name(expr)) {

    # Variable with implicit coefficient 1
    return(stats::setNames(1, as.character(expr)))

  } else if (is.numeric(expr)) {

    # Constant term, ignore it
    return(NULL)

  } else {

    stop("Unsupported expression type")

  }
}


get_coefficients <- function(formula, species) {

  coeffs <- parse_expr(formula[[3]])
  output_vector <- numeric(length(species))
  names(output_vector) <- species
  if (!is.null(coeffs)) {
    output_vector[names(coeffs)] <- coeffs
  }

  return(output_vector)

}


# Simple likelihood function
scale_likelihood <- function(pars, M, y, sig){

  ns <- ncol(M)
  v1 <- inv_logit(pars[1:ns])
  v2 <- c(1, exp(pars[(ns+1):length(pars)]))

  yest <- M%*%(v1*v2) / (M%*%v2)

  yest <- ifelse(
    yest == 1,
    0.99999,
    ifelse(yest == 0,
           0.00001,
           yest))


  ll <- stats::dnorm(y, logit(yest),
              sd = sig,
              log = TRUE)

  return(-sum(ll))


}


# To avoid heavy comps later, conveniently rename
# list elements of modeled_to_measured so that it is
# name of modeled species.
renameFormulas <- function(lst) {

  # Helper to extract the "first variable" used in a formula
  get_first_var <- function(f) {
    vars <- all.vars(f)
    if (length(vars) == 0) {
      stop("No variables found in formula.")
    }
    vars[1]
  }

  # For each named sublist, rename the formula objects
  lapply(lst, function(sublist) {
    # sublist is a list of formula objects
    new_names <- sapply(sublist, get_first_var)
    setNames(sublist, new_names)
  })
}
