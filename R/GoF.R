#' Assess goodness of fit of mixture, kinetic, or dynamical systems models
#'
#' `AssessGoF()` evaluates goodness of fit at either the read-level mutational mixture
#' modeling stage (`type = "fractions"`, assessing `EstimateFractions()`) or the
#' cross-sample kinetic / compartmental ODE modeling stage (`type = "averages"` or
#' `"dynamics"`, assessing `AverageAndRegularize()` or `EZDynamics()`):
#' \itemize{
#'   \item \strong{Mixture model goodness of fit (`type = "fractions"`)}: Evaluates how well
#'   each feature's observed read-level mutational distribution in `cB` (or `EZbakRArrowData`)
#'   agrees with the exact `nT`-conditioned Binomial mixture distribution implied by
#'   `EstimateFractions()` (supporting standard 2-component, hierarchical feature-specific
#'   `pnew`, and general `>2`-component multi-label mixtures, as well as `-label` background
#'   controls). Reports sample-size-invariant Total Variation Distance (`mean_tvd`, `max_tvd`,
#'   representing the fraction of reads displaced into the wrong mutational count bin),
#'   empirical-overdispersion-calibrated likelihood-ratio deviance (`chisq_gof`, `df_gof`,
#'   `pval_gof`, `padj_gof`), cross-sample failure frequency (`n_samples_flagged`,
#'   `frac_samples_flagged`), and mechanistic `dominant_failure_mode` tags. Input filtering
#'   (`min_reads`, `min_samples`, `feature_sample_counts`) determines which features enter the
#'   `AssessGoF()` output tables without modifying `obj$fractions` or any other table.
#'   \item \strong{Kinetic & ODE goodness of fit (`type = "averages"` or `"dynamics"`)}:
#'   Evaluates how well model-predicted sample-by-sample fractions (`fraction_high<muttype>`)
#'   agree with `EstimateFractions()` point estimates and 95\% profile-likelihood confidence
#'   intervals, orthogonally decomposing total discrepancy into \strong{replicate overdispersion}
#'   (`chisq_rep`, `pval_rep`, `padj_rep`) and \strong{structural model lack-of-fit}
#'   (`chisq_lof`, `pval_lof`, `padj_lof`), alongside probability-scale `wrmse_fn`.
#' }
#'
#' Summary statistics per feature are stored in `obj$gof[[output_name]]` and
#' sample-by-sample diagnostics/residuals are stored in `obj$gof_residuals[[output_name]]`
#' (when `store_residuals = TRUE`), both accessible via `EZget(obj, type = "gof")` and
#' `EZget(obj, type = "gof_residuals")`.
#'
#' @param obj An `EZbakRData` (or `EZbakRArrowData`) object containing output from
#' `EstimateFractions()` (and optionally `AverageAndRegularize()` or `EZDynamics()`).
#' @param type Which model class to assess: `"auto"` (default; selects `"dynamics"` or
#' `"averages"` when present and unambiguous, or `"fractions"` when only `EstimateFractions()`
#' has been run on a `cB`-backed object), `"fractions"`, `"averages"`, or `"dynamics"`.
#' @param features Optional feature column vector passed to `EZget()` to select the
#' upstream model table.
#' @param parameter Optional parameter name passed to `EZget()` when `type = "averages"`
#' or `"dynamics"`. When `type = "averages"` and `parameter = NULL`, `AssessGoF()`
#' automatically prioritizes `"log_kdeg"` / `"logit_fraction_*"` tables so co-existing
#' `"log_ksyn"` tables do not cause search ambiguity.
#' @param formula_mean Optional `formula_mean` query passed to `EZget()` for `"averages"`.
#' @param sd_grouping_factors Optional `sd_grouping_factors` query passed to `EZget()`.
#' @param fit_params Optional `fit_params` query passed to `EZget()`.
#' @param kstrat Optional kinetic strategy (`"standard"`, `"shortfeed"`, `"pulse-chase"`)
#' query passed to `EZget()` for `"averages"`.
#' @param populations Optional mutational populations query passed to `EZget()`.
#' @param fraction_design Optional fraction design table passed to `EZget()`.
#' @param sub_features Optional `sub_features` query passed to `EZget()` for `"dynamics"`.
#' @param grouping_features Optional `grouping_features` query passed to `EZget()` for `"dynamics"`.
#' @param sample_feature Optional `sample_feature` query passed to `EZget()` for `"dynamics"`.
#' @param modeled_to_measured Optional `modeled_to_measured` query passed to `EZget()` for `"dynamics"`.
#' @param graph Optional `graph` adjacency matrix query passed to `EZget()` for `"dynamics"`.
#' @param dynamics_design_factors Optional `dynamics_design_factors` query passed to `EZget()`.
#' @param scale_factors Optional `scale_factors` query passed to `EZget()` for `"dynamics"`.
#' @param feature_lengths Optional `feature_lengths` query passed to `EZget()`.
#' @param repeatID Optional numerical repeat ID of the upstream table.
#' @param fractions_repeatID Optional numerical repeat ID of the `EstimateFractions()` table.
#' @param exactMatch Logical; passed to `EZget()` (default `TRUE`).
#' @param fraction_col Optional name of the `fraction_high<muttype>` column in the
#' `fractions` table to evaluate (auto-detected when a single high-mutation population exists).
#' @param min_reads Minimum read count (`n > min_reads`) required for a feature-sample
#' observation to be included in `AssessGoF()`'s output table (without modifying any
#' upstream tables). Defaults to `10` when `type = "fractions"` (matching
#' `AverageAndRegularize()`) and `0` when `type = "averages"` or `"dynamics"`.
#' @param min_samples Minimum number of labeled samples with `n > min_reads` required
#' for a feature to be included in `AssessGoF(..., type = "fractions")`. When `NULL`
#' (default) and `feature_sample_counts = NULL`, requires all labeled samples (matching
#' `AverageAndRegularize()`).
#' @param feature_sample_counts Optional data frame with `<feature columns>` and `nsamps`
#' specifying feature-specific sample-count thresholds for `type = "fractions"`, matching
#' the `feature_sample_counts` convention of `AverageAndRegularize()`.
#' @param min_bin_count Minimum expected read count per mutation-count bin before adaptive
#' tail pooling in `type = "fractions"` (default `5`).
#' @param sample_padj_cutoff Per-sample BH-adjusted p-value threshold used to count
#' `n_samples_flagged` in `type = "fractions"` (default `0.05`).
#' @param sample_tvd_cutoff Per-sample Total Variation Distance threshold used to count
#' `n_samples_flagged` in `type = "fractions"` (default `0.05`, i.e., 5\% displaced reads).
#' @param dispersion_mode How baseline biological dispersion is calibrated:
#' \itemize{
#'   \item For `type = "averages"` / `"dynamics"`: `"shrunk"` / `"empirical"` (default;
#'   Empirical-Bayes trend-moderated replicate dispersion), `"global"`, or
#'   `"measurement_only"` / `"nominal"`.
#'   \item For `type = "fractions"`: `"shrunk"` / `"empirical"` / `"global"` (default;
#'   calibrates deviance against each sample's median per-read divergence across features
#'   so ultra-deep housekeeping genes are not flagged for negligible `pnew` drift), or
#'   `"measurement_only"` / `"nominal"` (strict nominal Binomial mixture test).
#' }
#' @param logit_sd_floor Baseline technical standard deviation floor on the logit fraction
#' scale for kinetic/ODE GoF (default `0.05`).
#' @param include_coverage_gof Logical (default `TRUE`); when assessing an `EZDynamics()`
#' fit that used `use_coverage = TRUE`, also evaluate steady-state normalized read coverage
#' goodness-of-fit (`wrmse_log10_cov`, `chisq_cov`, `pval_lof_joint`, `padj_lof_joint`).
#' @param store_residuals Logical (default `TRUE`); if `TRUE`, stores sample-by-sample
#' diagnostics and residuals in `obj$gof_residuals[[output_name]]` (required for
#' `EZFitPlot()` and sample-level `EZGoFPlot()` modes).
#' @param character_limit Maximum length of automatically generated output table names.
#' @param overwrite Logical (default `TRUE`); whether to overwrite an existing GoF table
#' with identical metadata.
#' @return Input `EZbakRData` object with class `"EZbakRGoF"` prepended, containing
#' feature-level goodness-of-fit metrics in `obj$gof` and (optionally) sample-level
#' diagnostics/residuals in `obj$gof_residuals`.
#' @import data.table
#' @importFrom magrittr %>%
#' @examples
#' simdata <- EZSimulate(30, nreps = 2)
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#' ezbdo <- EstimateFractions(ezbdo)
#' ezbdo <- AssessGoF(ezbdo, type = "fractions")
#' ezbdo <- EstimateKinetics(ezbdo)
#' ezbdo <- AverageAndRegularize(ezbdo)
#' ezbdo <- AssessGoF(ezbdo, type = "averages")
#' @export
AssessGoF <- function(obj,
                      type = c("auto", "fractions", "averages", "dynamics"),
                      features = NULL,
                      parameter = NULL,
                      formula_mean = NULL,
                      sd_grouping_factors = NULL,
                      fit_params = NULL,
                      kstrat = NULL,
                      populations = NULL,
                      fraction_design = NULL,
                      sub_features = NULL,
                      grouping_features = NULL,
                      sample_feature = NULL,
                      modeled_to_measured = NULL,
                      graph = NULL,
                      dynamics_design_factors = NULL,
                      scale_factors = NULL,
                      feature_lengths = NULL,
                      repeatID = NULL,
                      fractions_repeatID = NULL,
                      exactMatch = TRUE,
                      fraction_col = NULL,
                      min_reads = NULL,
                      min_samples = NULL,
                      feature_sample_counts = NULL,
                      min_bin_count = 5,
                      sample_padj_cutoff = 0.05,
                      sample_tvd_cutoff = 0.05,
                      dispersion_mode = c("shrunk", "empirical", "global", "measurement_only", "nominal"),
                      logit_sd_floor = 0.05,
                      include_coverage_gof = TRUE,
                      store_residuals = TRUE,
                      character_limit = 20,
                      overwrite = TRUE){

  if(!methods::is(obj, "EZbakRFractions") || length(obj[["fractions"]]) == 0L){
    stop("`obj` must contain `EstimateFractions()` output in `obj$fractions` before running `AssessGoF()`!")
  }

  type <- match.arg(type)
  dispersion_mode <- match.arg(dispersion_mode)

  if(!is.numeric(logit_sd_floor) || length(logit_sd_floor) != 1L ||
     !is.finite(logit_sd_floor) || logit_sd_floor < 0){
    stop("`logit_sd_floor` must be a single non-negative finite number!")
  }
  if(!is.null(min_reads) && (!is.numeric(min_reads) || length(min_reads) != 1L ||
                             !is.finite(min_reads) || min_reads < 0)){
    stop("`min_reads` must be a single non-negative finite number!")
  }
  if(!is.numeric(min_bin_count) || length(min_bin_count) != 1L ||
     !is.finite(min_bin_count) || min_bin_count <= 0){
    stop("`min_bin_count` must be a single positive finite number!")
  }

  if(!is.null(formula_mean)){
    environment(formula_mean) <- NULL
  }

  ### Step 1: Resolve upstream model type ("fractions", "averages", or "dynamics")
  has_dyn <- length(obj[["metadata"]][["dynamics"]]) > 0L
  has_avg <- length(obj[["metadata"]][["averages"]]) > 0L
  has_cb <- !is.null(obj[["cB"]]) || !is.null(obj[["cBds"]]) || methods::is(obj, "EZbakRArrowData")

  if(type == "auto"){
    dyn_hint <- !is.null(graph) || !is.null(sub_features) ||
      !is.null(grouping_features) || !is.null(modeled_to_measured) ||
      !is.null(dynamics_design_factors)
    avg_hint <- !is.null(formula_mean) || !is.null(kstrat) ||
      !is.null(sd_grouping_factors) || !is.null(fit_params)

    if(has_dyn && !has_avg){
      model_type <- "dynamics"
    }else if(has_avg && !has_dyn){
      model_type <- "averages"
    }else if(has_dyn && has_avg){
      if(dyn_hint && !avg_hint){
        model_type <- "dynamics"
      }else if(avg_hint && !dyn_hint){
        model_type <- "averages"
      }else{
        stop("Both `averages` and `dynamics` tables exist in `obj`! Please specify `type = 'fractions'`, `'averages'`, or `'dynamics'`.")
      }
    }else if(has_cb && length(obj[["mutation_rates"]]) > 0L){
      model_type <- "fractions"
    }else{
      stop("`obj` has neither `averages`/`dynamics` tables (for kinetic GoF) nor read-level `cB` + `mutation_rates` (for mixture GoF)!")
    }
  }else{
    model_type <- type
    if(model_type == "fractions"){
      if(!has_cb || length(obj[["mutation_rates"]]) == 0L){
        stop("`type = 'fractions'` mixture GoF requires read-level mutational counts (`obj$cB` or `EZbakRArrowData`) and `obj$mutation_rates` from `EstimateFractions()`!")
      }
    }else if(model_type == "averages" && !has_avg){
      stop("No `averages` tables found in `obj`! Run `AverageAndRegularize()` first.")
    }else if(model_type == "dynamics" && !has_dyn){
      stop("No `dynamics` tables found in `obj`! Run `EZDynamics()` first.")
    }
  }

  ### Step 2: Locate upstream table and compute GoF
  if(model_type == "fractions"){
    eff_disp <- if(dispersion_mode %in% c("measurement_only", "nominal")) "nominal" else "empirical"
    eff_fxn_rep <- repeatID %||% fractions_repeatID
    always_chk <- !is.null(features) || !is.null(populations) ||
      !is.null(fraction_design) || !is.null(eff_fxn_rep)

    source_table <- EZget(
      obj,
      type = "fractions",
      features = features,
      populations = populations,
      fraction_design = fraction_design,
      repeatID = eff_fxn_rep,
      exactMatch = exactMatch,
      alwaysCheck = always_chk,
      returnNameOnly = TRUE
    )

    if(is.null(source_table) || length(source_table) == 0L){
      stop("No `fractions` table matched your search criteria!")
    }
    if(length(source_table) > 1L){
      stop("More than one `fractions` table matched your search criteria (",
           paste(source_table, collapse = ", "),
           ")! Specify `features`, `populations`, `fraction_design`, or `repeatID` to disambiguate.")
    }

    pred_res <- compute_fractions_gof(
      obj = obj,
      source_table = source_table,
      min_reads = min_reads %||% 10,
      min_samples = min_samples,
      feature_sample_counts = feature_sample_counts,
      min_bin_count = min_bin_count,
      sample_padj_cutoff = sample_padj_cutoff,
      sample_tvd_cutoff = sample_tvd_cutoff,
      dispersion_mode = eff_disp
    )
    eff_disp_store <- eff_disp

  }else if(model_type == "averages"){
    eff_disp <- switch(dispersion_mode,
                       empirical = "shrunk",
                       nominal = "measurement_only",
                       dispersion_mode)
    eff_param <- parameter
    if(is.null(eff_param)){
      avg_metas <- obj[["metadata"]][["averages"]]
      cand_params <- unique(vapply(avg_metas, function(m) as.character(m[["parameter"]]), character(1)))
      gof_ok_params <- cand_params[cand_params == "log_kdeg" | grepl("^logit_fraction_", cand_params)]
      if(length(gof_ok_params) == 0L){
        stop("No GoF-compatible `averages` tables found! Only `log_kdeg` or `logit_fraction_*` models predict sample fractions (`",
             paste(cand_params, collapse = ", "), "` cannot predict fractions directly).")
      }else if(length(gof_ok_params) == 1L){
        eff_param <- gof_ok_params[1L]
      }else if("log_kdeg" %in% gof_ok_params){
        # When a user ran AverageAndRegularize on fractions (for EZDynamics prep) AND on log_kdeg,
        # require explicit `parameter` so we never guess silently between two valid models
        eff_param <- NULL
      }
    }else{
      if(!(eff_param == "log_kdeg" || grepl("^logit_fraction_", eff_param))){
        stop("`parameter = '", eff_param, "'` does not directly predict sample `fraction_high<muttype>` values! ",
             "Use `parameter = 'log_kdeg'` or a `logit_fraction_*` parameter.")
      }
    }

    always_chk <- !is.null(eff_param) || !is.null(features) || !is.null(formula_mean) ||
      !is.null(kstrat) || !is.null(sd_grouping_factors) || !is.null(fit_params) || !is.null(repeatID)

    source_table <- EZget(
      obj,
      type = "averages",
      features = features,
      parameter = eff_param,
      formula_mean = formula_mean,
      sd_grouping_factors = sd_grouping_factors,
      fit_params = fit_params,
      kstrat = kstrat,
      populations = populations,
      fraction_design = fraction_design,
      repeatID = repeatID,
      exactMatch = exactMatch,
      alwaysCheck = always_chk,
      returnNameOnly = TRUE
    )

    if(is.null(source_table) || length(source_table) == 0L){
      stop("No `averages` table matched your search criteria!")
    }
    if(length(source_table) > 1L){
      stop("More than one `averages` table matched your search criteria (",
           paste(source_table, collapse = ", "),
           ")! Specify `parameter`, `features`, `formula_mean`, or `kstrat` to disambiguate.")
    }

    pred_res <- compute_averages_gof(
      obj = obj,
      source_table = source_table,
      populations = populations,
      fraction_design = fraction_design,
      fractions_repeatID = fractions_repeatID,
      fraction_col = fraction_col,
      dispersion_mode = eff_disp,
      logit_sd_floor = logit_sd_floor
    )
    eff_disp_store <- eff_disp

  }else{
    eff_disp <- switch(dispersion_mode,
                       empirical = "shrunk",
                       nominal = "measurement_only",
                       dispersion_mode)
    always_chk <- !is.null(features) || !is.null(sub_features) || !is.null(grouping_features) ||
      !is.null(sample_feature) || !is.null(modeled_to_measured) || !is.null(graph) ||
      !is.null(dynamics_design_factors) || !is.null(parameter) || !is.null(repeatID)

    source_table <- EZget(
      obj,
      type = "dynamics",
      features = features,
      sub_features = sub_features,
      grouping_features = grouping_features,
      sample_feature = sample_feature,
      modeled_to_measured = modeled_to_measured,
      graph = graph,
      dynamics_design_factors = dynamics_design_factors,
      scale_factors = scale_factors,
      feature_lengths = feature_lengths,
      parameter = parameter,
      repeatID = repeatID,
      exactMatch = exactMatch,
      alwaysCheck = always_chk,
      returnNameOnly = TRUE
    )

    if(is.null(source_table) || length(source_table) == 0L){
      stop("No `dynamics` table matched your search criteria!")
    }
    if(length(source_table) > 1L){
      stop("More than one `dynamics` table matched your search criteria (",
           paste(source_table, collapse = ", "),
           ")! Specify `grouping_features`, `graph`, or `repeatID` to disambiguate.")
    }

    pred_res <- compute_dynamics_gof(
      obj = obj,
      source_table = source_table,
      populations = populations,
      fraction_design = fraction_design,
      fractions_repeatID = fractions_repeatID,
      fraction_col = fraction_col,
      dispersion_mode = eff_disp,
      logit_sd_floor = logit_sd_floor,
      include_coverage_gof = include_coverage_gof
    )
    eff_disp_store <- eff_disp
  }

  # Optional explicit min_reads filter on kinetic/ODE GoF tables
  if(model_type != "fractions" && !is.null(min_reads) && min_reads > 0){
    st_f <- pred_res$summary_tbl
    keep_k <- !is.na(st_f$mean_n) & (st_f$mean_n > min_reads)
    if(!any(keep_k)){
      stop("No features made it past `min_reads` filtering!")
    }
    st_f <- st_f[keep_k, , drop = FALSE]
    if("pval_rep" %in% colnames(st_f)){
      ok_r <- !is.na(st_f$pval_rep)
      st_f$padj_rep <- rep(NA_real_, nrow(st_f))
      if(any(ok_r)) st_f$padj_rep[ok_r] <- stats::p.adjust(st_f$pval_rep[ok_r], method = "BH")
    }
    if("pval_lof" %in% colnames(st_f)){
      ok_l <- !is.na(st_f$pval_lof)
      st_f$padj_lof <- rep(NA_real_, nrow(st_f))
      if(any(ok_l)) st_f$padj_lof[ok_l] <- stats::p.adjust(st_f$pval_lof[ok_l], method = "BH")
    }
    if("pval_lof_joint" %in% colnames(st_f)){
      ok_j <- !is.na(st_f$pval_lof_joint)
      st_f$padj_lof_joint <- rep(NA_real_, nrow(st_f))
      if(any(ok_j)) st_f$padj_lof_joint[ok_j] <- stats::p.adjust(st_f$pval_lof_joint[ok_j], method = "BH")
    }
    pred_res$summary_tbl <- st_f
    f_cols_k <- pred_res$metadata$features
    pred_res$residuals_tbl <- dplyr::semi_join(
      pred_res$residuals_tbl,
      st_f[, f_cols_k, drop = FALSE],
      by = f_cols_k
    )
  }

  ### Step 3: Name and store outputs in obj$gof and obj$gof_residuals
  proposed_name <- if(model_type == "fractions"){
    paste0("gof_fxn_", source_table)
  }else{
    paste0("gof_", source_table)
  }
  if(nchar(proposed_name) > character_limit){
    num_gof <- length(obj[["gof"]])
    proposed_name <- paste0("gof", num_gof + 1L)
  }

  gof_meta <- pred_res$metadata

  if(length(obj[["metadata"]][["gof"]]) > 0L){
    output_name <- decide_output(
      obj,
      proposed_name = proposed_name,
      type = "gof",
      source_type = gof_meta$source_type,
      source_table = gof_meta$source_table,
      features = gof_meta$features,
      populations = gof_meta$populations,
      fraction_design = gof_meta$fraction_design,
      parameter = gof_meta$parameter,
      kstrat = gof_meta$kstrat,
      overwrite = overwrite
    )
    if(output_name %in% names(obj[["metadata"]][["gof"]])){
      gof_repeatID <- obj[["metadata"]][["gof"]][[output_name]][["repeatID"]]
    }else{
      gof_repeatID <- length(EZget(
        obj,
        type = "gof",
        source_type = gof_meta$source_type,
        features = gof_meta$features,
        populations = gof_meta$populations,
        fraction_design = gof_meta$fraction_design,
        parameter = gof_meta$parameter,
        kstrat = gof_meta$kstrat,
        formula_mean = gof_meta$formula_mean,
        sub_features = gof_meta$sub_features,
        grouping_features = gof_meta$grouping_features,
        sample_feature = gof_meta$sample_feature,
        returnNameOnly = TRUE,
        exactMatch = TRUE,
        alwaysCheck = TRUE
      )) + 1
    }
  }else{
    output_name <- proposed_name
    gof_repeatID <- 1
  }

  gof_meta$repeatID <- as.numeric(gof_repeatID)
  gof_meta$dispersion_mode <- eff_disp_store
  gof_meta$logit_sd_floor <- logit_sd_floor
  gof_meta$store_residuals <- store_residuals

  obj[["gof"]][[output_name]] <- dplyr::as_tibble(pred_res$summary_tbl)
  obj[["metadata"]][["gof"]][[output_name]] <- gof_meta

  if(store_residuals){
    obj[["gof_residuals"]][[output_name]] <- dplyr::as_tibble(pred_res$residuals_tbl)
    obj[["metadata"]][["gof_residuals"]][[output_name]] <- gof_meta
  }else if(!is.null(obj[["gof_residuals"]][[output_name]])){
    obj[["gof_residuals"]][[output_name]] <- NULL
    obj[["metadata"]][["gof_residuals"]][[output_name]] <- NULL
  }

  if(!methods::is(obj, "EZbakRGoF")){
    class(obj) <- c("EZbakRGoF", class(obj))
  }

  return(obj)
}


# Null-coalescing helper
`%||%` <- function(x, y) if(is.null(x)) y else x


# Resolve target fraction_high<pop> column and companion uncertainty columns
resolve_gof_fraction_cols <- function(fxn_tbl, fraction_col = NULL, parameter_hint = NULL){
  fcols <- colnames(fxn_tbl)
  if(is.null(fraction_col)){
    if(!is.null(parameter_hint) && grepl("^logit_fraction_", parameter_hint)){
      cand <- sub("^logit_", "", parameter_hint)
      if(cand %in% fcols){
        fraction_col <- cand
      }
    }
  }

  if(is.null(fraction_col)){
    all_fc <- fcols[grepl("^fraction_", fcols)]
    nhigh <- lengths(regmatches(all_fc, gregexpr("high", all_fc)))
    cand_fc <- all_fc[nhigh == 1L]
    if(length(cand_fc) == 0L){
      stop("Could not find a `fraction_high<muttype>` column in the `fractions` table!")
    }
    if(length(cand_fc) > 1L){
      stop("Multiple single-label `fraction_high*` columns exist (",
           paste(cand_fc, collapse = ", "),
           ")! Please specify `fraction_col`.")
    }
    fraction_col <- cand_fc[1L]
  }else if(!(fraction_col %in% fcols)){
    stop("Specified `fraction_col = '", fraction_col, "'` does not exist in the `fractions` table!")
  }

  logit_col <- paste0("logit_", fraction_col)
  se_col <- paste0("se_logit_", fraction_col)
  ci_lo_col <- paste0("ci_lower_", fraction_col)
  ci_hi_col <- paste0("ci_upper_", fraction_col)

  if(!(logit_col %in% fcols)){
    stop("Required logit column `", logit_col, "` not found in `fractions` table!")
  }

  # Infer corresponding muttype (e.g. TC in fraction_highTC or fraction_highTC_lowGA)
  m <- regmatches(fraction_col, regexpr("high([A-Za-z]+)", fraction_col))
  muttype <- if(length(m) == 1L && nchar(m) > 4L) sub("^high", "", m) else NULL

  list(
    fraction_col = fraction_col,
    logit_col = logit_col,
    se_col = if(se_col %in% fcols) se_col else NULL,
    ci_lo_col = if(ci_lo_col %in% fcols) ci_lo_col else NULL,
    ci_hi_col = if(ci_hi_col %in% fcols) ci_hi_col else NULL,
    muttype = muttype
  )
}


# Resolve label-time column(s) in metadf for a given muttype and kstrat
resolve_gof_tl_cols <- function(metadf, muttype = NULL, kstrat = NULL, label_time_hint = NULL){
  mcols <- colnames(metadf)
  if(identical(kstrat, "pulse-chase")){
    tp_cand <- if(!is.null(muttype) && paste0("tpulse_", muttype) %in% mcols) paste0("tpulse_", muttype) else "tpulse"
    tc_cand <- if(!is.null(muttype) && paste0("tchase_", muttype) %in% mcols) paste0("tchase_", muttype) else "tchase"
    if(!all(c(tp_cand, tc_cand) %in% mcols)){
      stop("Pulse-chase GoF requires `", tp_cand, "` and `", tc_cand, "` columns in `metadf`!")
    }
    return(list(mode = "pulse-chase", tpulse_col = tp_cand, tchase_col = tc_cand, tl_col = tp_cand))
  }

  if(!is.null(label_time_hint) && label_time_hint %in% mcols){
    return(list(mode = "pulse", tl_col = label_time_hint))
  }
  if(!is.null(muttype) && paste0("tl_", muttype) %in% mcols){
    return(list(mode = "pulse", tl_col = paste0("tl_", muttype)))
  }
  if("tl" %in% mcols){
    return(list(mode = "pulse", tl_col = "tl"))
  }
  tl_any <- mcols[grepl("^tl($|_)", mcols)]
  if(length(tl_any) == 1L){
    return(list(mode = "pulse", tl_col = tl_any[1L]))
  }
  stop("Could not uniquely identify label time column in `obj$metadf`!")
}


# Reconstruct S_labeled x P design matrix X aligned to avg_meta$fit_params
build_averages_design_matrix <- function(metadf_lab, avg_meta){
  fit_params <- avg_meta[["fit_params"]]
  samps <- as.character(metadf_lab$sample)
  dm_saved <- avg_meta[["design_matrix"]]

  if(!is.null(dm_saved) && is.matrix(dm_saved) &&
     all(fit_params %in% colnames(dm_saved)) &&
     all(samps %in% rownames(dm_saved))){
    return(dm_saved[samps, fit_params, drop = FALSE])
  }

  f_mean <- avg_meta[["formula_mean"]]
  param_nm <- avg_meta[["parameter"]] %||% "y"

  if(length(fit_params) == 1L && identical(fit_params, "Intercept")){
    return(matrix(1, nrow = nrow(metadf_lab), ncol = 1L, dimnames = list(samps, fit_params)))
  }

  f_vars <- if(!is.null(f_mean)) setdiff(all.vars(f_mean), param_nm) else character(0)
  if(length(fit_params) == 1L && length(f_vars) <= 1L){
    if(length(f_vars) == 0L || length(unique(as.character(metadf_lab[[f_vars[1L]]]))) <= 1L){
      return(matrix(1, nrow = nrow(metadf_lab), ncol = 1L, dimnames = list(samps, fit_params)))
    }
  }

  df_tmp <- as.data.frame(metadf_lab)
  conv_tl <- if(is.null(avg_meta[["convert_tl_to_factor"]])) TRUE else isTRUE(avg_meta[["convert_tl_to_factor"]])
  tl_vars <- f_vars[grepl("^tl($|_)|^tpulse|^tchase", f_vars)]
  if(conv_tl && length(tl_vars) > 0L){
    for(tv in tl_vars){
      df_tmp[[tv]] <- as.factor(df_tmp[[tv]])
    }
  }
  for(v in f_vars){
    if(v %in% names(df_tmp) && (is.character(df_tmp[[v]]) || is.factor(df_tmp[[v]]))){
      df_tmp[[v]] <- droplevels(as.factor(df_tmp[[v]]))
    }
  }
  df_tmp[[param_nm]] <- 1

  X_raw <- stats::model.matrix(f_mean, df_tmp)
  rownames(X_raw) <- samps

  if(!all(fit_params %in% colnames(X_raw))){
    missing_p <- setdiff(fit_params, colnames(X_raw))
    stop("Could not align design matrix columns to `averages` parameters (missing: ",
         paste(missing_p, collapse = ", "), ")!")
  }

  X_raw[, fit_params, drop = FALSE]
}


# Boundary-safe, profile-CI-aware standardized residuals and effective pseudo-logit observations
compute_gof_obs_metrics <- function(obs_fn,
                                    obs_logit_fn,
                                    se_logit_fn,
                                    ci_lo_fn,
                                    ci_hi_fn,
                                    n_reads,
                                    pred_fn,
                                    pred_logit_fn){
  N <- length(obs_fn)
  if(N == 0L){
    return(list(
      std_resid = numeric(0),
      eta_tilde = numeric(0),
      sigma_tilde = numeric(0),
      prob_sd_dir = numeric(0),
      ci_lo_eff = numeric(0),
      ci_hi_eff = numeric(0)
    ))
  }

  z95 <- stats::qnorm(0.975)
  eps_b <- 1e-6

  pred_fn_c <- pmin(pmax(as.numeric(pred_fn), 1e-7), 1 - 1e-7)
  pred_eta <- as.numeric(pred_logit_fn)
  bad_peta <- !is.finite(pred_eta)
  if(any(bad_peta)){
    pred_eta[bad_peta] <- logit(pred_fn_c[bad_peta])
  }

  obs_fn_num <- pmin(pmax(as.numeric(obs_fn), 0), 1)
  n_safe <- pmax(as.numeric(n_reads), 1)

  # Beta(1,1) fallback logit SE when se_logit_fn is missing or non-finite
  a_fb <- obs_fn_num * n_safe + 1
  b_fb <- (1 - obs_fn_num) * n_safe + 1
  f_fb <- a_fb / (a_fb + b_fb)
  se_fb <- sqrt(a_fb * b_fb / (((a_fb + b_fb)^2) * (a_fb + b_fb + 1))) / pmax(f_fb * (1 - f_fb), 1e-6)

  if(is.null(se_logit_fn)){
    se_num <- se_fb
  }else{
    se_num <- as.numeric(se_logit_fn)
    bad_se <- !is.finite(se_num) | se_num <= 0
    if(any(bad_se)) se_num[bad_se] <- se_fb[bad_se]
  }

  # Impute Wilks/Clopper-Pearson bounds if CI columns are absent in custom EZbakRFractions inputs
  if(is.null(ci_lo_fn) || is.null(ci_hi_fn)){
    eta_fin <- ifelse(obs_fn_num <= eps_b | obs_fn_num >= 1 - eps_b, logit(f_fb), as.numeric(obs_logit_fn))
    bad_ef <- !is.finite(eta_fin)
    if(any(bad_ef)) eta_fin[bad_ef] <- logit(pmin(pmax(obs_fn_num[bad_ef], 1e-6), 1 - 1e-6))

    ci_lo <- inv_logit(eta_fin - z95 * se_num)
    ci_hi <- inv_logit(eta_fin + z95 * se_num)
    at_zero <- obs_fn_num <= eps_b
    at_one <- obs_fn_num >= 1 - eps_b
    if(any(at_zero)){
      ci_lo[at_zero] <- 0
      ci_hi[at_zero] <- 1 - 0.05^(1 / n_safe[at_zero])
    }
    if(any(at_one)){
      ci_hi[at_one] <- 1
      ci_lo[at_one] <- 0.05^(1 / n_safe[at_one])
    }
  }else{
    ci_lo <- pmin(pmax(as.numeric(ci_lo_fn), 0), 1)
    ci_hi <- pmin(pmax(as.numeric(ci_hi_fn), 0), 1)
    na_ci <- is.na(ci_lo) | is.na(ci_hi)
    if(any(na_ci)){
      ci_lo[na_ci] <- 0
      ci_hi[na_ci] <- 1
    }
  }

  is_uninf <- (ci_lo <= eps_b) & (ci_hi >= 1 - eps_b)
  is_low_b <- (!is_uninf) & ((obs_fn_num <= eps_b) | (ci_lo <= eps_b))
  is_high_b <- (!is_uninf) & (!is_low_b) & ((obs_fn_num >= 1 - eps_b) | (ci_hi >= 1 - eps_b))
  is_int <- !(is_uninf | is_low_b | is_high_b)

  z_meas <- numeric(N)
  sigma_tilde <- rep(10, N)
  eta_tilde <- pred_eta
  prob_sd_dir <- rep(0.5, N)

  # 1. One-sided lower-boundary profile interval [0, U_i] (ci_lower_fn == 0 or obs_fn == 0):
  # Avoids penalizing predictions <= EB-shrunk point estimate when 0 is inside the 95% profile CI!
  if(any(is_low_b)){
    idx <- which(is_low_b)
    u_i <- pmax(ci_hi[idx], 1e-8)
    o_i <- obs_fn_num[idx]
    p_i <- pred_fn_c[idx]
    o_mle <- ifelse(o_i <= eps_b, 0, pmin(o_i, 0.25 * u_i))

    z_b <- ifelse(
      p_i <= o_mle,
      0,
      -z95 * sqrt(pmax(p_i - o_mle, 0) / pmax(u_i - o_mle, 1e-12))
    )
    s_prob <- pmax(u_i / z95, 1e-4)
    jac <- pmax(p_i * (1 - p_i), 1e-4)
    sig_pos <- pmax(sqrt(pmax(p_i - o_mle, 1e-12) * pmax(u_i - o_mle, 1e-12)) / (z95 * jac), 0.25)
    sig_zero <- pmax(abs(logit(pmin(u_i, 1 - 1e-7)) - pred_eta[idx]) / z95, 0.5)
    sig_eq <- ifelse(p_i <= o_mle, sig_zero, sig_pos)

    z_meas[idx] <- z_b
    sigma_tilde[idx] <- sig_eq
    eta_tilde[idx] <- pred_eta[idx] + z_b * sig_eq
    prob_sd_dir[idx] <- s_prob
  }

  # 2. One-sided upper-boundary profile interval [L_i, 1] (ci_upper_fn == 1 or obs_fn == 1):
  # Avoids penalizing predictions >= EB-shrunk point estimate when 1 is inside the 95% profile CI!
  if(any(is_high_b)){
    idx <- which(is_high_b)
    l_i <- pmin(ci_lo[idx], 1 - 1e-8)
    o_i <- obs_fn_num[idx]
    p_i <- pred_fn_c[idx]
    o_mle <- ifelse(o_i >= 1 - eps_b, 1, pmax(o_i, 1 - 0.25 * (1 - l_i)))

    z_b <- ifelse(
      p_i >= o_mle,
      0,
      z95 * sqrt(pmax(o_mle - p_i, 0) / pmax(o_mle - l_i, 1e-12))
    )
    s_prob <- pmax((1 - l_i) / z95, 1e-4)
    jac <- pmax(p_i * (1 - p_i), 1e-4)
    sig_pos <- pmax(sqrt(pmax(o_mle - p_i, 1e-12) * pmax(o_mle - l_i, 1e-12)) / (z95 * jac), 0.25)
    sig_zero <- pmax(abs(pred_eta[idx] - logit(pmax(l_i, 1e-7))) / z95, 0.5)
    sig_eq <- ifelse(p_i >= o_mle, sig_zero, sig_pos)

    z_meas[idx] <- z_b
    sigma_tilde[idx] <- sig_eq
    eta_tilde[idx] <- pred_eta[idx] + z_b * sig_eq
    prob_sd_dir[idx] <- s_prob
  }

  # 3. Two-sided interior profile intervals (0 < ci_lower_fn < ci_upper_fn < 1)
  if(any(is_int)){
    idx <- which(is_int)
    o_i <- obs_fn_num[idx]
    p_i <- pred_fn_c[idx]
    lo_i <- ci_lo[idx]
    hi_i <- ci_hi[idx]

    eta_lo <- logit(pmax(lo_i, 1e-7))
    eta_hi <- logit(pmin(hi_i, 1 - 1e-7))
    eta_o <- as.numeric(obs_logit_fn[idx])
    bad_eo <- !is.finite(eta_o)
    if(any(bad_eo)) eta_o[bad_eo] <- logit(o_i[bad_eo])
    # Guard against extreme EB shrinkage pulling point estimate outside unpenalized profile CI
    out_ci <- (eta_o <= eta_lo) | (eta_o >= eta_hi)
    if(any(out_ci)){
      eta_o[out_ci] <- 0.5 * (eta_lo[out_ci] + eta_hi[out_ci])
      o_i[out_ci] <- inv_logit(eta_o[out_ci])
    }

    eta_p <- pred_eta[idx]
    se_i <- se_num[idx]
    diff_eta <- eta_o - eta_p
    pred_above <- eta_p >= eta_o

    s_hi <- pmax((eta_hi - eta_o) / z95, se_i, 1e-4)
    s_lo <- pmax((eta_o - eta_lo) / z95, se_i, 1e-4)
    sig_dir <- ifelse(pred_above, s_hi, s_lo)
    s_prob <- ifelse(pred_above, pmax((hi_i - o_i) / z95, 1e-4), pmax((o_i - lo_i) / z95, 1e-4))

    z_i <- diff_eta / pmax(sig_dir, 1e-8)
    z_meas[idx] <- z_i
    sigma_tilde[idx] <- sig_dir
    eta_tilde[idx] <- eta_o
    prob_sd_dir[idx] <- s_prob
  }

  list(
    std_resid = z_meas,
    eta_tilde = eta_tilde,
    sigma_tilde = sigma_tilde,
    prob_sd_dir = prob_sd_dir,
    ci_lo_eff = ci_lo,
    ci_hi_eff = ci_hi
  )
}


# Vectorized feature-level orthogonal ANOVA GoF summary (Pure Replicate Error + Structural LoF)
summarize_feature_gof <- function(feat_keys_df,
                                  feat_id,
                                  cell_key_str,
                                  obs_fn,
                                  pred_fn,
                                  pred_logit_fn,
                                  n_reads,
                                  obs_metrics,
                                  n_params_vec,
                                  n_at_bounds_vec = NULL,
                                  dispersion_mode = "shrunk",
                                  logit_sd_floor = 0.05){
  F_tot <- nrow(feat_keys_df)
  if(F_tot == 0L){
    return(feat_keys_df)
  }

  if(is.null(n_at_bounds_vec)){
    n_at_bounds_vec <- integer(F_tot)
  }

  # Build unique integer cell_id per (feat_id, experimental_cell)
  combo_str <- paste(feat_id, cell_key_str, sep = "\r")
  u_combos <- unique(combo_str)
  cell_id <- match(combo_str, u_combos)
  first_obs_of_cell <- match(u_combos, combo_str)
  feat_of_cell <- feat_id[first_obs_of_cell]
  pred_eta_cell <- as.numeric(pred_logit_fn[first_obs_of_cell])
  bad_pec <- !is.finite(pred_eta_cell)
  if(any(bad_pec)){
    pred_eta_cell[bad_pec] <- logit(pmin(pmax(pred_fn[first_obs_of_cell[bad_pec]], 1e-7), 1 - 1e-7))
  }

  n_obs_vec <- tabulate(feat_id, nbins = F_tot)
  n_groups_vec <- tabulate(feat_of_cell, nbins = F_tot)
  df_rep_vec <- pmax(n_obs_vec - n_groups_vec, 0L)

  eta_tilde <- obs_metrics$eta_tilde
  sig_tilde <- obs_metrics$sigma_tilde

  # Helper to compute orthogonal (chisq_rep, chisq_lof) for a given observation weight vector
  orth_decomp <- function(w_obs){
    W_cell <- as.numeric(rowsum(w_obs, cell_id, reorder = FALSE))
    num_cell <- as.numeric(rowsum(w_obs * eta_tilde, cell_id, reorder = FALSE))
    eta_bar <- num_cell / pmax(W_cell, 1e-15)

    ss_rep_i <- w_obs * (eta_tilde - eta_bar[cell_id])^2
    ss_rep_f <- numeric(F_tot)
    rs_rep <- rowsum(ss_rep_i, feat_id, reorder = FALSE)
    ss_rep_f[as.integer(rownames(rs_rep))] <- as.numeric(rs_rep)

    ss_lof_c <- W_cell * (eta_bar - pred_eta_cell)^2
    ss_lof_f <- numeric(F_tot)
    rs_lof <- rowsum(ss_lof_c, feat_of_cell, reorder = FALSE)
    ss_lof_f[as.integer(rownames(rs_lof))] <- as.numeric(rs_lof)

    list(rep = ss_rep_f, lof = ss_lof_f)
  }

  # 1. Pure measurement-scale replicate overdispersion test
  w_meas <- 1 / pmax(sig_tilde^2, 1e-12)
  dec_meas <- orth_decomp(w_meas)
  chisq_rep_meas <- dec_meas$rep
  pval_rep <- ifelse(df_rep_vec > 0L,
                     stats::pchisq(chisq_rep_meas, df = df_rep_vec, lower.tail = FALSE),
                     NA_real_)
  padj_rep <- rep(NA_real_, F_tot)
  ok_rep <- !is.na(pval_rep)
  if(any(ok_rep)){
    padj_rep[ok_rep] <- stats::p.adjust(pval_rep[ok_rep], method = "BH")
  }

  # 2. Structural lack-of-fit scale (incorporating logit_sd_floor + dispersion_mode)
  w_base <- 1 / pmax(sig_tilde^2 + logit_sd_floor^2, 1e-12)
  dec_base <- orth_decomp(w_base)

  # Feature-level summary covariates for Empirical Bayes dispersion trend
  rs_cov <- rowsum(as.numeric(n_reads), feat_id, reorder = FALSE)
  mean_cov_vec <- numeric(F_tot)
  mean_cov_vec[as.integer(rownames(rs_cov))] <- as.numeric(rs_cov) / pmax(n_obs_vec[as.integer(rownames(rs_cov))], 1L)

  rs_abs_eta <- rowsum(abs(pred_logit_fn), feat_id, reorder = FALSE)
  mean_abs_eta <- numeric(F_tot)
  mean_abs_eta[as.integer(rownames(rs_abs_eta))] <- as.numeric(rs_abs_eta) / pmax(n_obs_vec[as.integer(rownames(rs_abs_eta))], 1L)

  phi_eff <- rep(1, F_tot)
  d0_prior <- Inf
  idx_r <- which(df_rep_vec >= 1L)

  if(dispersion_mode == "global" && length(idx_r) >= 1L){
    nu_r <- df_rep_vec[idx_r]
    med_null <- stats::qchisq(0.5, df = nu_r) / nu_r
    phi_glob <- max(1, stats::median((dec_base$rep[idx_r] / nu_r) / med_null, na.rm = TRUE))
    if(is.finite(phi_glob)){
      phi_eff <- rep(phi_glob, F_tot)
    }
  }else if(dispersion_mode == "shrunk"){
    if(length(idx_r) >= 10L){
      nu_r <- df_rep_vec[idx_r]
      raw_ratio <- pmax(dec_base$rep[idx_r] / nu_r, 1e-3)
      # Exact bias & sampling variance of log(chisq_nu / nu)
      y_unb <- log(raw_ratio) - (digamma(nu_r / 2) - log(nu_r / 2))
      v_samp <- trigamma(nu_r / 2)

      tr <- fit_ezbakR_se_trend(
        y_logse = y_unb,
        cov_x = log10(pmax(mean_cov_vec[idx_r], 1)),
        mean_x = mean_abs_eta[idx_r]
      )
      v_prior <- max(tr$fit_var - mean(v_samp), 0.05)
      d0_prior <- max(2 / v_prior, 2)

      y_post_r <- (y_unb / v_samp + tr$fitted_values / v_prior) / (1 / v_samp + 1 / v_prior)
      phi_eff[idx_r] <- pmax(exp(y_post_r), 1)

      idx_no_r <- setdiff(seq_len(F_tot), idx_r)
      if(length(idx_no_r) > 0L){
        phi_med <- stats::median(phi_eff[idx_r], na.rm = TRUE)
        phi_eff[idx_no_r] <- max(phi_med, 1)
      }
    }else if(length(idx_r) >= 1L){
      nu_r <- df_rep_vec[idx_r]
      med_null <- stats::qchisq(0.5, df = nu_r) / nu_r
      phi_glob <- max(1, stats::median((dec_base$rep[idx_r] / nu_r) / med_null, na.rm = TRUE))
      if(is.finite(phi_glob)){
        phi_eff <- rep(phi_glob, F_tot)
      }
    }
  }

  chisq_lof <- dec_base$lof / phi_eff
  chisq_total <- (dec_base$rep + dec_base$lof) / phi_eff

  p_free_vec <- pmax(n_params_vec - n_at_bounds_vec, 0L)
  df_lof_vec <- pmax(n_groups_vec - p_free_vec, 0L)
  # Active box-boundary inequality constraint violation even when classical G - P == 0
  bound_viol <- (df_lof_vec == 0L) & (n_at_bounds_vec > 0L) & (chisq_lof > 1e-4)
  if(any(bound_viol)){
    df_lof_vec[bound_viol] <- 1L
  }
  df_total_vec <- df_rep_vec + df_lof_vec
  reduced_chisq <- ifelse(df_total_vec > 0L, chisq_total / df_total_vec, NA_real_)

  pval_lof <- rep(NA_real_, F_tot)
  has_lof <- df_lof_vec > 0L
  if(any(has_lof)){
    if(is.finite(d0_prior) && any(df_rep_vec[has_lof] > 0L)){
      f_stat <- chisq_lof[has_lof] / df_lof_vec[has_lof]
      df2_mod <- df_rep_vec[has_lof] + d0_prior
      pval_lof[has_lof] <- stats::pf(f_stat, df1 = df_lof_vec[has_lof], df2 = df2_mod, lower.tail = FALSE)
    }else{
      pval_lof[has_lof] <- stats::pchisq(chisq_lof[has_lof], df = df_lof_vec[has_lof], lower.tail = FALSE)
    }
  }

  padj_lof <- rep(NA_real_, F_tot)
  if(any(!is.na(pval_lof))){
    ok_l <- !is.na(pval_lof)
    padj_lof[ok_l] <- stats::p.adjust(pval_lof[ok_l], method = "BH")
  }

  # 3. Probability-scale effect size metrics (wrmse_fn, mae_fn, max_abs_resid_fn, mean_signed_resid_fn)
  resid_fn <- as.numeric(obs_fn) - as.numeric(pred_fn)
  w_prob <- 1 / (obs_metrics$prob_sd_dir^2 + 0.01^2)

  rs_w <- rowsum(w_prob, feat_id, reorder = FALSE)
  rs_wr2 <- rowsum(w_prob * (resid_fn^2), feat_id, reorder = FALSE)
  rs_wr1 <- rowsum(w_prob * resid_fn, feat_id, reorder = FALSE)
  rs_abs <- rowsum(abs(resid_fn), feat_id, reorder = FALSE)

  wrmse_fn <- rep(NA_real_, F_tot)
  mean_signed_resid_fn <- rep(NA_real_, F_tot)
  mae_fn <- rep(NA_real_, F_tot)
  f_present <- as.integer(rownames(rs_w))
  wrmse_fn[f_present] <- sqrt(as.numeric(rs_wr2) / pmax(as.numeric(rs_w), 1e-15))
  mean_signed_resid_fn[f_present] <- as.numeric(rs_wr1) / pmax(as.numeric(rs_w), 1e-15)
  mae_fn[f_present] <- as.numeric(rs_abs) / pmax(n_obs_vec[f_present], 1L)

  # Fast per-feature max(|resid_fn|)
  ord_abs <- order(feat_id, -abs(resid_fn))
  first_in_feat <- !duplicated(feat_id[ord_abs])
  max_abs_resid_fn <- rep(NA_real_, F_tot)
  max_abs_resid_fn[feat_id[ord_abs][first_in_feat]] <- abs(resid_fn[ord_abs][first_in_feat])

  aic_vec <- chisq_total + 2 * n_params_vec
  bic_vec <- chisq_total + log(pmax(n_obs_vec, 1L)) * n_params_vec

  out <- dplyr::bind_cols(
    feat_keys_df,
    dplyr::tibble(
      n_obs = as.integer(n_obs_vec),
      n_groups = as.integer(n_groups_vec),
      n_params = as.integer(n_params_vec),
      n_at_bounds = as.integer(n_at_bounds_vec),
      wrmse_fn = wrmse_fn,
      mae_fn = mae_fn,
      max_abs_resid_fn = max_abs_resid_fn,
      mean_signed_resid_fn = mean_signed_resid_fn,
      chisq_total = chisq_total,
      df_total = as.integer(df_total_vec),
      reduced_chisq = reduced_chisq,
      chisq_lof = chisq_lof,
      df_lof = as.integer(df_lof_vec),
      pval_lof = pval_lof,
      padj_lof = padj_lof,
      chisq_rep = chisq_rep_meas,
      df_rep = as.integer(df_rep_vec),
      pval_rep = pval_rep,
      padj_rep = padj_rep,
      aic = aic_vec,
      bic = bic_vec,
      dispersion_factor = phi_eff,
      mean_coverage = mean_cov_vec
    )
  )

  out
}


# Forward GoF engine for AverageAndRegularize outputs
compute_averages_gof <- function(obj,
                                 source_table,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 fractions_repeatID = NULL,
                                 fraction_col = NULL,
                                 dispersion_mode = "shrunk",
                                 logit_sd_floor = 0.05){
  avg_tbl <- obj[["averages"]][[source_table]]
  avg_meta <- obj[["metadata"]][["averages"]][[source_table]]
  feat_cols <- avg_meta[["features"]]
  param_nm <- avg_meta[["parameter"]]
  fit_params <- avg_meta[["fit_params"]]
  k_strat <- avg_meta[["kstrat"]]

  is_direct_fn <- grepl("^logit_fraction_", param_nm)
  if(!is_direct_fn){
    if(is.null(k_strat)){
      src_kin <- avg_meta[["source_table"]]
      if(!is.null(src_kin) && !is.null(obj[["metadata"]][["kinetics"]][[src_kin]])){
        k_strat <- obj[["metadata"]][["kinetics"]][[src_kin]][["kstrat"]]
      }else{
        kin_cands <- EZget(obj, type = "kinetics", features = feat_cols,
                           exactMatch = TRUE, returnNameOnly = TRUE)
        if(length(kin_cands) == 1L){
          k_strat <- obj[["metadata"]][["kinetics"]][[kin_cands]][["kstrat"]]
        }
      }
    }
    k_strat <- k_strat %||% "standard"
    if(!(k_strat %in% c("standard", "shortfeed", "pulse-chase"))){
      stop("`AssessGoF()` supports `kstrat` in c('standard', 'shortfeed', 'pulse-chase'), not '",
           k_strat, "'.")
    }
  }

  fractions_name <- EZget(
    obj,
    type = "fractions",
    features = feat_cols,
    populations = populations %||% avg_meta[["populations"]],
    fraction_design = fraction_design %||% avg_meta[["fraction_design"]],
    repeatID = fractions_repeatID,
    exactMatch = TRUE,
    returnNameOnly = TRUE
  )
  if(is.null(fractions_name) || length(fractions_name) != 1L){
    stop("Could not uniquely locate the upstream `fractions` table for `", source_table, "`!")
  }

  fxn_tbl <- obj[["fractions"]][[fractions_name]]
  fcol_info <- resolve_gof_fraction_cols(fxn_tbl, fraction_col = fraction_col, parameter_hint = param_nm)
  tl_info <- resolve_gof_tl_cols(obj$metadf, muttype = fcol_info$muttype, kstrat = k_strat)

  metadf <- dplyr::as_tibble(obj$metadf)
  if(tl_info$mode == "pulse-chase"){
    metadf_lab <- metadf[as.numeric(metadf[[tl_info$tpulse_col]]) > 0, , drop = FALSE]
  }else{
    metadf_lab <- metadf[as.numeric(metadf[[tl_info$tl_col]]) > 0, , drop = FALSE]
  }
  if(nrow(metadf_lab) == 0L){
    stop("No labeled samples found in `obj$metadf`!")
  }

  X <- build_averages_design_matrix(metadf_lab, avg_meta)
  samps_lab <- rownames(X)
  S_lab <- length(samps_lab)

  mean_cols <- paste0("mean_", fit_params)
  if(!all(mean_cols %in% colnames(avg_tbl))){
    stop("Expected parameter columns missing in `averages` table `", source_table, "`!")
  }

  B_mat <- as.matrix(avg_tbl[, mean_cols, drop = FALSE])
  # Handle any NA coefficients cleanly (e.g. if singular design left NAs)
  B_clean <- B_mat
  B_clean[!is.finite(B_clean)] <- 0
  L_hat <- tcrossprod(B_clean, X) # F_avg x S_lab

  if(is_direct_fn){
    pred_logit_mat <- L_hat
    pred_fn_mat <- inv_logit(pred_logit_mat)
    tl_vec <- as.numeric(metadf_lab[[tl_info$tl_col]])
    cell_sig_samps <- apply(cbind(round(X, 8), tl_vec), 1L, paste, collapse = "|")
  }else if(k_strat == "standard"){
    tl_vec <- as.numeric(metadf_lab[[tl_info$tl_col]])
    kt_mat <- exp(L_hat) * rep(tl_vec, each = nrow(L_hat))
    pred_fn_mat <- -expm1(-pmax(kt_mat, 0))
    pred_logit_mat <- logit(pmin(pmax(pred_fn_mat, 1e-7), 1 - 1e-7))
    cell_sig_samps <- apply(cbind(round(X, 8), tl_vec), 1L, paste, collapse = "|")
  }else if(k_strat == "shortfeed"){
    tl_vec <- as.numeric(metadf_lab[[tl_info$tl_col]])
    kt_mat <- exp(L_hat) * rep(tl_vec, each = nrow(L_hat))
    pred_fn_mat <- pmin(pmax(kt_mat, 1e-7), 1 - 1e-7)
    pred_logit_mat <- logit(pred_fn_mat)
    cell_sig_samps <- apply(cbind(round(X, 8), tl_vec), 1L, paste, collapse = "|")
  }else{
    tp_vec <- as.numeric(metadf_lab[[tl_info$tpulse_col]])
    tc_vec <- as.numeric(metadf_lab[[tl_info$tchase_col]])
    tl_vec <- tp_vec + tc_vec
    k_mat <- exp(L_hat)
    tp_mat <- k_mat * rep(tp_vec, each = nrow(k_mat))
    tc_mat <- k_mat * rep(tc_vec, each = nrow(k_mat))
    pred_fn_mat <- (-expm1(-pmax(tp_mat, 0))) * exp(-pmax(tc_mat, 0))
    pred_logit_mat <- logit(pmin(pmax(pred_fn_mat, 1e-7), 1 - 1e-7))
    cell_sig_samps <- apply(cbind(round(X, 8), tp_vec, tc_vec), 1L, paste, collapse = "|")
  }

  # Match observed rows in fxn_tbl to (feature_row, sample_col)
  feat_keys_df <- dplyr::as_tibble(avg_tbl[, feat_cols, drop = FALSE])
  F_tot <- nrow(feat_keys_df)

  fxn_sub <- fxn_tbl[fxn_tbl$sample %in% samps_lab, , drop = FALSE]
  avg_key_df <- feat_keys_df
  avg_key_df$.ez_feat_id <- seq_len(F_tot)

  joined <- dplyr::inner_join(dplyr::as_tibble(fxn_sub), avg_key_df, by = feat_cols)
  if(nrow(joined) == 0L){
    stop("No overlapping `(sample, feature)` observations found between `", source_table,
         "` and `fractions` table `", fractions_name, "`!")
  }

  samp_col_idx <- match(as.character(joined$sample), samps_lab)
  feat_row_idx <- joined$.ez_feat_id
  lin_idx <- feat_row_idx + (samp_col_idx - 1L) * F_tot

  pred_fn_obs <- pred_fn_mat[lin_idx]
  pred_logit_obs <- pred_logit_mat[lin_idx]

  obs_fn <- joined[[fcol_info$fraction_col]]
  obs_logit <- joined[[fcol_info$logit_col]]
  se_logit <- if(!is.null(fcol_info$se_col)) joined[[fcol_info$se_col]] else NULL
  ci_lo <- if(!is.null(fcol_info$ci_lo_col)) joined[[fcol_info$ci_lo_col]] else NULL
  ci_hi <- if(!is.null(fcol_info$ci_hi_col)) joined[[fcol_info$ci_hi_col]] else NULL

  obs_metrics <- compute_gof_obs_metrics(
    obs_fn = obs_fn,
    obs_logit_fn = obs_logit,
    se_logit_fn = se_logit,
    ci_lo_fn = ci_lo,
    ci_hi_fn = ci_hi,
    n_reads = joined$n,
    pred_fn = pred_fn_obs,
    pred_logit_fn = pred_logit_obs
  )

  # Effective rank P_f per feature (exact qr(X)$rank when complete; per-subset rank if missing samples)
  rank_full <- qr(X)$rank
  obs_counts <- tabulate(feat_row_idx, nbins = F_tot)
  if(all(obs_counts == S_lab)){
    n_params_vec <- rep(rank_full, F_tot)
  }else{
    n_params_vec <- rep(rank_full, F_tot)
    inc_feats <- which(obs_counts < S_lab & obs_counts > 0L)
    if(length(inc_feats) > 0L){
      s_by_f <- split(samp_col_idx, feat_row_idx)
      for(fi in inc_feats){
        s_i <- unique(s_by_f[[as.character(fi)]])
        n_params_vec[fi] <- qr(X[s_i, , drop = FALSE])$rank
      }
    }
  }

  summary_tbl <- summarize_feature_gof(
    feat_keys_df = feat_keys_df,
    feat_id = feat_row_idx,
    cell_key_str = cell_sig_samps[samp_col_idx],
    obs_fn = obs_fn,
    pred_fn = pred_fn_obs,
    pred_logit_fn = pred_logit_obs,
    n_reads = joined$n,
    obs_metrics = obs_metrics,
    n_params_vec = n_params_vec,
    n_at_bounds_vec = integer(F_tot),
    dispersion_mode = dispersion_mode,
    logit_sd_floor = logit_sd_floor
  )

  meta_extra_cols <- setdiff(colnames(metadf_lab), "sample")
  meta_attach <- metadf_lab[samp_col_idx, meta_extra_cols, drop = FALSE]

  residuals_tbl <- dplyr::bind_cols(
    dplyr::tibble(sample = as.character(joined$sample)),
    joined[, feat_cols, drop = FALSE],
    dplyr::tibble(tl = tl_vec[samp_col_idx]),
    meta_attach[, setdiff(colnames(meta_attach), "tl"), drop = FALSE],
    dplyr::tibble(
      obs_fn = as.numeric(obs_fn),
      pred_fn = as.numeric(pred_fn_obs),
      resid_fn = as.numeric(obs_fn) - as.numeric(pred_fn_obs),
      obs_logit_fn = as.numeric(obs_logit),
      pred_logit_fn = as.numeric(pred_logit_obs),
      resid_logit_fn = as.numeric(obs_logit) - as.numeric(pred_logit_obs),
      std_resid = obs_metrics$std_resid,
      ci_lower_fn = obs_metrics$ci_lo_eff,
      ci_upper_fn = obs_metrics$ci_hi_eff,
      n = as.numeric(joined$n)
    )
  )

  list(
    summary_tbl = summary_tbl,
    residuals_tbl = residuals_tbl,
    metadata = list(
      source_type = "averages",
      source_table = source_table,
      fractions_table = fractions_name,
      fraction_col = fcol_info$fraction_col,
      tl_info = tl_info,
      features = feat_cols,
      parameter = param_nm,
      formula_mean = avg_meta[["formula_mean"]],
      sd_grouping_factors = avg_meta[["sd_grouping_factors"]],
      fit_params = fit_params,
      kstrat = if(is_direct_fn) NULL else k_strat,
      populations = avg_meta[["populations"]],
      fraction_design = avg_meta[["fraction_design"]]
    )
  )
}


# Assign measured_specie column identically to EZDynamics()
infer_measured_species <- function(df, sub_features, unassigned_name = "__no_feature"){
  feature_mat <- as.matrix(df[, sub_features, drop = FALSE]) != unassigned_name
  fm_cols <- colnames(feature_mat)
  always_assigned <- fm_cols[which(colSums(feature_mat) == nrow(feature_mat))]

  if(length(always_assigned) > 1L){
    stop("Multiple `sub_features` never equal `unassigned_name`; cannot infer `measured_specie`!")
  }else if(length(always_assigned) == 0L){
    sometimes_assigned <- fm_cols
  }else{
    sometimes_assigned <- fm_cols[fm_cols != always_assigned]
  }

  sub_mat <- feature_mat[, !(fm_cols %in% always_assigned), drop = FALSE]
  if(ncol(sub_mat) == 0L){
    rep(always_assigned, nrow(feature_mat))
  }else{
    has_any <- rowSums(sub_mat) > 0L
    first_col <- max.col(sub_mat, ties.method = "first")
    ms <- character(nrow(feature_mat))
    ms[has_any] <- sometimes_assigned[first_col[has_any]]
    if(any(!has_any)){
      ms[!has_any] <- always_assigned
    }
    ms
  }
}


# General stoichiometry-aware forward ODE evaluator for arbitrary (sample_feature, measured_specie, tl)
evaluate_ode_forward <- function(par_vec,
                                 graph,
                                 formula_list,
                                 tls,
                                 sample_features,
                                 feature_types,
                                 use_coverage = TRUE){
  n_g <- nrow(graph)
  zero_index <- which(colnames(graph) == "0")
  non_zero_idx <- seq_len(n_g)[-zero_index]
  species_names <- rownames(graph)[non_zero_idx]
  K <- length(species_names)

  ksyn_index <- unique(graph["0", ])
  ksyn_index <- ksyn_index[ksyn_index != 0]
  max_par <- max(graph)

  if(!use_coverage){
    full_slot_to_vec <- integer(max_par + 1L)
    full_slot_to_vec[1L] <- 1L
    rem_slots <- setdiff(seq_len(max_par), ksyn_index)
    full_slot_to_vec[ksyn_index + 1L] <- 2L
    full_slot_to_vec[rem_slots + 1L] <- seq_along(rem_slots) + 2L
  }else{
    full_slot_to_vec <- c(1L, seq_len(max_par) + 2L)
  }
  graph_map <- full_slot_to_vec[as.integer(graph) + 1L]

  param_vec <- c(0, 1, exp(as.numeric(par_vec)))
  param_graph <- matrix(param_vec[graph_map], nrow = n_g, ncol = n_g)

  sub_pg <- param_graph[non_zero_idx, , drop = FALSE]
  A <- t(sub_pg[, non_zero_idx, drop = FALSE])
  diag(A) <- diag(A) - rowSums(sub_pg)
  b_vec <- -param_graph[zero_index, non_zero_idx]

  Rss <- tryCatch(solve(A, b_vec), error = function(e) NULL)
  n_obs <- length(tls)
  if(is.null(Rss) || any(!is.finite(Rss)) || any(Rss <= 0)){
    return(list(
      pred_fn = rep(NA_real_, n_obs),
      pred_logit_fn = rep(NA_real_, n_obs),
      pred_ss = rep(NA_real_, n_obs),
      keep_mask = rep(FALSE, n_obs)
    ))
  }

  keep_mask <- logical(n_obs)
  C_rows <- vector("list", n_obs)
  for(i in seq_len(n_obs)){
    sf <- sample_features[i]
    ft <- feature_types[i]
    f_obj <- formula_list[[sf]][[ft]]
    if(!is.null(f_obj)){
      keep_mask[i] <- TRUE
      C_rows[[i]] <- get_coefficients(stats::as.formula(f_obj), species_names)
    }
  }

  keep_idx <- which(keep_mask)
  pred_fn <- rep(NA_real_, n_obs)
  pred_logit <- rep(NA_real_, n_obs)
  pred_ss <- rep(NA_real_, n_obs)
  if(length(keep_idx) == 0L){
    return(list(pred_fn = pred_fn, pred_logit_fn = pred_logit, pred_ss = pred_ss, keep_mask = keep_mask))
  }

  kept_tls <- as.numeric(tls[keep_idx])
  u_tls <- unique(kept_tls)
  tl_map <- match(kept_tls, u_tls)

  R_new_mat <- matrix(0, nrow = K, ncol = length(u_tls))
  for(u in seq_len(length(u_tls))){
    if(u_tls[u] <= 0){
      R_new_mat[, u] <- 0
    }else{
      E_u <- tryCatch(expm_ode(u_tls[u] * A), error = function(e) NULL)
      if(is.null(E_u) || any(!is.finite(E_u))){
        return(list(
          pred_fn = rep(NA_real_, n_obs),
          pred_logit_fn = rep(NA_real_, n_obs),
          pred_ss = rep(NA_real_, n_obs),
          keep_mask = keep_mask
        ))
      }
      R_new_mat[, u] <- pmax(Rss - as.vector(E_u %*% Rss), 0)
    }
  }

  C_mat <- do.call(rbind, C_rows[keep_idx])
  ss_kept <- as.vector(C_mat %*% Rss)
  if(K == 1L){
    new_kept <- C_mat[, 1L] * R_new_mat[1L, tl_map]
  }else{
    new_kept <- rowSums(C_mat * t(R_new_mat[, tl_map, drop = FALSE]))
  }

  fn_raw <- pmin(pmax(new_kept / pmax(ss_kept, 1e-15), 0), 1)
  fn_clamped <- pmin(pmax(fn_raw, 1e-4), 1 - 1e-4)

  pred_fn[keep_idx] <- fn_raw
  pred_logit[keep_idx] <- ifelse(kept_tls <= 0, -Inf, logit(fn_clamped))
  pred_ss[keep_idx] <- ss_kept

  list(
    pred_fn = pred_fn,
    pred_logit_fn = pred_logit,
    pred_ss = pred_ss,
    keep_mask = keep_mask
  )
}


# Forward GoF engine for EZDynamics outputs
compute_dynamics_gof <- function(obj,
                                 source_table,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 fractions_repeatID = NULL,
                                 fraction_col = NULL,
                                 dispersion_mode = "shrunk",
                                 logit_sd_floor = 0.05,
                                 include_coverage_gof = TRUE){
  dyn_tbl <- obj[["dynamics"]][[source_table]]
  dyn_meta <- obj[["metadata"]][["dynamics"]][[source_table]]

  graph <- dyn_meta[["graph"]]
  sub_features <- dyn_meta[["sub_features"]]
  grouping_features <- dyn_meta[["grouping_features"]]
  sample_feature <- dyn_meta[["sample_feature"]]
  formula_list <- dyn_meta[["modeled_to_measured"]]
  design_factors <- dyn_meta[["dynamics_design_factors"]]
  scale_df <- dyn_meta[["scale_factors"]]
  feat_cols <- dyn_meta[["features"]]
  unassigned_name <- dyn_meta[["unassigned_name"]] %||% "__no_feature"

  # Infer parameter_names and use_coverage cleanly for both new and legacy objects
  param_names <- dyn_meta[["parameter_names"]]
  if(is.null(param_names)){
    dcols <- colnames(dyn_tbl)
    param_names <- dcols[paste0("se_", dcols) %in% dcols]
  }
  ksyn_index <- unique(graph["0", ])
  ksyn_index <- ksyn_index[ksyn_index != 0]
  use_cov <- dyn_meta[["use_coverage"]]
  if(is.null(use_cov)){
    use_cov <- (length(param_names) == max(graph))
  }

  from_fractions <- ("sample" %in% colnames(dyn_tbl)) || identical(dyn_meta[["source_type"]], "fractions")
  dyn_group_cols <- if(from_fractions) c("sample", grouping_features) else c(grouping_features, design_factors)

  # Locate upstream fractions table
  if(from_fractions && !is.null(dyn_meta[["source_table"]]) &&
     dyn_meta[["source_table"]] %in% names(obj[["fractions"]]) &&
     is.null(fractions_repeatID)){
    fractions_name <- dyn_meta[["source_table"]]
  }else{
    fractions_name <- EZget(
      obj,
      type = "fractions",
      features = feat_cols,
      populations = populations %||% dyn_meta[["populations"]],
      fraction_design = fraction_design %||% dyn_meta[["fraction_design"]],
      repeatID = fractions_repeatID,
      exactMatch = TRUE,
      returnNameOnly = TRUE
    )
  }
  if(is.null(fractions_name) || length(fractions_name) != 1L){
    stop("Could not uniquely locate the upstream `fractions` table for `", source_table, "`!")
  }

  fxn_tbl <- dplyr::as_tibble(obj[["fractions"]][[fractions_name]])
  fcol_info <- resolve_gof_fraction_cols(fxn_tbl, fraction_col = fraction_col, parameter_hint = dyn_meta[["parameter"]])
  tl_info <- resolve_gof_tl_cols(
    obj$metadf,
    muttype = fcol_info$muttype,
    label_time_hint = dyn_meta[["label_time_name"]]
  )

  metadf_lab <- dplyr::as_tibble(obj$metadf)
  metadf_lab <- metadf_lab[as.numeric(metadf_lab[[tl_info$tl_col]]) > 0, , drop = FALSE]
  if(nrow(metadf_lab) == 0L){
    stop("No labeled samples found in `obj$metadf`!")
  }

  is_imputed_sf <- identical(sample_feature, "imputed_sample_feature")
  if(is_imputed_sf){
    metadf_lab$imputed_sample_feature <- "total"
  }else{
    metadf_lab[[sample_feature]] <- as.character(metadf_lab[[sample_feature]])
  }
  if(length(design_factors) > 0L){
    for(df_col in design_factors){
      metadf_lab[[df_col]] <- as.character(metadf_lab[[df_col]])
    }
  }

  # Attach metadf + measured_specie to labeled fractions rows
  fxn_lab <- dplyr::inner_join(fxn_tbl, metadf_lab, by = "sample")
  if(nrow(fxn_lab) == 0L){
    stop("No labeled observations found in `fractions` table `", fractions_name, "`!")
  }
  fxn_lab$measured_specie <- infer_measured_species(fxn_lab, sub_features = sub_features, unassigned_name = unassigned_name)

  # Optional normalized coverage lookup
  do_cov_gof <- isTRUE(include_coverage_gof) && isTRUE(use_cov)
  if(do_cov_gof){
    if(from_fractions){
      # EZDynamics(type = "fractions") fits steady-state abundances directly to raw `n`
      fxn_lab$.ez_obs_cov <- as.numeric(fxn_lab$n)
    }else{
      rc_dt <- get_normalized_read_counts(
        obj = obj,
        features_to_analyze = feat_cols,
        fractions_name = fractions_name,
        feature_lengths = dyn_meta[["feature_lengths"]]
      )
      rc_tbl <- dplyr::as_tibble(rc_dt)[, c("sample", feat_cols, "normalized_reads"), drop = FALSE]
      fxn_lab <- dplyr::inner_join(fxn_lab, rc_tbl, by = c("sample", feat_cols))
      if(!is.null(scale_df)){
        join_sf_cols <- intersect(colnames(scale_df), colnames(fxn_lab))
        sd_clean <- dplyr::as_tibble(scale_df)
        for(jc in join_sf_cols) sd_clean[[jc]] <- as.character(sd_clean[[jc]])
        fxn_lab <- dplyr::inner_join(fxn_lab, sd_clean, by = join_sf_cols)
        fxn_lab$.ez_obs_cov <- as.numeric(fxn_lab$normalized_reads) * as.numeric(fxn_lab$scale)
      }else{
        fxn_lab$.ez_obs_cov <- as.numeric(fxn_lab$normalized_reads)
      }
    }
  }

  # Match observations to rows of dyn_tbl
  dyn_keys_df <- dplyr::as_tibble(dyn_tbl[, dyn_group_cols, drop = FALSE])
  for(dc in dyn_group_cols){
    dyn_keys_df[[dc]] <- as.character(dyn_keys_df[[dc]])
    fxn_lab[[dc]] <- as.character(fxn_lab[[dc]])
  }
  F_tot <- nrow(dyn_keys_df)
  dyn_keys_id <- dyn_keys_df
  dyn_keys_id$.ez_dyn_id <- seq_len(F_tot)

  joined <- dplyr::inner_join(fxn_lab, dyn_keys_id, by = dyn_group_cols)
  if(nrow(joined) == 0L){
    stop("No matching observations found between `", source_table, "` and `fractions` table `", fractions_name, "`!")
  }

  # Precompute unique (sample_feature, measured_specie, tl) combos across the entire dataset!
  joined_tl <- as.numeric(joined[[tl_info$tl_col]])
  joined_sf <- as.character(joined[[sample_feature]])
  joined_ft <- as.character(joined$measured_specie)

  combo_key <- paste(joined_sf, joined_ft, joined_tl, sep = "\r")
  u_combos <- unique(combo_key)
  u_first <- match(u_combos, combo_key)
  u_sf <- joined_sf[u_first]
  u_ft <- joined_ft[u_first]
  u_tl <- joined_tl[u_first]
  obs_to_u <- match(combo_key, u_combos)
  M_u <- length(u_combos)

  # Evaluate ODE forward predictions across all F_tot rows for the M_u unique condition/time tuples
  par_mat <- as.matrix(dyn_tbl[, param_names, drop = FALSE])
  npars <- ncol(par_mat)

  pred_fn_u <- matrix(NA_real_, nrow = F_tot, ncol = M_u)
  pred_logit_u <- matrix(NA_real_, nrow = F_tot, ncol = M_u)
  pred_ss_u <- matrix(NA_real_, nrow = F_tot, ncol = M_u)

  # Box bounds check for active inequality constraint detection
  lb <- dyn_meta[["lower_bounds"]]
  ub <- dyn_meta[["upper_bounds"]]
  if(is.null(lb) || is.null(ub) || length(lb) != npars || length(ub) != npars){
    bo <- if(max(graph) > 1L) seq(from = -0.01, to = 0.01, length.out = max(graph)) else 0
    lb_full <- -10 + bo
    ub_full <- 10 + bo
    if(use_cov && !from_fractions){
      lb_full[ksyn_index] <- 0
      ub_full[ksyn_index] <- 15
    }
    if(!use_cov){
      lb <- lb_full[-ksyn_index]
      ub <- ub_full[-ksyn_index]
    }else{
      lb <- lb_full
      ub <- ub_full
    }
  }

  # Which parameters actually shape fraction trajectories (exclude single-source ksyn)
  if(use_cov && length(ksyn_index) == 1L){
    fn_par_mask <- seq_len(npars) != ksyn_index
  }else{
    fn_par_mask <- rep(TRUE, npars)
  }
  n_params_fn <- sum(fn_par_mask)

  n_at_bounds_fn <- integer(F_tot)
  n_at_bounds_joint <- integer(F_tot)

  for(g in seq_len(F_tot)){
    p_g <- par_mat[g, ]
    at_b <- (p_g <= lb + 1e-3) | (p_g >= ub - 1e-3)
    at_b[is.na(at_b)] <- FALSE
    n_at_bounds_joint[g] <- sum(at_b)
    n_at_bounds_fn[g] <- sum(at_b[fn_par_mask])

    fwd <- evaluate_ode_forward(
      par_vec = p_g,
      graph = graph,
      formula_list = formula_list,
      tls = u_tl,
      sample_features = u_sf,
      feature_types = u_ft,
      use_coverage = use_cov
    )
    pred_fn_u[g, ] <- fwd$pred_fn
    pred_logit_u[g, ] <- fwd$pred_logit_fn
    pred_ss_u[g, ] <- fwd$pred_ss
  }

  lin_idx <- joined$.ez_dyn_id + (obs_to_u - 1L) * F_tot
  pred_fn_obs <- pred_fn_u[lin_idx]
  pred_logit_obs <- pred_logit_u[lin_idx]
  pred_ss_obs <- pred_ss_u[lin_idx]

  # Drop rows whose (sample_feature, measured_specie) has no formula in modeled_to_measured
  valid_obs <- is.finite(pred_fn_obs)
  if(!any(valid_obs)){
    stop("Forward ODE evaluation produced no valid predictions for `", source_table, "`!")
  }
  if(!all(valid_obs)){
    joined <- joined[valid_obs, , drop = FALSE]
    pred_fn_obs <- pred_fn_obs[valid_obs]
    pred_logit_obs <- pred_logit_obs[valid_obs]
    pred_ss_obs <- pred_ss_obs[valid_obs]
    combo_key <- combo_key[valid_obs]
  }

  obs_fn <- joined[[fcol_info$fraction_col]]
  obs_logit <- joined[[fcol_info$logit_col]]
  se_logit <- if(!is.null(fcol_info$se_col)) joined[[fcol_info$se_col]] else NULL
  ci_lo <- if(!is.null(fcol_info$ci_lo_col)) joined[[fcol_info$ci_lo_col]] else NULL
  ci_hi <- if(!is.null(fcol_info$ci_hi_col)) joined[[fcol_info$ci_hi_col]] else NULL

  obs_metrics <- compute_gof_obs_metrics(
    obs_fn = obs_fn,
    obs_logit_fn = obs_logit,
    se_logit_fn = se_logit,
    ci_lo_fn = ci_lo,
    ci_hi_fn = ci_hi,
    n_reads = joined$n,
    pred_fn = pred_fn_obs,
    pred_logit_fn = pred_logit_obs
  )

  # Restore original types on dyn_keys_df from dyn_tbl for clean output
  dyn_keys_orig <- dplyr::as_tibble(dyn_tbl[, dyn_group_cols, drop = FALSE])

  summary_tbl <- summarize_feature_gof(
    feat_keys_df = dyn_keys_orig,
    feat_id = joined$.ez_dyn_id,
    cell_key_str = combo_key,
    obs_fn = obs_fn,
    pred_fn = pred_fn_obs,
    pred_logit_fn = pred_logit_obs,
    n_reads = joined$n,
    obs_metrics = obs_metrics,
    n_params_vec = rep(n_params_fn, F_tot),
    n_at_bounds_vec = n_at_bounds_fn,
    dispersion_mode = dispersion_mode,
    logit_sd_floor = logit_sd_floor
  )

  # Optional companion steady-state coverage GoF when use_coverage == TRUE
  cov_cols_res <- NULL
  if(do_cov_gof){
    obs_cov_lin <- pmax(as.numeric(joined$.ez_obs_cov), 1e-8)
    obs_l10_cov <- log10(obs_cov_lin)
    pred_l10_cov <- log10(pmax(pred_ss_obs, 1e-8))
    resid_l10_cov <- obs_l10_cov - pred_l10_cov
    # Match dynamics_likelihood coverage SD (ODEs.R line 1379)
    sd_l10_cov <- log(10) * sqrt(1 / pmax(obs_cov_lin, 1) + 1 / 100)
    z_cov <- resid_l10_cov / sd_l10_cov

    # Cell-level coverage lack-of-fit (grouping across replicates of same sf + measured_specie + tl)
    cov_cell_str <- paste(joined$.ez_dyn_id, combo_key, sep = "\r")
    u_cc <- unique(cov_cell_str)
    cc_id <- match(cov_cell_str, u_cc)
    first_cc <- match(u_cc, cov_cell_str)
    feat_cc <- joined$.ez_dyn_id[first_cc]
    pred_cc <- pred_l10_cov[first_cc]

    w_cov <- 1 / (sd_l10_cov^2)
    W_cc <- as.numeric(rowsum(w_cov, cc_id, reorder = FALSE))
    bar_cc <- as.numeric(rowsum(w_cov * obs_l10_cov, cc_id, reorder = FALSE)) / pmax(W_cc, 1e-15)
    ss_lof_cc <- W_cc * (bar_cc - pred_cc)^2

    chisq_cov_lof <- numeric(F_tot)
    rs_cl <- rowsum(ss_lof_cc, feat_cc, reorder = FALSE)
    chisq_cov_lof[as.integer(rownames(rs_cl))] <- as.numeric(rs_cl)

    rs_wc <- rowsum(w_cov, joined$.ez_dyn_id, reorder = FALSE)
    rs_wc2 <- rowsum(w_cov * (resid_l10_cov^2), joined$.ez_dyn_id, reorder = FALSE)
    wrmse_cov <- rep(NA_real_, F_tot)
    f_c_pres <- as.integer(rownames(rs_wc))
    wrmse_cov[f_c_pres] <- sqrt(as.numeric(rs_wc2) / pmax(as.numeric(rs_wc), 1e-15))

    n_cov_groups <- tabulate(feat_cc, nbins = F_tot)
    p_joint_free <- pmax(npars - n_at_bounds_joint, 0L)
    df_lof_joint <- pmax((summary_tbl$n_groups + n_cov_groups) - p_joint_free, 0L)
    chisq_lof_joint <- summary_tbl$chisq_lof + (chisq_cov_lof / summary_tbl$dispersion_factor)
    j_bound_viol <- (df_lof_joint == 0L) & (n_at_bounds_joint > 0L) & (chisq_lof_joint > 1e-4)
    if(any(j_bound_viol)){
      df_lof_joint[j_bound_viol] <- 1L
    }
    pval_lof_joint <- ifelse(df_lof_joint > 0L,
                             stats::pchisq(chisq_lof_joint, df = df_lof_joint, lower.tail = FALSE),
                             NA_real_)
    padj_lof_joint <- rep(NA_real_, F_tot)
    ok_j <- !is.na(pval_lof_joint)
    if(any(ok_j)){
      padj_lof_joint[ok_j] <- stats::p.adjust(pval_lof_joint[ok_j], method = "BH")
    }

    summary_tbl$wrmse_log10_cov <- wrmse_cov
    summary_tbl$chisq_lof_cov <- chisq_cov_lof
    summary_tbl$chisq_lof_joint <- chisq_lof_joint
    summary_tbl$df_lof_joint <- as.integer(df_lof_joint)
    summary_tbl$pval_lof_joint <- pval_lof_joint
    summary_tbl$padj_lof_joint <- padj_lof_joint

    cov_cols_res <- dplyr::tibble(
      obs_log10_cov = obs_l10_cov,
      pred_log10_cov = pred_l10_cov,
      resid_log10_cov = resid_l10_cov,
      std_resid_cov = z_cov
    )
  }

  # Assemble tidy residuals table
  id_cols_keep <- unique(c(feat_cols, grouping_features))
  id_cols_keep <- id_cols_keep[id_cols_keep %in% colnames(joined)]
  sf_df <- if(!is_imputed_sf) dplyr::tibble(!!sample_feature := as.character(joined[[sample_feature]])) else NULL
  df_cond <- if(length(design_factors) > 0L) joined[, design_factors, drop = FALSE] else NULL

  residuals_tbl <- dplyr::bind_cols(
    dplyr::tibble(sample = as.character(joined$sample)),
    joined[, id_cols_keep, drop = FALSE],
    dplyr::tibble(measured_specie = as.character(joined$measured_specie)),
    sf_df,
    dplyr::tibble(tl = as.numeric(joined[[tl_info$tl_col]])),
    df_cond,
    dplyr::tibble(
      obs_fn = as.numeric(obs_fn),
      pred_fn = as.numeric(pred_fn_obs),
      resid_fn = as.numeric(obs_fn) - as.numeric(pred_fn_obs),
      obs_logit_fn = as.numeric(obs_logit),
      pred_logit_fn = as.numeric(pred_logit_obs),
      resid_logit_fn = as.numeric(obs_logit) - as.numeric(pred_logit_obs),
      std_resid = obs_metrics$std_resid,
      ci_lower_fn = obs_metrics$ci_lo_eff,
      ci_upper_fn = obs_metrics$ci_hi_eff,
      n = as.numeric(joined$n)
    ),
    cov_cols_res
  )

  list(
    summary_tbl = summary_tbl,
    residuals_tbl = residuals_tbl,
    metadata = list(
      source_type = "dynamics",
      source_table = source_table,
      fractions_table = fractions_name,
      fraction_col = fcol_info$fraction_col,
      tl_info = tl_info,
      features = dyn_group_cols,
      raw_features = feat_cols,
      sub_features = sub_features,
      grouping_features = grouping_features,
      sample_feature = sample_feature,
      modeled_to_measured = formula_list,
      graph = graph,
      parameter = dyn_meta[["parameter"]],
      parameter_names = param_names,
      use_coverage = use_cov,
      dynamics_design_factors = design_factors,
      scale_factors = scale_df,
      feature_lengths = dyn_meta[["feature_lengths"]]
    )
  )
}


# Evaluate smooth theoretical trajectories on a continuous time grid for EZFitPlot()
predict_gof_trajectories <- function(obj, gof_name, selected_keys_df, n_grid = 150L){
  gof_meta <- obj[["metadata"]][["gof"]][[gof_name]]
  res_tbl <- obj[["gof_residuals"]][[gof_name]]
  src_type <- gof_meta[["source_type"]]
  src_table <- gof_meta[["source_table"]]
  feat_cols <- gof_meta[["features"]]

  sub_res <- dplyr::inner_join(res_tbl, selected_keys_df, by = feat_cols)
  max_tl <- max(sub_res$tl, na.rm = TRUE)
  if(!is.finite(max_tl) || max_tl <= 0) max_tl <- 1

  if(src_type == "averages"){
    avg_tbl <- obj[["averages"]][[src_table]]
    avg_meta <- obj[["metadata"]][["averages"]][[src_table]]
    fit_params <- avg_meta[["fit_params"]]
    param_nm <- avg_meta[["parameter"]]
    k_strat <- gof_meta[["kstrat"]]
    tl_info <- gof_meta[["tl_info"]]
    is_direct <- grepl("^logit_fraction_", param_nm)

    metadf <- dplyr::as_tibble(obj$metadf)
    if(tl_info$mode == "pulse-chase"){
      metadf_lab <- metadf[as.numeric(metadf[[tl_info$tpulse_col]]) > 0, , drop = FALSE]
    }else{
      metadf_lab <- metadf[as.numeric(metadf[[tl_info$tl_col]]) > 0, , drop = FALSE]
    }
    X <- build_averages_design_matrix(metadf_lab, avg_meta)

    cond_cols <- setdiff(colnames(metadf_lab), c("sample", "tl", tl_info$tl_col,
                                                 tl_info$tpulse_col %||% "", tl_info$tchase_col %||% ""))
    # Unique condition profiles (unique rows of X paired with human-readable condition label)
    if(length(cond_cols) > 0L){
      cond_lbl_all <- apply(metadf_lab[, cond_cols, drop = FALSE], 1L, function(r) paste(r, collapse = ":"))
    }else{
      cond_lbl_all <- rep("All samples", nrow(metadf_lab))
    }
    x_sig <- apply(round(X, 8), 1L, paste, collapse = "|")
    u_x <- !duplicated(paste(x_sig, cond_lbl_all, sep = "\r"))
    X_u <- X[u_x, , drop = FALSE]
    cond_u <- cond_lbl_all[u_x]
    tl_u <- as.numeric(metadf_lab[[tl_info$tl_col]])[u_x]

    avg_sub <- dplyr::inner_join(avg_tbl, selected_keys_df, by = feat_cols)
    B_sub <- as.matrix(avg_sub[, paste0("mean_", fit_params), drop = FALSE])
    B_sub[!is.finite(B_sub)] <- 0
    L_u <- tcrossprod(B_sub, X_u) # n_sel x n_cond

    curve_rows <- list()
    if(identical(k_strat, "pulse-chase")){
      tp_val <- stats::median(as.numeric(metadf_lab[[tl_info$tpulse_col]]), na.rm = TRUE)
      tc_max <- max(as.numeric(metadf_lab[[tl_info$tchase_col]]), na.rm = TRUE)
      t_pulse_grid <- seq(0, tp_val, length.out = max(20L, round(n_grid * 0.35)))
      t_chase_grid <- seq(tp_val, tp_val + tc_max * 1.05, length.out = max(30L, round(n_grid * 0.65)))
      t_grid <- unique(c(t_pulse_grid, t_chase_grid))

      for(i in seq_len(nrow(avg_sub))){
        fk_i <- avg_sub[i, feat_cols, drop = FALSE]
        for(c_idx in seq_len(ncol(L_u))){
          k_val <- exp(L_u[i, c_idx])
          fn_grid <- ifelse(
            t_grid <= tp_val,
            -expm1(-k_val * t_grid),
            (-expm1(-k_val * tp_val)) * exp(-k_val * (t_grid - tp_val))
          )
          curve_rows[[length(curve_rows) + 1L]] <- dplyr::bind_cols(
            fk_i[rep(1L, length(t_grid)), , drop = FALSE],
            dplyr::tibble(tl = t_grid, pred_fn = fn_grid, curve_group = cond_u[c_idx])
          )
        }
      }
    }else if(is_direct && any(duplicated(cond_u))){
      # Direct logit_fraction model where formula_mean depends on tl: connect (0, 0) -> (tl, pred_fn) per condition
      u_cond_lvls <- unique(cond_u)
      for(i in seq_len(nrow(avg_sub))){
        fk_i <- avg_sub[i, feat_cols, drop = FALSE]
        for(clv in u_cond_lvls){
          c_idxs <- which(cond_u == clv)
          ord_t <- order(tl_u[c_idxs])
          c_idxs <- c_idxs[ord_t]
          t_pts <- c(0, tl_u[c_idxs])
          fn_pts <- c(0, inv_logit(L_u[i, c_idxs]))
          curve_rows[[length(curve_rows) + 1L]] <- dplyr::bind_cols(
            fk_i[rep(1L, length(t_pts)), , drop = FALSE],
            dplyr::tibble(tl = t_pts, pred_fn = fn_pts, curve_group = clv)
          )
        }
      }
    }else{
      t_grid <- seq(0, max_tl * 1.06, length.out = n_grid)
      for(i in seq_len(nrow(avg_sub))){
        fk_i <- avg_sub[i, feat_cols, drop = FALSE]
        for(c_idx in seq_len(ncol(L_u))){
          if(is_direct){
            fn_grid <- rep(inv_logit(L_u[i, c_idx]), length(t_grid))
          }else if(identical(k_strat, "shortfeed")){
            fn_grid <- pmin(exp(L_u[i, c_idx]) * t_grid, 1)
          }else{
            fn_grid <- -expm1(-exp(L_u[i, c_idx]) * t_grid)
          }
          curve_rows[[length(curve_rows) + 1L]] <- dplyr::bind_cols(
            fk_i[rep(1L, length(t_grid)), , drop = FALSE],
            dplyr::tibble(tl = t_grid, pred_fn = fn_grid, curve_group = cond_u[c_idx])
          )
        }
      }
    }
    return(dplyr::bind_rows(curve_rows))
  }

  # Dynamics trajectories
  dyn_tbl <- obj[["dynamics"]][[src_table]]
  graph <- gof_meta[["graph"]]
  formula_list <- gof_meta[["modeled_to_measured"]]
  param_names <- gof_meta[["parameter_names"]]
  use_cov <- gof_meta[["use_coverage"]]
  sf_col <- gof_meta[["sample_feature"]]
  is_imputed_sf <- identical(sf_col, "imputed_sample_feature")
  design_factors <- gof_meta[["dynamics_design_factors"]]

  dyn_sub <- dplyr::inner_join(dyn_tbl, selected_keys_df, by = feat_cols)
  t_grid <- seq(0, max_tl * 1.06, length.out = n_grid)

  # All valid (sample_feature, measured_specie) pairs in formula_list
  sf_names <- names(formula_list)
  sf_ft_pairs <- list()
  for(sf in sf_names){
    for(ft in names(formula_list[[sf]])){
      sf_ft_pairs[[length(sf_ft_pairs) + 1L]] <- list(sf = sf, ft = ft)
    }
  }

  eval_sf <- rep(vapply(sf_ft_pairs, `[[`, character(1), "sf"), each = length(t_grid))
  eval_ft <- rep(vapply(sf_ft_pairs, `[[`, character(1), "ft"), each = length(t_grid))
  eval_tl <- rep(t_grid, times = length(sf_ft_pairs))

  curve_rows <- vector("list", nrow(dyn_sub))
  for(i in seq_len(nrow(dyn_sub))){
    p_i <- as.numeric(dyn_sub[i, param_names])
    fwd_i <- evaluate_ode_forward(
      par_vec = p_i,
      graph = graph,
      formula_list = formula_list,
      tls = eval_tl,
      sample_features = eval_sf,
      feature_types = eval_ft,
      use_coverage = use_cov
    )
    fk_i <- dyn_sub[i, feat_cols, drop = FALSE]
    lbl_spec <- if(!is_imputed_sf && length(sf_names) > 1L){
      paste0(eval_sf, ":", eval_ft)
    }else{
      eval_ft
    }
    if(length(design_factors) > 0L){
      cond_str <- paste(as.character(dyn_sub[i, design_factors]), collapse = ":")
    }else{
      cond_str <- "All"
    }
    df_c <- dplyr::bind_cols(
      fk_i[rep(1L, length(eval_tl)), , drop = FALSE],
      dplyr::tibble(
        tl = eval_tl,
        pred_fn = fwd_i$pred_fn,
        measured_specie = eval_ft,
        curve_group = lbl_spec,
        condition_group = cond_str
      )
    )
    if(!is_imputed_sf){
      df_c[[sf_col]] <- eval_sf
    }
    curve_rows[[i]] <- df_c[is.finite(df_c$pred_fn), , drop = FALSE]
  }

  dplyr::bind_rows(curve_rows)
}


################################################################################
# MIXTURE MODEL GOODNESS OF FIT ENGINE (type = "fractions")
################################################################################

# Fast vectorized pure-data MLE for 2-component mixtures (skips profile CI loops)
fast_2comp_mle_vec <- function(grp_id, grp_end, n, delta, N_grp = NULL,
                               eta_min = -25, eta_max = 25){
  ngroups <- length(grp_end)
  if(ngroups == 0L){
    return(list(
      pi_mle = numeric(0),
      eta_mle = numeric(0),
      at_0 = logical(0),
      at_1 = logical(0)
    ))
  }
  if(is.null(N_grp)) N_grp <- grp_sum_contig(n, grp_end)

  delta_clean <- delta
  delta_clean[is.nan(delta_clean)] <- 0
  delta_clean <- pmin(pmax(delta_clean, -300), 300)

  log_n <- log(pmax(n, 1e-300))
  log_N <- log(pmax(N_grp, 1e-300))

  lse0 <- grp_logsumexp(log_n + delta_clean, grp_id, grp_end)
  lse1 <- grp_logsumexp(log_n - delta_clean, grp_id, grp_end)

  at_0 <- lse0 <= log_N + 1e-13
  at_1 <- !at_0 & (lse1 <= log_N + 1e-13)

  eta_mle <- numeric(ngroups)
  eta_mle[at_0] <- -Inf
  eta_mle[at_1] <- Inf

  eval_grad <- function(eta_v, sg, se, sn, sdel, sN){
    theta <- stats::plogis(eta_v)
    w <- stats::plogis(eta_v[sg] + sdel)
    nw <- sn * w
    g1 <- grp_sum_contig(nw, se) - sN * theta
    info <- sN * theta * (1 - theta) - grp_sum_contig(nw * (1 - w), se)
    list(g1 = g1, info = info)
  }

  interior <- which(!(at_0 | at_1))
  if(length(interior) > 0L){
    active <- interior
    act_map <- integer(ngroups)
    act_map[active] <- seq_along(active)
    rk0 <- act_map[grp_id] > 0L
    sg <- act_map[grp_id[rk0]]
    gsz <- diff(c(0L, grp_end))[active]
    se <- cumsum(gsz)
    sn <- n[rk0]
    sdel <- delta_clean[rk0]
    sN <- N_grp[active]

    lo_a <- rep(eta_min, length(active))
    hi_a <- rep(eta_max, length(active))
    eta_a <- rep(0, length(active))

    g_lo <- eval_grad(lo_a, sg, se, sn, sdel, sN)$g1
    g_hi <- eval_grad(hi_a, sg, se, sn, sdel, sN)$g1
    pin_lo <- g_lo <= 0
    pin_hi <- !pin_lo & (g_hi >= 0)
    if(any(pin_lo | pin_hi)){
      if(any(pin_lo)) eta_mle[active[pin_lo]] <- eta_min
      if(any(pin_hi)) eta_mle[active[pin_hi]] <- eta_max
      keep <- which(!(pin_lo | pin_hi))
      if(length(keep) > 0L){
        remap <- integer(length(active))
        remap[keep] <- seq_along(keep)
        rk <- remap[sg] > 0L
        sg <- remap[sg[rk]]
        gsz <- gsz[keep]
        se <- cumsum(gsz)
        sn <- sn[rk]; sdel <- sdel[rk]; sN <- sN[keep]
        active <- active[keep]
        lo_a <- lo_a[keep]; hi_a <- hi_a[keep]; eta_a <- eta_a[keep]
      }else{
        active <- integer(0)
      }
    }

    if(length(active) > 0L){
      for(iter in seq_len(25L)){
        ev <- eval_grad(eta_a, sg, se, sn, sdel, sN)
        pos <- ev$g1 > 0
        lo_a[pos] <- pmax(lo_a[pos], eta_a[pos])
        hi_a[!pos] <- pmin(hi_a[!pos], eta_a[!pos])

        step <- ifelse(ev$info > 1e-12, ev$g1 / ev$info, NA_real_)
        cand <- eta_a + step
        bad <- is.na(cand) | cand <= lo_a | cand >= hi_a
        cand[bad] <- 0.5 * (lo_a[bad] + hi_a[bad])

        d_eta <- abs(cand - eta_a)
        eta_a <- cand
        conv <- (d_eta < 1e-9 & abs(ev$g1) < 1e-7) | ((hi_a - lo_a) < 1e-9)
        if(all(conv)){
          eta_mle[active] <- eta_a
          break
        }
        if(iter %% 3L == 0L && sum(conv) > 0.25 * length(conv)){
          done <- which(conv)
          eta_mle[active[done]] <- eta_a[done]
          rem <- which(!conv)
          active <- active[rem]
          eta_a <- eta_a[rem]; lo_a <- lo_a[rem]; hi_a <- hi_a[rem]
          sN <- sN[rem]; gsz <- gsz[rem]
          remap <- integer(length(conv))
          remap[rem] <- seq_along(rem)
          rk <- remap[sg] > 0L
          sg <- remap[sg[rk]]; se <- cumsum(gsz)
          sn <- sn[rk]; sdel <- sdel[rk]
        }
        if(iter == 25L) eta_mle[active] <- eta_a
      }
    }
  }

  pi_mle <- stats::plogis(eta_mle)
  pi_mle[at_0] <- 0
  pi_mle[at_1] <- 1

  list(
    pi_mle = pi_mle,
    eta_mle = eta_mle,
    at_0 = at_0,
    at_1 = at_1
  )
}


# Vectorized adaptive bin pooling guaranteeing E_pool >= min_bin_count in active bins
pool_gof_bins_vec <- function(O_mat, E_mat, min_bin_count = 5){
  nr <- nrow(E_mat)
  nc <- ncol(E_mat)
  O_pool <- O_mat
  E_pool <- E_mat
  active <- matrix(TRUE, nrow = nr, ncol = nc)

  if(nr == 0L || nc <= 1L){
    return(list(
      O_pool = O_pool,
      E_pool = E_pool,
      active = active,
      C_act = rowSums(active)
    ))
  }

  # Sweep 1: pool sparse upper tails downward (c -> c - 1)
  for(c_idx in seq.int(nc, 2L, by = -1L)){
    merge_down <- active[, c_idx] & (E_pool[, c_idx] < min_bin_count)
    if(any(merge_down)){
      E_pool[merge_down, c_idx - 1L] <- E_pool[merge_down, c_idx - 1L] + E_pool[merge_down, c_idx]
      O_pool[merge_down, c_idx - 1L] <- O_pool[merge_down, c_idx - 1L] + O_pool[merge_down, c_idx]
      E_pool[merge_down, c_idx] <- 0
      O_pool[merge_down, c_idx] <- 0
      active[merge_down, c_idx] <- FALSE
    }
  }

  # Sweep 2: pool sparse lower bins upward into the next surviving active bin > c_idx
  for(c_idx in seq_len(nc - 1L)){
    need_up <- which(active[, c_idx] & (E_pool[, c_idx] < min_bin_count))
    if(length(need_up) > 0L){
      act_right <- active[need_up, (c_idx + 1L):nc, drop = FALSE]
      has_right <- rowSums(act_right) > 0L
      if(any(has_right)){
        rows_u <- need_up[has_right]
        ar_u <- act_right[has_right, , drop = FALSE]
        # First TRUE column in (c_idx + 1):nc for each row
        next_col <- c_idx + max.col(ar_u, ties.method = "first")
        lin_src <- cbind(rows_u, rep.int(c_idx, length(rows_u)))
        lin_dst <- cbind(rows_u, next_col)
        E_pool[lin_dst] <- E_pool[lin_dst] + E_pool[lin_src]
        O_pool[lin_dst] <- O_pool[lin_dst] + O_pool[lin_src]
        E_pool[lin_src] <- 0
        O_pool[lin_src] <- 0
        active[lin_src] <- FALSE
      }
    }
  }

  list(
    O_pool = O_pool,
    E_pool = E_pool,
    active = active,
    C_act = rowSums(active)
  )
}


# Extract exact conditional Binomial read counts (sample, features, mutcols, basecols, n)
# from either in-memory EZbakRData$cB or on-disk EZbakRArrowData
extract_cb_for_gof <- function(obj, samples_keep, samples_nolabel,
                               feat_cols, mutcols, basecols,
                               fxn_meta, features_keep_df = NULL){
  sample <- n <- NULL
  `.` <- list

  remove_features <- fxn_meta[["remove_features"]] %||% c("NA", "__no_feature")
  filter_cols <- fxn_meta[["filter_cols"]] %||% feat_cols
  if(length(filter_cols) == 1L && identical(filter_cols, "all")){
    filter_cols <- feat_cols
  }
  filter_cols <- intersect(filter_cols, feat_cols)
  remove_mf <- isTRUE(fxn_meta[["remove_multi_features"]])
  split_mf <- isTRUE(fxn_meta[["split_multi_features"]])
  mf_cols <- fxn_meta[["multi_feature_cols"]]
  mf_sep <- fxn_meta[["multi_feature_sep"]] %||% "+"
  mf_cols_rm <- if(is.null(mf_cols) || (length(mf_cols) == 1L && identical(mf_cols, "all"))){
    feat_cols
  }else{
    intersect(mf_cols, feat_cols)
  }

  cols_grp <- unique(c("sample", feat_cols, mutcols, basecols))
  pold_nolabel <- list()
  nl_keep <- intersect(samples_keep, samples_nolabel)

  keep_dt <- NULL
  if(!is.null(features_keep_df)){
    keep_dt <- data.table::as.data.table(dplyr::distinct(features_keep_df[, feat_cols, drop = FALSE]))
    for(fc in feat_cols) keep_dt[[fc]] <- as.character(keep_dt[[fc]])
    data.table::setkeyv(keep_dt, feat_cols)
  }

  filter_and_split_dt <- function(dt_in){
    dt_w <- data.table::as.data.table(dt_in)
    if(nrow(dt_w) == 0L) return(dt_w)
    if(length(filter_cols) > 0L && length(remove_features) > 0L){
      bad_row <- Reduce(`|`, lapply(dt_w[, filter_cols, with = FALSE], `%in%`, remove_features))
      dt_w <- dt_w[!bad_row]
    }
    if(nrow(dt_w) == 0L) return(dt_w)
    if(remove_mf && length(mf_cols_rm) > 0L){
      dt_w <- remove_multi_feature_rows(dt_w, mf_cols_rm, sep = mf_sep)
    }
    if(nrow(dt_w) == 0L) return(dt_w)
    if(split_mf && !is.null(mf_cols) && length(mf_cols) > 0L){
      dt_w <- data.table::as.data.table(split_features(dt_w, mf_cols, cols_grp, sep = mf_sep))
    }
    for(fc in feat_cols) dt_w[[fc]] <- as.character(dt_w[[fc]])
    if(!is.null(keep_dt)){
      dt_w <- dt_w[keep_dt, on = feat_cols, nomatch = NULL]
    }
    dt_w
  }

  if(methods::is(obj, "EZbakRArrowData")){
    ds <- get_arrow_ds(obj)
    sel_cols <- unique(c(cols_grp, "n"))
    parts <- vector("list", length(samples_keep))

    for(si in seq_along(samples_keep)){
      s_cur <- samples_keep[si]
      s_raw <- ds %>%
        dplyr::filter(sample == s_cur) %>%
        dplyr::select(dplyr::all_of(sel_cols)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(cols_grp))) %>%
        dplyr::summarise(n = sum(n), .groups = "drop") %>%
        dplyr::collect() %>%
        dplyr::mutate(
          n = as.numeric(n),
          dplyr::across(dplyr::all_of(unique(c(mutcols, basecols))), as.numeric)
        )

      if(s_cur %in% nl_keep && nrow(s_raw) > 0L){
        rates_nl <- vapply(seq_along(mutcols), function(j){
          num_j <- sum(as.numeric(s_raw[[mutcols[j]]]) * s_raw$n)
          den_j <- max(sum(as.numeric(s_raw[[basecols[j]]]) * s_raw$n), 1)
          pmin(pmax(num_j / den_j, 1e-6), 0.25)
        }, numeric(1))
        names(rates_nl) <- mutcols
        pold_nolabel[[s_cur]] <- rates_nl
      }

      parts[[si]] <- filter_and_split_dt(s_raw)
    }

    cb_dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  }else{
    cB_raw <- data.table::as.data.table(obj[["cB"]])
    cB_sub <- cB_raw[sample %in% samples_keep, .(n = as.numeric(sum(n))), by = cols_grp]

    for(s_cur in nl_keep){
      s_rows <- cB_sub[sample == s_cur]
      if(nrow(s_rows) > 0L){
        rates_nl <- vapply(seq_along(mutcols), function(j){
          num_j <- sum(as.numeric(s_rows[[mutcols[j]]]) * s_rows$n)
          den_j <- max(sum(as.numeric(s_rows[[basecols[j]]]) * s_rows$n), 1)
          pmin(pmax(num_j / den_j, 1e-6), 0.25)
        }, numeric(1))
        names(rates_nl) <- mutcols
        pold_nolabel[[s_cur]] <- rates_nl
      }
    }

    cb_dt <- filter_and_split_dt(cB_sub)
  }

  list(
    cb_dt = cb_dt,
    pold_nolabel = pold_nolabel
  )
}


# Evaluate mixture-model goodness of fit for EstimateFractions() output
compute_fractions_gof <- function(obj,
                                  source_table,
                                  min_reads = 10,
                                  min_samples = NULL,
                                  feature_sample_counts = NULL,
                                  min_bin_count = 5,
                                  sample_padj_cutoff = 0.05,
                                  sample_tvd_cutoff = 0.05,
                                  dispersion_mode = "empirical"){
  sample <- n <- nsamps <- tvd <- kl_div <- NULL
  failure_mode <- is_nolabel <- cnt <- sum_cnt <- NULL
  .raw_mode <- .g2_s <- .df_s <- NULL

  fxn_tbl <- as.data.frame(obj[["fractions"]][[source_table]])
  fxn_meta <- obj[["metadata"]][["fractions"]][[source_table]]
  feat_cols <- fxn_meta[["features"]]
  fraction_design <- as.data.frame(fxn_meta[["fraction_design"]])
  mutrate_design <- fraction_design[as.logical(fraction_design$present),
                                    setdiff(colnames(fraction_design), "present"),
                                    drop = FALSE]

  mutcols <- colnames(mutrate_design)
  basecols <- paste0("n", substr(mutcols, 1L, 1L))
  K <- nrow(mutrate_design)
  J <- ncol(mutrate_design)
  md_mat <- as.matrix(mutrate_design)

  # Identify strategy ("standard", "hierarchical", or "general")
  hier_cand_name <- paste0(paste(feat_cols, collapse = "_"), "_", mutcols[1L])
  if(J == 1L && K == 2L){
    strat_meta <- fxn_meta[["strategy"]]
    if(identical(strat_meta, "hierarchical") ||
       (is.null(strat_meta) && (hier_cand_name %in% names(obj[["mutation_rates"]])))){
      eff_strategy <- "hierarchical"
    }else{
      eff_strategy <- "standard"
    }
  }else{
    eff_strategy <- "general"
  }
  hier_table_name <- fxn_meta[["hier_mutrate_table"]] %||% hier_cand_name

  # Identify labeled vs -label samples from metadf
  metadf <- as.data.frame(obj[["metadf"]])
  metadf$sample <- as.character(metadf$sample)
  tl_cols_possible <- if(J == 1L){
    c("tl", "tpulse")
  }else{
    c(paste0("tl_", mutcols), paste0("tpulse_", mutcols), "tl", "tpulse")
  }
  tl_cols <- intersect(tl_cols_possible, colnames(metadf))
  if(length(tl_cols) == 0L){
    stop("Could not find label-time columns (`tl` / `tpulse`) in `obj$metadf`!")
  }

  tl_mat <- as.matrix(metadf[, tl_cols, drop = FALSE])
  storage.mode(tl_mat) <- "double"
  is_nl_meta <- apply(tl_mat == 0, 1L, all)
  samples_nolabel <- metadf$sample[is_nl_meta]
  samples_labeled <- metadf$sample[!is_nl_meta]

  tl_disp_vec <- if("tl" %in% colnames(metadf)){
    as.numeric(metadf$tl)
  }else if("tpulse" %in% colnames(metadf)){
    as.numeric(metadf$tpulse)
  }else{
    apply(tl_mat, 1L, max)
  }
  names(tl_disp_vec) <- metadf$sample

  cond_cols <- setdiff(colnames(metadf), c("sample", "tl"))

  # Coerce feature columns to character on our local copy (never touching obj$fractions)
  fxn_tbl$sample <- as.character(fxn_tbl$sample)
  for(fc in feat_cols) fxn_tbl[[fc]] <- as.character(fxn_tbl[[fc]])

  fxn_labeled <- fxn_tbl[fxn_tbl$sample %in% samples_labeled, , drop = FALSE]
  num_labeled_samps <- length(unique(fxn_labeled$sample))
  if(num_labeled_samps == 0L){
    stop("No labeled samples found in the selected `fractions` table!")
  }

  ### Step 1: Filter features strictly for AssessGoF output (mirrors AverageAndRegularize)
  cnt_labeled <- fxn_labeled %>%
    dplyr::filter(n > min_reads) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(feat_cols))) %>%
    dplyr::summarise(cnt = dplyr::n(), .groups = "drop")

  if(!is.null(feature_sample_counts)){
    if(!is.data.frame(feature_sample_counts) || !("nsamps" %in% colnames(feature_sample_counts))){
      stop("`feature_sample_counts` must be a data frame containing `<feature columns>` and `nsamps`!")
    }
    fsc_df <- as.data.frame(feature_sample_counts)
    filter_features <- setdiff(colnames(fsc_df), "nsamps")
    if(length(filter_features) == 0L || !all(filter_features %in% feat_cols)){
      stop("`feature_sample_counts` feature columns must be a subset of `features` (",
           paste(feat_cols, collapse = ", "), ")!")
    }
    for(fc in filter_features) fsc_df[[fc]] <- as.character(fsc_df[[fc]])

    keep_groups <- cnt_labeled %>%
      dplyr::inner_join(fsc_df, by = filter_features) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(filter_features, "nsamps")))) %>%
      dplyr::summarise(sum_cnt = sum(cnt), .groups = "drop") %>%
      dplyr::filter(sum_cnt >= nsamps) %>%
      dplyr::select(dplyr::all_of(filter_features))

    distinct_feats <- dplyr::distinct(fxn_labeled[, feat_cols, drop = FALSE])
    features_to_keep <- dplyr::inner_join(distinct_feats, keep_groups, by = filter_features)
  }else{
    req_samps <- if(is.null(min_samples)){
      num_labeled_samps
    }else{
      as.integer(min_samples)
    }
    if(length(req_samps) != 1L || is.na(req_samps) || req_samps < 1L){
      stop("`min_samples` must be a positive integer!")
    }
    features_to_keep <- cnt_labeled %>%
      dplyr::filter(cnt >= req_samps) %>%
      dplyr::select(dplyr::all_of(feat_cols))
  }

  if(nrow(features_to_keep) == 0L){
    stop("No features made it past read-count filtering (`min_reads = ", min_reads, "`) in `AssessGoF()`!")
  }

  fxn_eval <- fxn_tbl %>%
    dplyr::inner_join(features_to_keep, by = feat_cols) %>%
    dplyr::filter(n > min_reads)

  samples_eval <- unique(fxn_eval$sample)

  ### Step 2: Extract read-level conditional Binomial counts for surviving features
  ext <- extract_cb_for_gof(
    obj = obj,
    samples_keep = samples_eval,
    samples_nolabel = samples_nolabel,
    feat_cols = feat_cols,
    mutcols = mutcols,
    basecols = basecols,
    fxn_meta = fxn_meta,
    features_keep_df = features_to_keep
  )
  cb_dt <- ext$cb_dt
  pold_nolabel <- ext$pold_nolabel

  if(nrow(cb_dt) == 0L){
    stop("No matching read-level rows found in `cB` for the selected `fractions` table!")
  }

  # Helper for 1D mechanistic failure-mode tagging
  classify_1d_failure <- function(P_obs, P_exp){
    nr <- nrow(P_obs)
    if(nr == 0L) return(character(0))
    diff_p <- P_obs - P_exp
    exc_inter <- pmax(diff_p[, 2L], 0) + pmax(diff_p[, 3L], 0)
    exc_tail <- rowSums(pmax(diff_p[, 5:7, drop = FALSE], 0))
    m_mid <- c(0, 1, 2, 3, 4, 5, 6.5)
    mu_obs <- as.vector(P_obs %*% m_mid)
    mu_exp <- as.vector(P_exp %*% m_mid)
    var_obs <- pmax(as.vector(P_obs %*% (m_mid * m_mid)) - mu_obs * mu_obs, 0)
    var_exp <- pmax(as.vector(P_exp %*% (m_mid * m_mid)) - mu_exp * mu_exp, 0)

    out_mode <- rep("intermediate_spike", nr)
    is_overdisp <- (diff_p[, 1L] > 0.02 & exc_tail > 0.01) | (var_obs > 1.25 * var_exp & diff_p[, 1L] > 0.005)
    out_mode[is_overdisp] <- "overdispersed_muts"
    is_tail <- (exc_tail >= exc_inter) & (exc_tail > 0.02)
    out_mode[is_tail] <- "high_mut_tail"
    out_mode
  }

  ### Step 3: Evaluate sample-by-sample exact conditional Binomial mixture GoF
  sample_res_list <- vector("list", length(samples_eval))
  disp_records <- vector("list", length(samples_eval))

  # Precompute population column names for general mixture
  if(eff_strategy == "general"){
    pop_labels <- vapply(seq_len(K), function(k_i){
      pv <- ifelse(as.logical(md_mat[k_i, ]), "high", "low")
      paste(paste0(pv, mutcols), collapse = "_")
    }, character(1))
    frac_cols_gen <- paste0("fraction_", pop_labels)
    low_rows_gen <- which(rowSums(md_mat) == 0)
    has_all_low <- length(low_rows_gen) > 0L
    low_idx_gen <- if(has_all_low) low_rows_gen[1L] else NA_integer_
  }else{
    frac_col_1d <- paste0("fraction_high", mutcols[1L])
  }

  for(si in seq_along(samples_eval)){
    s_cur <- samples_eval[si]
    s_is_nl <- s_cur %in% samples_nolabel

    s_dt <- cb_dt[sample == s_cur]
    if(nrow(s_dt) == 0L) next

    fxn_s_dt <- data.table::as.data.table(fxn_eval[fxn_eval$sample == s_cur, , drop = FALSE])
    data.table::setkeyv(fxn_s_dt, feat_cols)
    s_dt <- s_dt[fxn_s_dt[, feat_cols, with = FALSE], on = feat_cols, nomatch = NULL]
    if(nrow(s_dt) == 0L) next

    data.table::setorderv(s_dt, feat_cols)
    grp_id <- data.table::rleidv(s_dt, cols = feat_cols)
    n_rows <- nrow(s_dt)
    grp_end <- if(n_rows > 1L){
      c(which(grp_id[-1L] != grp_id[-n_rows]), n_rows)
    }else{
      1L
    }
    ng <- length(grp_end)
    n_vec <- as.numeric(s_dt[["n"]])
    N_grp <- grp_sum_contig(n_vec, grp_end)

    grp_keys_dt <- s_dt[grp_end, feat_cols, with = FALSE]
    fxn_matched <- fxn_s_dt[grp_keys_dt, on = feat_cols]

    if(eff_strategy %in% c("standard", "hierarchical")){
      # 1D 7-bin exact conditional Binomial mixture (bins: 0, 1, 2, 3, 4, 5, 6+)
      C_bins <- 7L
      m_vec <- as.integer(round(as.numeric(s_dt[[mutcols[1L]]])))
      b_vec <- as.integer(round(as.numeric(s_dt[[basecols[1L]]])))
      bin_idx <- pmin(pmax(m_vec, 0L), 6L) + 1L

      O_mat <- matrix(0, nrow = ng, ncol = C_bins)
      for(c_i in seq_len(C_bins)){
        O_mat[, c_i] <- grp_sum_contig(ifelse(bin_idx == c_i, n_vec, 0), grp_end)
      }

      est_fn <- as.numeric(fxn_matched[[frac_col_1d]])
      max_b <- max(b_vec, 0L)
      b_grid <- 0:max_b

      if(s_is_nl){
        pold_s <- pold_nolabel[[s_cur]][1L] %||% 1e-4
        pold_s <- pmin(pmax(pold_s, 1e-7), 1 - 1e-7)
        P_old_lut <- matrix(0, nrow = max_b + 1L, ncol = C_bins)
        for(km in 0:5) P_old_lut[, km + 1L] <- stats::dbinom(km, b_grid, pold_s)
        P_old_lut[, 7L] <- stats::pbinom(5L, b_grid, pold_s, lower.tail = FALSE)

        E_mat <- matrix(0, nrow = ng, ncol = C_bins)
        b_idx1 <- b_vec + 1L
        for(c_i in seq_len(C_bins)){
          E_mat[, c_i] <- grp_sum_contig(n_vec * P_old_lut[b_idx1, c_i], grp_end)
        }
        mle_fn <- rep(0, ng)
        p_eff <- rep(0, ng)
        at_upper_bound <- rep(FALSE, ng)

      }else if(eff_strategy == "standard"){
        mr_tbl <- as.data.frame(obj[["mutation_rates"]][[mutcols[1L]]])
        mr_row <- mr_tbl[as.character(mr_tbl$sample) == s_cur, , drop = FALSE]
        pnew_s <- pmin(pmax(as.numeric(mr_row$pnew[1L]), 1e-7), 1 - 1e-7)
        pold_s <- pmin(pmax(as.numeric(mr_row$pold[1L]), 1e-7), 1 - 1e-7)

        ll_new <- stats::dbinom(m_vec, b_vec, pnew_s, log = TRUE)
        ll_old <- stats::dbinom(m_vec, b_vec, pold_s, log = TRUE)
        mle_res <- fast_2comp_mle_vec(grp_id, grp_end, n_vec, ll_new - ll_old, N_grp)
        mle_fn <- mle_res$pi_mle
        at_upper_bound <- mle_res$at_1 | (mle_fn >= 1 - 1e-6)

        P_new_lut <- matrix(0, nrow = max_b + 1L, ncol = C_bins)
        P_old_lut <- matrix(0, nrow = max_b + 1L, ncol = C_bins)
        for(km in 0:5){
          P_new_lut[, km + 1L] <- stats::dbinom(km, b_grid, pnew_s)
          P_old_lut[, km + 1L] <- stats::dbinom(km, b_grid, pold_s)
        }
        P_new_lut[, 7L] <- stats::pbinom(5L, b_grid, pnew_s, lower.tail = FALSE)
        P_old_lut[, 7L] <- stats::pbinom(5L, b_grid, pold_s, lower.tail = FALSE)

        E_new <- matrix(0, nrow = ng, ncol = C_bins)
        E_old <- matrix(0, nrow = ng, ncol = C_bins)
        b_idx1 <- b_vec + 1L
        for(c_i in seq_len(C_bins)){
          E_new[, c_i] <- grp_sum_contig(n_vec * P_new_lut[b_idx1, c_i], grp_end)
          E_old[, c_i] <- grp_sum_contig(n_vec * P_old_lut[b_idx1, c_i], grp_end)
        }
        E_mat <- mle_fn * E_new + (1 - mle_fn) * E_old
        p_eff <- rep(1, ng)

      }else{
        # Hierarchical feature-specific pnew 2-component model
        hmr <- data.table::as.data.table(obj[["mutation_rates"]][[hier_table_name]])
        hmr$sample <- as.character(hmr$sample)
        for(fc in feat_cols) hmr[[fc]] <- as.character(hmr[[fc]])
        hmr_s <- hmr[sample == s_cur]
        data.table::setkeyv(hmr_s, feat_cols)
        hmr_matched <- hmr_s[grp_keys_dt, on = feat_cols]

        pnew_g <- pmin(pmax(as.numeric(hmr_matched$pnew), 1e-7), 1 - 1e-7)
        pold_s <- pmin(pmax(as.numeric(hmr_matched$pold[1L]), 1e-7), 1 - 1e-7)
        pnew_r <- pnew_g[grp_id]

        ll_new <- stats::dbinom(m_vec, b_vec, pnew_r, log = TRUE)
        ll_old <- stats::dbinom(m_vec, b_vec, pold_s, log = TRUE)
        mle_res <- fast_2comp_mle_vec(grp_id, grp_end, n_vec, ll_new - ll_old, N_grp)
        mle_fn <- mle_res$pi_mle
        at_upper_bound <- mle_res$at_1 | (mle_fn >= 1 - 1e-6)

        P_old_lut <- matrix(0, nrow = max_b + 1L, ncol = C_bins)
        for(km in 0:5) P_old_lut[, km + 1L] <- stats::dbinom(km, b_grid, pold_s)
        P_old_lut[, 7L] <- stats::pbinom(5L, b_grid, pold_s, lower.tail = FALSE)

        E_new <- matrix(0, nrow = ng, ncol = C_bins)
        E_old <- matrix(0, nrow = ng, ncol = C_bins)
        b_idx1 <- b_vec + 1L
        for(km in 0:5){
          E_new[, km + 1L] <- grp_sum_contig(n_vec * stats::dbinom(km, b_vec, pnew_r), grp_end)
          E_old[, km + 1L] <- grp_sum_contig(n_vec * P_old_lut[b_idx1, km + 1L], grp_end)
        }
        E_new[, 7L] <- grp_sum_contig(n_vec * stats::pbinom(5L, b_vec, pnew_r, lower.tail = FALSE), grp_end)
        E_old[, 7L] <- grp_sum_contig(n_vec * P_old_lut[b_idx1, 7L], grp_end)

        E_mat <- mle_fn * E_new + (1 - mle_fn) * E_old

        # Effective parameter count p_eff = 1 + I_22 / (I_22 + tau_2) in [1, 2]
        hp_df <- fxn_meta[["hier_prior"]]
        prior_sd_s <- NA_real_
        if(!is.null(hp_df) && is.data.frame(hp_df) && "pnew_prior_sd" %in% colnames(hp_df)){
          idx_hp <- which(as.character(hp_df$sample) == s_cur)
          if(length(idx_hp) > 0L) prior_sd_s <- as.numeric(hp_df$pnew_prior_sd[idx_hp[1L]])
        }
        if(!is.finite(prior_sd_s) || prior_sd_s <= 0){
          lp_all <- stats::qlogis(pmin(pmax(as.numeric(hmr_s$pnew), 1e-5), 1 - 1e-5))
          prior_sd_s <- if(length(lp_all) >= 2L) pmin(pmax(stats::sd(lp_all), 0.05), 0.5) else 0.15
        }
        tau_2 <- 1 / (prior_sd_s * prior_sd_s)

        eta_clamped <- stats::qlogis(pmin(pmax(mle_fn, 1e-12), 1 - 1e-12))
        w_r <- stats::plogis(eta_clamped[grp_id] + pmin(pmax(ll_new - ll_old, -300), 300))
        w_r[mle_fn[grp_id] <= 0] <- 0
        w_r[mle_fn[grp_id] >= 1] <- 1
        d1_r <- m_vec - b_vec * pnew_r
        info_22 <- pmax(
          grp_sum_contig(n_vec * (w_r * b_vec * pnew_r * (1 - pnew_r) - w_r * (1 - w_r) * d1_r * d1_r), grp_end),
          0
        )
        p_eff <- 1 + (info_22 / (info_22 + tau_2))
      }

      raw_diag_mode <- if(s_is_nl) rep("nolabel_anomaly", ng) else classify_1d_failure(O_mat / N_grp, E_mat / N_grp)

    }else{
      # General K-component / multi-channel exact conditional Binomial mixture
      B1 <- if(J == 1L) 7L else if(J == 2L) 4L else if(J == 3L) 3L else 2L
      cap_m <- B1 - 1L
      C_bins <- as.integer(B1^J)

      m_list <- vector("list", J)
      b_list <- vector("list", J)
      cell_idx <- rep(1L, n_rows)
      stride <- 1L
      for(j in seq_len(J)){
        mj <- as.integer(round(as.numeric(s_dt[[mutcols[j]]])))
        bj <- as.integer(round(as.numeric(s_dt[[basecols[j]]])))
        m_list[[j]] <- mj
        b_list[[j]] <- bj
        cell_idx <- cell_idx + pmin(pmax(mj, 0L), cap_m) * stride
        stride <- stride * B1
      }

      O_mat <- matrix(0, nrow = ng, ncol = C_bins)
      for(c_i in seq_len(C_bins)){
        O_mat[, c_i] <- grp_sum_contig(ifelse(cell_idx == c_i, n_vec, 0), grp_end)
      }

      frac_est_mat <- as.matrix(fxn_matched[, frac_cols_gen, with = FALSE])
      storage.mode(frac_est_mat) <- "double"

      if(s_is_nl){
        pold_s <- pold_nolabel[[s_cur]]
        if(is.null(pold_s)){
          pold_s <- rep(1e-4, J)
          names(pold_s) <- mutcols
        }
        pold_s <- pmin(pmax(pold_s, 1e-7), 1 - 1e-7)
        pnew_s <- pold_s

        pi_mle <- matrix(0, nrow = ng, ncol = K)
        if(has_all_low){
          pi_mle[, low_idx_gen] <- 1
        }else{
          pi_mle[, 1L] <- 1
        }
        p_eff <- rep(0, ng)
        est_fn <- rep(0, ng)
        mle_fn <- rep(0, ng)
        at_upper_bound <- rep(FALSE, ng)
      }else{
        pnew_s <- numeric(J)
        pold_s <- numeric(J)
        for(j in seq_len(J)){
          nm_j <- mutcols[j]
          mr_j <- as.data.frame(obj[["mutation_rates"]][[nm_j]])
          mr_row <- mr_j[as.character(mr_j$sample) == s_cur, , drop = FALSE]
          if(nrow(mr_row) > 0L && is.finite(mr_row$pnew[1L]) && is.finite(mr_row$pold[1L])){
            pnew_s[j] <- pmin(pmax(as.numeric(mr_row$pnew[1L]), 1e-7), 1 - 1e-7)
            pold_s[j] <- pmin(pmax(as.numeric(mr_row$pold[1L]), 1e-7), 1 - 1e-7)
          }else{
            emp_r <- sum(m_list[[j]] * n_vec) / max(sum(b_list[[j]] * n_vec), 1)
            pold_s[j] <- pmin(pmax(emp_r, 1e-6), 0.25)
            pnew_s[j] <- pold_s[j]
          }
        }

        # Exact conditional Binomial log-likelihood matrix (n_rows x K)
        L_mat <- matrix(0, nrow = n_rows, ncol = K)
        for(j in seq_len(J)){
          ll_hi <- stats::dbinom(m_list[[j]], b_list[[j]], pnew_s[j], log = TRUE)
          ll_lo <- stats::dbinom(m_list[[j]], b_list[[j]], pold_s[j], log = TRUE)
          for(k_i in seq_len(K)){
            L_mat[, k_i] <- L_mat[, k_i] + if(as.logical(md_mat[k_i, j])) ll_hi else ll_lo
          }
        }

        # 12 vectorized EM updates warm-started from EstimateFractions weights
        pi_mle <- pmax(frac_est_mat, 1e-6)
        pi_mle <- pi_mle / rowSums(pi_mle)
        for(em_it in seq_len(12L)){
          log_pi <- log(pmax(pi_mle, 1e-300))
          Z <- log_pi[grp_id, , drop = FALSE] + L_mat
          z_mx <- do.call(pmax, unname(split(Z, col(Z))))
          z_mx_s <- ifelse(is.finite(z_mx), z_mx, 0)
          exp_Z <- exp(Z - z_mx_s)
          sum_ez <- rowSums(exp_Z)
          resp <- exp_Z / pmax(sum_ez, 1e-300)
          for(k_i in seq_len(K)){
            pi_mle[, k_i] <- grp_sum_contig(n_vec * resp[, k_i], grp_end) / N_grp
          }
          pi_mle <- pmax(pi_mle, 0)
          pi_mle <- pi_mle / rowSums(pi_mle)
        }

        p_eff <- pmax(rowSums(pi_mle > 1e-4) - 1, 0)
        if(has_all_low){
          est_fn <- 1 - frac_est_mat[, low_idx_gen]
          mle_fn <- 1 - pi_mle[, low_idx_gen]
        }else{
          est_fn <- frac_est_mat[, 1L]
          mle_fn <- pi_mle[, 1L]
        }
        at_upper_bound <- rep(FALSE, ng)
      }

      # Lookup tables per channel j: P_hi_list[[j]] and P_lo_list[[j]] of size (max_bj + 1) x B1
      P_hi_list <- vector("list", J)
      P_lo_list <- vector("list", J)
      for(j in seq_len(J)){
        max_bj <- max(b_list[[j]], 0L)
        bg_j <- 0:max_bj
        Phi <- matrix(0, nrow = max_bj + 1L, ncol = B1)
        Plo <- matrix(0, nrow = max_bj + 1L, ncol = B1)
        for(km in 0:(cap_m - 1L)){
          Phi[, km + 1L] <- stats::dbinom(km, bg_j, pnew_s[j])
          Plo[, km + 1L] <- stats::dbinom(km, bg_j, pold_s[j])
        }
        Phi[, B1] <- stats::pbinom(cap_m - 1L, bg_j, pnew_s[j], lower.tail = FALSE)
        Plo[, B1] <- stats::pbinom(cap_m - 1L, bg_j, pold_s[j], lower.tail = FALSE)
        P_hi_list[[j]] <- Phi
        P_lo_list[[j]] <- Plo
      }

      # Multi-index grid across C_bins = B1^J
      grid_args <- replicate(J, seq_len(B1), simplify = FALSE)
      cell_grid <- as.matrix(do.call(expand.grid, grid_args))

      E_mat <- matrix(0, nrow = ng, ncol = C_bins)
      for(k_i in seq_len(K)){
        # Skip strictly inactive component across all features (e.g., -label controls)
        if(all(pi_mle[, k_i] == 0)) next
        E_k <- matrix(0, nrow = ng, ncol = C_bins)
        for(c_i in seq_len(C_bins)){
          prob_r <- rep(1, n_rows)
          for(j in seq_len(J)){
            lut_kj <- if(as.logical(md_mat[k_i, j])) P_hi_list[[j]] else P_lo_list[[j]]
            prob_r <- prob_r * lut_kj[b_list[[j]] + 1L, cell_grid[c_i, j]]
          }
          E_k[, c_i] <- grp_sum_contig(n_vec * prob_r, grp_end)
        }
        E_mat <- E_mat + pi_mle[, k_i] * E_k
      }

      raw_diag_mode <- rep(if(s_is_nl) "nolabel_anomaly" else "multilabel_mismatch", ng)
    }

    ### Distributional effect sizes (unpooled, depth-invariant)
    P_obs <- O_mat / N_grp
    P_exp <- E_mat / N_grp
    tvd_vec <- 0.5 * rowSums(abs(P_obs - P_exp))
    ratio_oe <- ifelse(O_mat > 0, O_mat / pmax(E_mat, 1e-300), 1)
    kl_vec <- pmax(rowSums(ifelse(O_mat > 0, P_obs * log(ratio_oe), 0)), 0)

    ### Adaptive bin pooling and likelihood-ratio deviance G^2
    pooled <- pool_gof_bins_vec(O_mat, E_mat, min_bin_count = min_bin_count)
    O_pool <- pooled$O_pool
    E_pool <- pooled$E_pool
    act_bins <- pooled$active
    C_act <- pooled$C_act

    term_g2 <- ifelse(act_bins & (O_pool > 0),
                      2 * O_pool * log(O_pool / pmax(E_pool, 1e-300)),
                      0)
    G2_nom <- pmax(rowSums(term_g2), 0)

    df_vec <- pmax(C_act - 1 - p_eff, 0)
    # Rescue active upper-boundary pinning when 2 pooled bins remain (0 vs 1+)
    if(eff_strategy == "standard" && !s_is_nl){
      # Find highest active bin for C_act == 2 rows
      c2_mask <- (C_act == 2L) & at_upper_bound
      if(any(c2_mask)){
        rows_c2 <- which(c2_mask)
        hi_col <- max.col(act_bins[rows_c2, , drop = FALSE], ties.method = "last")
        lin_hi <- cbind(rows_c2, hi_col)
        pinned_excess <- O_pool[lin_hi] > E_pool[lin_hi]
        if(any(pinned_excess)){
          df_vec[rows_c2[pinned_excess]] <- 1
        }
      }
    }
    df_vec[df_vec < 0.25] <- 0

    ### Empirical sample-wide per-read overdispersion calibration
    elig <- (df_vec >= 0.5) & (N_grp >= max(25, min_reads))
    phi_hat_s <- 0
    if(identical(dispersion_mode, "empirical") && sum(elig) >= 10L){
      med_chisq <- stats::qchisq(0.5, df = df_vec[elig])
      phi_cand <- (df_vec[elig] / N_grp[elig]) * ((G2_nom[elig] / pmax(med_chisq, 1e-12)) - 1)
      phi_hat_s <- max(0, stats::median(phi_cand, na.rm = TRUE))
    }

    scale_fac <- 1 + phi_hat_s * (N_grp / pmax(df_vec, 1e-6))
    G2_adj <- ifelse(df_vec > 0, G2_nom / scale_fac, NA_real_)
    std_dev_vec <- ifelse(df_vec > 0, sqrt(G2_adj / df_vec), NA_real_)
    pval_vec <- ifelse(df_vec > 0, stats::pchisq(G2_adj, df = df_vec, lower.tail = FALSE), NA_real_)
    padj_vec <- rep(NA_real_, ng)
    ok_p <- !is.na(pval_vec)
    if(any(ok_p)){
      padj_vec[ok_p] <- stats::p.adjust(pval_vec[ok_p], method = "BH")
    }

    is_flagged <- !is.na(padj_vec) & (padj_vec < sample_padj_cutoff) & (tvd_vec >= sample_tvd_cutoff)
    f_mode <- ifelse(is_flagged, raw_diag_mode, "ok")

    meta_row_s <- metadf[metadf$sample == s_cur, cond_cols, drop = FALSE]
    res_s <- dplyr::bind_cols(
      dplyr::tibble(sample = rep(s_cur, ng)),
      dplyr::as_tibble(grp_keys_dt),
      dplyr::tibble(tl = rep(unname(tl_disp_vec[s_cur]), ng))
    )
    if(length(cond_cols) > 0L && nrow(meta_row_s) >= 1L){
      res_s <- dplyr::bind_cols(res_s, meta_row_s[rep(1L, ng), , drop = FALSE])
    }
    res_s <- dplyr::bind_cols(
      res_s,
      dplyr::tibble(
        n = N_grp,
        est_fn = est_fn,
        mle_fn = mle_fn,
        tvd = tvd_vec,
        kl_div = kl_vec,
        chisq_gof = G2_adj,
        chisq_nominal = ifelse(df_vec > 0, G2_nom, NA_real_),
        df_gof = df_vec,
        std_dev = std_dev_vec,
        pval_gof = pval_vec,
        padj_gof = padj_vec,
        failure_mode = f_mode,
        .raw_mode = raw_diag_mode,
        is_nolabel = rep(s_is_nl, ng)
      )
    )

    sample_res_list[[si]] <- res_s
    disp_records[[si]] <- dplyr::tibble(
      sample = s_cur,
      is_nolabel = s_is_nl,
      phi_dispersion = phi_hat_s,
      n_features_calibrated = sum(elig)
    )
  }

  all_res <- dplyr::bind_rows(sample_res_list)
  all_res$.g2_s <- all_res$chisq_gof
  all_res$.df_s <- all_res$df_gof
  sample_disp_tbl <- dplyr::bind_rows(disp_records)

  ### Step 4: Aggregate across samples into feature-level summary_tbl
  has_nl_rows <- any(all_res$is_nolabel)

  summary_tbl <- all_res %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(feat_cols))) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_labeled_samples = sum(!is_nolabel),
      n_samples_flagged = sum(failure_mode != "ok"),
      frac_samples_flagged = sum(failure_mode != "ok") / dplyr::n(),
      mean_n = mean(n[!is_nolabel], na.rm = TRUE),
      mean_tvd = mean(tvd[!is_nolabel], na.rm = TRUE),
      max_tvd = if(any(!is_nolabel)) max(tvd[!is_nolabel], na.rm = TRUE) else NA_real_,
      mean_kl_div = mean(kl_div[!is_nolabel], na.rm = TRUE),
      chisq_gof = {
        ok_c <- (!is_nolabel) & (.df_s > 0) & is.finite(.g2_s)
        if(any(ok_c)) sum(.g2_s[ok_c]) else NA_real_
      },
      df_gof = {
        ok_c <- (!is_nolabel) & (.df_s > 0) & is.finite(.g2_s)
        if(any(ok_c)) sum(.df_s[ok_c]) else 0
      },
      pval_gof = {
        ok_c <- (!is_nolabel) & (.df_s > 0) & is.finite(.g2_s)
        if(any(ok_c)) stats::pchisq(sum(.g2_s[ok_c]), df = sum(.df_s[ok_c]), lower.tail = FALSE) else NA_real_
      },
      .raw_dom_cand = {
        lb_idx <- which(!is_nolabel)
        if(length(lb_idx) > 0L) .raw_mode[lb_idx[which.max(tvd[lb_idx])[1L]]] else .raw_mode[1L]
      },
      .flagged_dom = {
        fm <- failure_mode[failure_mode != "ok"]
        if(length(fm) > 0L){
          tb <- table(fm)
          names(tb)[which.max(tb)]
        }else{
          NA_character_
        }
      },
      nolabel_tvd = if(any(is_nolabel)) mean(tvd[is_nolabel], na.rm = TRUE) else NA_real_,
      nolabel_pval = {
        ok_nl <- is_nolabel & (.df_s > 0) & is.finite(.g2_s)
        if(any(ok_nl)) stats::pchisq(sum(.g2_s[ok_nl]), df = sum(.df_s[ok_nl]), lower.tail = FALSE) else NA_real_
      },
      .groups = "drop"
    )

  summary_tbl$padj_gof <- rep(NA_real_, nrow(summary_tbl))
  ok_pg <- !is.na(summary_tbl$pval_gof)
  if(any(ok_pg)){
    summary_tbl$padj_gof[ok_pg] <- stats::p.adjust(summary_tbl$pval_gof[ok_pg], method = "BH")
  }

  if(has_nl_rows){
    summary_tbl$nolabel_padj <- rep(NA_real_, nrow(summary_tbl))
    ok_pnl <- !is.na(summary_tbl$nolabel_pval)
    if(any(ok_pnl)){
      summary_tbl$nolabel_padj[ok_pnl] <- stats::p.adjust(summary_tbl$nolabel_pval[ok_pnl], method = "BH")
    }
  }else{
    summary_tbl$nolabel_tvd <- NULL
    summary_tbl$nolabel_pval <- NULL
  }

  comb_flagged <- !is.na(summary_tbl$padj_gof) &
    (summary_tbl$padj_gof < sample_padj_cutoff) &
    is.finite(summary_tbl$mean_tvd) &
    (summary_tbl$mean_tvd >= sample_tvd_cutoff)

  summary_tbl$dominant_failure_mode <- ifelse(
    !is.na(summary_tbl$.flagged_dom),
    summary_tbl$.flagged_dom,
    ifelse(comb_flagged, summary_tbl$.raw_dom_cand, "ok")
  )
  summary_tbl$.raw_dom_cand <- NULL
  summary_tbl$.flagged_dom <- NULL

  # Order columns cleanly
  base_sum_cols <- c(feat_cols, "n_samples", "n_labeled_samples",
                     "n_samples_flagged", "frac_samples_flagged",
                     "mean_n", "mean_tvd", "max_tvd", "mean_kl_div",
                     "chisq_gof", "df_gof", "pval_gof", "padj_gof")
  nl_sum_cols <- if(has_nl_rows) c("nolabel_tvd", "nolabel_pval", "nolabel_padj") else character(0)
  summary_tbl <- summary_tbl[, c(base_sum_cols, nl_sum_cols, "dominant_failure_mode"), drop = FALSE]

  residuals_tbl <- all_res %>%
    dplyr::select(-dplyr::all_of(c(".raw_mode", "is_nolabel", ".g2_s", ".df_s")))

  list(
    summary_tbl = summary_tbl,
    residuals_tbl = residuals_tbl,
    metadata = list(
      source_type = "fractions",
      source_table = source_table,
      features = feat_cols,
      populations = fxn_meta[["populations"]],
      fraction_design = fxn_meta[["fraction_design"]],
      strategy = eff_strategy,
      min_reads = min_reads,
      min_bin_count = min_bin_count,
      sample_padj_cutoff = sample_padj_cutoff,
      sample_tvd_cutoff = sample_tvd_cutoff,
      sample_dispersion = sample_disp_tbl
    )
  )
}


# Reconstruct exact observed & fitted mutation-count histograms on the fly for EZFitPlot()
predict_mixture_histograms <- function(obj,
                                       gof_name,
                                       selected_keys_df,
                                       samples = NULL,
                                       rank_by = "worst_lof",
                                       max_mut_plot = 10L){
  sample <- NULL

  gof_meta <- obj[["metadata"]][["gof"]][[gof_name]]
  res_tbl <- as.data.frame(obj[["gof_residuals"]][[gof_name]])
  source_table <- gof_meta[["source_table"]]
  fxn_tbl <- as.data.frame(obj[["fractions"]][[source_table]])
  fxn_meta <- obj[["metadata"]][["fractions"]][[source_table]]
  feat_cols <- gof_meta[["features"]]
  eff_strategy <- gof_meta[["strategy"]] %||% "standard"

  fraction_design <- as.data.frame(fxn_meta[["fraction_design"]])
  mutrate_design <- fraction_design[as.logical(fraction_design$present),
                                    setdiff(colnames(fraction_design), "present"),
                                    drop = FALSE]
  mutcols <- colnames(mutrate_design)
  basecols <- paste0("n", substr(mutcols, 1L, 1L))
  K <- nrow(mutrate_design)
  J <- ncol(mutrate_design)
  md_mat <- as.matrix(mutrate_design)

  metadf <- as.data.frame(obj[["metadf"]])
  metadf$sample <- as.character(metadf$sample)
  tl_cols_possible <- if(J == 1L){
    c("tl", "tpulse")
  }else{
    c(paste0("tl_", mutcols), paste0("tpulse_", mutcols), "tl", "tpulse")
  }
  tl_cols <- intersect(tl_cols_possible, colnames(metadf))
  tl_mat <- as.matrix(metadf[, tl_cols, drop = FALSE])
  storage.mode(tl_mat) <- "double"
  samples_nolabel <- metadf$sample[apply(tl_mat == 0, 1L, all)]

  for(fc in feat_cols){
    res_tbl[[fc]] <- as.character(res_tbl[[fc]])
    selected_keys_df[[fc]] <- as.character(selected_keys_df[[fc]])
    fxn_tbl[[fc]] <- as.character(fxn_tbl[[fc]])
  }
  fxn_tbl$sample <- as.character(fxn_tbl$sample)
  res_tbl$sample <- as.character(res_tbl$sample)

  res_sub <- dplyr::inner_join(res_tbl, selected_keys_df, by = feat_cols)
  if(nrow(res_sub) == 0L){
    stop("No residual records found for the requested features!")
  }

  # Choose which (feature, sample) pairs to plot
  if(!is.null(samples)){
    if(!identical(samples, "all")){
      samps_req <- as.character(samples)
      res_sub <- res_sub[res_sub$sample %in% samps_req, , drop = FALSE]
      if(nrow(res_sub) == 0L){
        stop("None of the requested `samples` were found for the selected features!")
      }
    }
    pairs_df <- res_sub
  }else{
    # Pick 1 representative sample per selected feature
    pair_rows <- vector("list", nrow(selected_keys_df))
    for(i in seq_len(nrow(selected_keys_df))){
      rk <- selected_keys_df[i, feat_cols, drop = FALSE]
      cand <- dplyr::inner_join(res_sub, rk, by = feat_cols)
      if(rank_by == "best_fit"){
        cand_lb <- cand[!(cand$sample %in% samples_nolabel), , drop = FALSE]
        if(nrow(cand_lb) > 0L) cand <- cand_lb
        idx_pick <- order(cand$tvd, -cand$n)[1L]
      }else if(rank_by == "worst_rep" && any(cand$sample %in% samples_nolabel)){
        cand_nl <- cand[cand$sample %in% samples_nolabel, , drop = FALSE]
        idx_pick <- order(-cand_nl$tvd, cand_nl$pval_gof)[1L]
        cand <- cand_nl
      }else{
        cand_lb <- cand[!(cand$sample %in% samples_nolabel), , drop = FALSE]
        base_c <- if(nrow(cand_lb) > 0L) cand_lb else cand
        flagged_c <- base_c[base_c$failure_mode != "ok", , drop = FALSE]
        pool_c <- if(nrow(flagged_c) > 0L) flagged_c else base_c
        idx_pick <- order(-pool_c$tvd, pool_c$pval_gof)[1L]
        cand <- pool_c
      }
      pair_rows[[i]] <- cand[idx_pick, , drop = FALSE]
    }
    pairs_df <- dplyr::bind_rows(pair_rows)
  }

  ext <- extract_cb_for_gof(
    obj = obj,
    samples_keep = unique(pairs_df$sample),
    samples_nolabel = samples_nolabel,
    feat_cols = feat_cols,
    mutcols = mutcols,
    basecols = basecols,
    fxn_meta = fxn_meta,
    features_keep_df = selected_keys_df
  )
  cb_dt <- ext$cb_dt
  pold_nolabel <- ext$pold_nolabel

  hier_cand_name <- fxn_meta[["hier_mutrate_table"]] %||%
    paste0(paste(feat_cols, collapse = "_"), "_", mutcols[1L])

  if(eff_strategy == "general"){
    pop_labels <- vapply(seq_len(K), function(k_i){
      pv <- ifelse(as.logical(md_mat[k_i, ]), "high", "low")
      paste(paste0(pv, mutcols), collapse = "_")
    }, character(1))
    frac_cols_gen <- paste0("fraction_", pop_labels)
  }

  obs_out <- list()
  curve_out <- list()

  for(r_i in seq_len(nrow(pairs_df))){
    pr <- pairs_df[r_i, , drop = FALSE]
    s_cur <- pr$sample[1L]
    s_is_nl <- s_cur %in% samples_nolabel
    fk_df <- pr[, feat_cols, drop = FALSE]

    sub_dt <- data.table::as.data.table(fk_df)
    dt_pair <- cb_dt[sample == s_cur][sub_dt, on = feat_cols, nomatch = NULL]
    if(nrow(dt_pair) == 0L) next
    n_r <- as.numeric(dt_pair[["n"]])
    N_tot <- sum(n_r)

    for(j in seq_len(J)){
      mut_j <- mutcols[j]
      base_j <- basecols[j]
      mj <- as.integer(round(as.numeric(dt_pair[[mut_j]])))
      bj <- as.integer(round(as.numeric(dt_pair[[base_j]])))

      m_cap <- as.integer(max(6L, min(as.integer(max_mut_plot), max(mj, 6L))))
      n_bins <- m_cap + 1L
      bin_labels <- c(as.character(0:(m_cap - 1L)), paste0(m_cap, "+"))
      m_binned <- pmin(pmax(mj, 0L), m_cap)

      obs_prop <- vapply(0:m_cap, function(mm){
        sum(n_r[m_binned == mm]) / N_tot
      }, numeric(1))

      if(s_is_nl){
        po_j <- if(!is.null(pold_nolabel[[s_cur]])) pold_nolabel[[s_cur]][mut_j] else 1e-4
        po_j <- pmin(pmax(po_j, 1e-7), 1 - 1e-7)
        pn_j <- po_j
        pi_hi <- 0
      }else if(eff_strategy == "hierarchical" && J == 1L){
        hmr <- as.data.frame(obj[["mutation_rates"]][[hier_cand_name]])
        hmr$sample <- as.character(hmr$sample)
        for(fc in feat_cols) hmr[[fc]] <- as.character(hmr[[fc]])
        hmr_m <- dplyr::inner_join(hmr[hmr$sample == s_cur, , drop = FALSE], fk_df, by = feat_cols)
        pn_j <- pmin(pmax(as.numeric(hmr_m$pnew[1L]), 1e-7), 1 - 1e-7)
        po_j <- pmin(pmax(as.numeric(hmr_m$pold[1L]), 1e-7), 1 - 1e-7)
        pi_hi <- as.numeric(pr$mle_fn[1L])
      }else if(eff_strategy == "standard" && J == 1L){
        mr_j <- as.data.frame(obj[["mutation_rates"]][[mut_j]])
        mr_m <- mr_j[as.character(mr_j$sample) == s_cur, , drop = FALSE]
        pn_j <- pmin(pmax(as.numeric(mr_m$pnew[1L]), 1e-7), 1 - 1e-7)
        po_j <- pmin(pmax(as.numeric(mr_m$pold[1L]), 1e-7), 1 - 1e-7)
        pi_hi <- as.numeric(pr$mle_fn[1L])
      }else{
        mr_j <- as.data.frame(obj[["mutation_rates"]][[mut_j]])
        mr_m <- mr_j[as.character(mr_j$sample) == s_cur, , drop = FALSE]
        if(nrow(mr_m) > 0L && is.finite(mr_m$pnew[1L])){
          pn_j <- pmin(pmax(as.numeric(mr_m$pnew[1L]), 1e-7), 1 - 1e-7)
          po_j <- pmin(pmax(as.numeric(mr_m$pold[1L]), 1e-7), 1 - 1e-7)
        }else{
          po_j <- pmin(pmax(sum(mj * n_r) / max(sum(bj * n_r), 1), 1e-6), 0.25)
          pn_j <- po_j
        }
        fxn_row <- dplyr::inner_join(fxn_tbl[fxn_tbl$sample == s_cur, , drop = FALSE], fk_df, by = feat_cols)
        hi_k <- which(as.logical(md_mat[, j]))
        pi_hi <- if(length(hi_k) > 0L && nrow(fxn_row) > 0L){
          sum(as.numeric(fxn_row[1L, frac_cols_gen[hi_k]]))
        }else{
          0
        }
      }

      p_old_vec <- numeric(n_bins)
      p_new_vec <- numeric(n_bins)
      for(mm in 0:(m_cap - 1L)){
        p_old_vec[mm + 1L] <- sum(n_r * stats::dbinom(mm, bj, po_j)) / N_tot
        p_new_vec[mm + 1L] <- sum(n_r * stats::dbinom(mm, bj, pn_j)) / N_tot
      }
      p_old_vec[n_bins] <- sum(n_r * stats::pbinom(m_cap - 1L, bj, po_j, lower.tail = FALSE)) / N_tot
      p_new_vec[n_bins] <- sum(n_r * stats::pbinom(m_cap - 1L, bj, pn_j, lower.tail = FALSE)) / N_tot

      comp_old <- (1 - pi_hi) * p_old_vec
      comp_new <- pi_hi * p_new_vec
      comp_tot <- comp_old + comp_new

      base_info <- dplyr::bind_cols(
        fk_df[rep(1L, n_bins), , drop = FALSE],
        dplyr::tibble(
          sample = rep(s_cur, n_bins),
          muttype = rep(mut_j, n_bins),
          mut_bin = 0:m_cap,
          mut_label = factor(bin_labels, levels = bin_labels),
          tvd = rep(pr$tvd[1L], n_bins),
          failure_mode = rep(pr$failure_mode[1L], n_bins)
        )
      )

      obs_out[[length(obs_out) + 1L]] <- dplyr::bind_cols(
        base_info,
        dplyr::tibble(obs_prop = obs_prop)
      )

      if(s_is_nl){
        curve_out[[length(curve_out) + 1L]] <- dplyr::bind_cols(
          base_info,
          dplyr::tibble(
            component = factor(rep("Total mixture", n_bins),
                               levels = c("Total mixture", "Old component", "New component")),
            pred_prop = comp_tot
          )
        )
      }else{
        curve_out[[length(curve_out) + 1L]] <- dplyr::bind_rows(
          dplyr::bind_cols(base_info, dplyr::tibble(
            component = "Total mixture", pred_prop = comp_tot
          )),
          dplyr::bind_cols(base_info, dplyr::tibble(
            component = "Old component", pred_prop = comp_old
          )),
          dplyr::bind_cols(base_info, dplyr::tibble(
            component = "New component", pred_prop = comp_new
          ))
        )
        curve_out[[length(curve_out)]]$component <- factor(
          curve_out[[length(curve_out)]]$component,
          levels = c("Total mixture", "Old component", "New component")
        )
      }
    }
  }

  list(
    obs_df = dplyr::bind_rows(obs_out),
    curve_df = dplyr::bind_rows(curve_out),
    pairs_df = pairs_df,
    J = J
  )
}

