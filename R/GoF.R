#' Assess goodness of fit of kinetic or dynamical systems models
#'
#' `AssessGoF()` evaluates how well model-predicted sample-by-sample fractions of
#' high-mutation reads (`fraction_high<muttype>`) from either `AverageAndRegularize()`
#' (e.g., standard steady-state turnover \eqn{dR/dt = k_{\text{syn}} - k_{\text{deg}} R})
#' or `EZDynamics()` (generalized linear ODE reaction networks) agree with the
#' observed sample-level estimates and profile-likelihood uncertainties from
#' `EstimateFractions()`.
#'
#' For each feature, `AssessGoF()` computes both interpretable probability-scale
#' effect sizes (precision-weighted root-mean-square error `wrmse_fn`, mean absolute
#' error `mae_fn`, and signed bias `mean_signed_resid_fn`) and uncertainty-standardized
#' goodness-of-fit statistics that cleanly handle boundary observations
#' (\eqn{\hat{\theta} \in \{0, 1\}}) using the 95\% profile-likelihood confidence bounds
#' from `EstimateFractions()`. When experimental cells contain biological replicates,
#' total lack-of-fit is orthogonally decomposed into:
#' \itemize{
#'   \item \strong{Replicate overdispersion} (`chisq_rep`, `df_rep`, `pval_rep`, `padj_rep`):
#'   tests whether biological replicates within identical experimental conditions and
#'   label times scatter more than expected from `EstimateFractions()` measurement
#'   uncertainty alone.
#'   \item \strong{Structural model lack-of-fit} (`chisq_lof`, `df_lof`, `pval_lof`, `padj_lof`):
#'   tests whether condition/timepoint/compartment means systematically deviate from
#'   the fitted kinetic trajectory (for example, multi-exponential decay, RNA processing
#'   or export delays, premature nuclear degradation, or non-additive covariate effects),
#'   calibrated against Empirical-Bayes-moderated replicate dispersion by default.
#' }
#'
#' Summary statistics per feature are stored in `obj$gof[[output_name]]` and
#' sample-by-sample observed vs. predicted values and residuals are stored in
#' `obj$gof_residuals[[output_name]]` (when `store_residuals = TRUE`), both accessible
#' via `EZget(obj, type = "gof")` and `EZget(obj, type = "gof_residuals")`.
#'
#' @param obj An `EZbakRData` (or `EZbakRArrowData`) object containing output from
#' `EstimateFractions()` as well as either `AverageAndRegularize()` or `EZDynamics()`.
#' @param type Which upstream model class to assess: `"auto"` (default; automatically
#' selects `"dynamics"` or `"averages"` when unambiguous), `"averages"`, or `"dynamics"`.
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
#' @param repeatID Optional numerical repeat ID of the upstream `"averages"` or `"dynamics"` table.
#' @param fractions_repeatID Optional numerical repeat ID of the `EstimateFractions()` table to compare against.
#' @param exactMatch Logical; passed to `EZget()` (default `TRUE`).
#' @param fraction_col Optional name of the `fraction_high<muttype>` column in the
#' `fractions` table to evaluate (auto-detected when a single high-mutation population exists).
#' @param dispersion_mode How replicate dispersion is accounted for when testing structural
#' lack-of-fit (`pval_lof`): `"shrunk"` (default; Empirical-Bayes trend-moderated dispersion
#' across replicates when available), `"global"` (single dataset-wide median replicate
#' overdispersion factor), or `"measurement_only"` (uses only `EstimateFractions()` profile
#' uncertainty plus `logit_sd_floor`).
#' @param logit_sd_floor Baseline technical standard deviation floor on the logit fraction
#' scale (default `0.05`) preventing ultra-high-coverage housekeeping genes from being
#' flagged for negligible technical jitter.
#' @param include_coverage_gof Logical (default `TRUE`); when assessing an `EZDynamics()`
#' fit that used `use_coverage = TRUE`, also evaluate steady-state normalized read coverage
#' goodness-of-fit (`wrmse_log10_cov`, `chisq_cov`, `pval_lof_joint`, `padj_lof_joint`).
#' @param store_residuals Logical (default `TRUE`); if `TRUE`, stores sample-by-sample
#' predictions and residuals in `obj$gof_residuals[[output_name]]` (required for
#' `EZFitPlot()` and sample-level `EZGoFPlot()` modes).
#' @param character_limit Maximum length of automatically generated output table names.
#' @param overwrite Logical (default `TRUE`); whether to overwrite an existing GoF table
#' with identical metadata.
#' @return Input `EZbakRData` object with class `"EZbakRGoF"` prepended, containing
#' feature-level goodness-of-fit metrics in `obj$gof` and (optionally) sample-level
#' predictions/residuals in `obj$gof_residuals`.
#' @import data.table
#' @importFrom magrittr %>%
#' @examples
#' simdata <- EZSimulate(30, nreps = 2)
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#' ezbdo <- EstimateFractions(ezbdo)
#' ezbdo <- EstimateKinetics(ezbdo)
#' ezbdo <- AverageAndRegularize(ezbdo)
#' ezbdo <- AssessGoF(ezbdo)
#' gof_tbl <- EZget(ezbdo, type = "gof")
#' @export
AssessGoF <- function(obj,
                      type = c("auto", "averages", "dynamics"),
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
                      dispersion_mode = c("shrunk", "global", "measurement_only"),
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

  if(!is.null(formula_mean)){
    environment(formula_mean) <- NULL
  }

  ### Step 1: Resolve upstream model type ("averages" vs "dynamics")
  has_dyn <- length(obj[["metadata"]][["dynamics"]]) > 0L
  has_avg <- length(obj[["metadata"]][["averages"]]) > 0L

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
        stop("Both `averages` and `dynamics` tables exist in `obj`! Please specify `type = 'averages'` or `type = 'dynamics'`.")
      }
    }else{
      stop("`obj` has neither `averages` nor `dynamics` tables! Run `AverageAndRegularize()` or `EZDynamics()` before `AssessGoF()`.")
    }
  }else{
    model_type <- type
    if(model_type == "averages" && !has_avg){
      stop("No `averages` tables found in `obj`! Run `AverageAndRegularize()` first.")
    }
    if(model_type == "dynamics" && !has_dyn){
      stop("No `dynamics` tables found in `obj`! Run `EZDynamics()` first.")
    }
  }

  ### Step 2: Locate upstream table and run forward predictions
  if(model_type == "averages"){
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
      dispersion_mode = dispersion_mode,
      logit_sd_floor = logit_sd_floor
    )

  }else{
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
      dispersion_mode = dispersion_mode,
      logit_sd_floor = logit_sd_floor,
      include_coverage_gof = include_coverage_gof
    )
  }

  ### Step 3: Name and store outputs in obj$gof and obj$gof_residuals
  proposed_name <- paste0("gof_", source_table)
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
  gof_meta$dispersion_mode <- dispersion_mode
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
