#' Estimate kinetic parameters
#'
#' Simple kinetic parameter (kdeg and ksyn) estimation using fraction estimates
#' from `EstimateFractions()`. Several slightly different kinetic parameter
#' inference strategies are implemented.
#'
#' `EstimateKinetics()` estimates the kinetics and RNA synthesis and degradation
#' from standard single-label NR-seq data. It also technically supports analyses
#' of the dual-label TILAC experiment, but that functionality is not as actively
#' tested or maintained as the more standard analyses.
#'
#' `EstimateKinetics()` assumes a simple linear ODE model of RNA dynamics:
#' \deqn{\frac{\text{dR}}{\text{dt}} = k_{\text{syn}} - k_{\text{deg}}*\text{R}}
#' where:
#'  - R = Amount of RNA
#'  - \eqn{k_{\text{syn}}} = Synthesis rate constant
#'  - \eqn{k_{\text{deg}}} = Degradation rate constant
#' and for which the general solution is:
#' \deqn{\text{R(t)} = \text{R(0)}e^{-k_{\text{deg}}\text{t}} + (1 - \frac{k_{\text{syn}}}{k_{\text{deg}}})e^{-k_{\text{deg}}\text{t}}  }
#' where \eqn{\text{R(0)}} is the initial RNA level.
#'
#' The default kinetic parameter estimation strategy (`strategy == "standard"`)
#' assumes that for labeled RNA (or more precisely RNA synthesized in the presence
#' of label) \eqn{\text{R(0)} = 0}. Thus, it assumes a pulse-label rather than
#' pulse-chase experimental design (I've written several places, like [here](https://simonlabcode.github.io/bakR/articles/Troubleshooting.html#experimental-suggestions)
#' and [here](https://isaacvock.github.io/EZbakR/articles/EstimateKinetics.html#pulse-chase-analyses)
#' for example, about why pulse-label designs are superior to pulse-chases in
#' almost all settings).
#'
#' For unlabeled RNA, it assumes that \eqn{\text{R(0)} = \frac{k_{\text{syn}}}{k_{\text{deg}}}}.
#' This is known as the steady-state assumption and is the key assumption of this method.
#' More broadly, assuming steady-state means that this method assumes that RNA levels
#' from a given feature are constant for the entire metabolic labeling duration.
#' As EZbakR is designed for analyses of bulk NR-seq data, "constant" means that the
#' average RNA levels across all profiled cells is constant (thus asynchronous populations
#' of dividing cells still count as steady-state, even if the RNA expression landscapes
#' in individual cells are quite dynamic). This assumption may be violated in cases
#' where labeling coincides with or closely follows a perturbation (e.g., drug treatment).
#' When the steady-state assumption hold, there is a simple relationship between the
#' fraction of reads from labeled RNA (\eqn{\theta}) and the turnover kinetics of the RNA:
#' \deqn{\theta = 1 - e^{-k_{\text{deg}}t}}
#' The power in this approach is thus that turnover kinetics are estimated from
#' an "internally normalized" quantity: \eqn{\theta} (termed the "fraction new",
#' "fraction labeled", "fraction high mutation content", or new-to-total ratio (NTR), depending
#' on where you look or who you ask). "Internally normalized" means that a normalization
#' scale factor does not need to be estimated in order to accurately infer turnover
#' kinetics. \eqn{\theta} represents a ratio of read counts from the same feature
#' in the same library, thus any scale factor would appear in both the numerator
#' and denominator of this estimate and cancel out.
#'
#' When this assumption is valid, the `strategy = "NSS"` approach may be preferable.
#' This approach relies on the same model of RNA dynamics, but now assumes that the
#' initial RNA levels (\eqn{\text{R(0)}} for the unlabeled RNA) are not at the steady-state
#' levels dictated by the current synthesis and turnover kinetics. The idea for this
#' strategy was first presented in [Narain et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8354102/).
#' To run this approach, you need to be able to estimate the initial RNA levels of
#' each label-fed sample, as the fraction of reads from labeled RNA will no longer.
#' only reflect the turnover kinetics (as the old RNA is assumed to potentially
#' not have reached the new steady-state levels yet). This means that for each
#' label-fed sample, you need to have a sample whose assayed RNA population represents
#' the starting RNA population for the label-fed sample. EZbakR will use these
#' reference samples to infer \eqn{\text{R(0)}} for unlabeled RNAs and estimate
#' turnover kinetics from the ratio of this quantity to the remaining unlabeled RNA levels
#' after labeling:
#' \deqn{\theta_{\text{NSS}} = \frac{\text{R(tl)}}{\text{R(0)}}}
#' While I really like the idea of this strategy, it is not without some severe limitations.
#' For one, the major benefit of NR-seq, internally normalized estimation of turnover
#' kinetics, is out the window. Thus, read counts need to be normalized between
#' the relevant reference and label-fed sample pairs. In addition, the variance
#' patterns of this ratio are somewhat unfortunate. Its variance
#' is incredibly high for more stable RNAs, severely limiting the effective
#' dynamic range of this approach relative to steady-state analyses. I continue
#' to work to refine EZbakR's implementation of this approach, but for now I
#' urge caution in its use and interpretation. See my discussion of an alternative
#' approach for navigating non-steady-state data [here](https://simonlabcode.github.io/bakR/articles/NSS.html#going-beyond-the-steady-state-assumption).
#'
#' As an aside, you may wonder how this strategy deals with dynamic regulation of
#' synthesis and degradation rate constants *during labeling*. To answer this,
#' you have to realize that the duration of metabolic labeling represents an
#' integration time over which the best we can do is assess average kinetics. Thus,
#' if rate constants are changing during the labeling, this strategy can still be
#' thought of as providing a an estimate of the time averaged synthesis and turnover
#' kinetics during the label time.
#'
#' Eventually, I will get around to implementing pulse-chase support in `EstimateKinetics()`.
#' I haven't yet because I don't like pulse-chase experiments and think they are
#' way more popular than they should be for purely historical reasons. But lots
#' of people keep doing pulse-chase NR-seq so c'est la vie.
#'
#' @param obj An `EZbakRFractions` object, which is an `EZbakRData` object on
#' which `EstimateFractions()` has been run.
#' @param strategy Kinetic parameter estimation strategy.
#' Options include:
#' \itemize{
#'  \item standard: Estimate a single new read and old read mutation rate for each
#'  sample. This is done via a binomial mixture model aggregating over
#'  \item NSS: Use strategy similar to that presented in Narain et al., 2021 that
#'  assumes you have data which provides a reference for how much RNA was present
#'  at the start of each labeling. In this case, `grouping_factor` and `reference_factor`
#'  must also be set.
#'  \item shortfeed: Estimate kinetic parameters assuming no degradation of labeled
#'  RNA, most appropriate if the metabolic label feed time is much shorter than
#'  the average half-life of an RNA in your system.
#'  \item tilac: Estimate TILAC-ratio as described in Courvan et al., 2022.
#'  \item pulse-chase: Estimate kdeg for a pulse-chase experiment. By default kdeg will be estimated
#'  for each time point at which label was present. This includes any pulse-only samples,
#'  as well as all samples including a chase after the pulse. If you don't want
#'  to include the estimates from the pulse-only samples in the final output,
#'  set `exclude_pulse_estimates` to `TRUE`. Pulse-chases are suboptimal for a number
#'  of experimental reasons, so we urge users to avoid performing this kind of experiment
#'  whenever possible (favoring instead a pulse-label design). One of the challenges of
#'  analyzing pulse-chase data is that the fraction labeled after the pulse must be compared
#'  to that after each chase. Due to a number of technical reasons, it is possible for the
#'  estimated labeling after the chase to be **higher* than that after the pulse. It is
#'  impossible to estimate kinetic parameters in this case. Thus, when this arises,
#'  a conservative kdeg estimate is imputed, equal to -log(1 - 1/(n+2))/tchase,
#'  which is the kdeg estimate you would get if you had no detected reads from labeled
#'  RNA and an uninformative prior on the fraction new (i.e., the estimated fraction
#'  new is (number of reads from labeled RNA + 1) / (number of reads + 2)).
#' }
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
#' @param repeatID If multiple `fractions` tables exist with the same metadata,
#' then this is the numerical index by which they are distinguished.
#' @param exactMatch If TRUE, then `features` and `populations` have to exactly match
#' those for a given fractions table for that table to be used. Means that you can't
#' specify a subset of features or populations by default, since this is TRUE
#' by default.
#' @param exclude_pulse_estimates If `strategy = "pulse-chase"`, would you like
#' to exclude the pulse-only kdeg and ksyn estimates from the final output? This
#' is a good idea if you used a very long pulse with the goal of labeling close
#' to 100% of all RNA.
#' @param grouping_factor Which sample-detail columns in the metadf should be used
#' to group samples by for calculating the average RPM (`strategy = "NSS"`) or
#' pulse fraction labeled (`strategy = "pulse-chase"`) at a particular time point?
#' Only relevant if `strategy = "NSS"` or `strategy = "pulse-chase"`.
#' @param reference_factor Which sample-detail column in the metafd should be used to
#' determine which group of samples provide information about the starting levels of RNA
#' for each sample? This should have values that match those in `grouping_factor`.
#' #' Only relevant if `strategy = "NSS"`.
#' @param character_limit Maximum number of characters for naming out fractions output. EZbakR
#' will try to name this as a "_" separated character vector of all of the features analyzed.
#' If this name is greater than `character_limit`, then it will default to "fraction#", where
#' "#" represents a simple numerical ID for the table.
#' @param feature_lengths Table of effective lengths for each feature combination in your
#' data. For example, if your analysis includes features named GF and XF, this
#' should be a data frame with columns GF, XF, and length.
#' @param scale_factors Data frame with two columns, "sample" and "scale_factor".
#' sample should be the sample name, and scale_factor is the factor which read counts
#' in that sample should be **divided** to get the normalized read counts. This should match
#' the convention from DESeq2 and thus you can provide, for example, scale factors derived
#' from DESeq2 for this. By default, EZbakR uses DESeq2's median of ratios normalization
#' method, but specifying `scale_factors` is useful when you have spike-ins.
#' @param overwrite If TRUE and a fractions estimate output already exists that
#' would possess the same metadata (features analyzed, populations analyzed,
#' and fraction_design), then it will get overwritten with the new output. Else,
#' it will be saved as a separate output with the same name + "_#" where "#" is a
#' numerical ID to distinguish the similar outputs.
#' @return `EZbakRKinetics` object, which is just an `EZbakRData` object with
#' a "kinetics" slot. This includes tables of kinetic parameter estimates for
#' each feature in each sample for which kinetic parameters can be estimated.
#' @import data.table
#'
#' @examples
#'
#' # Simulate data to analyze
#' simdata <- SimulateOneRep(30)
#'
#' # Create EZbakR input
#' metadf <- data.frame(sample = "sampleA", tl = 2)
#' ezbdo <- EZbakRData(simdata$cB, metadf)
#'
#' # Estimate Fractions
#' ezbdo <- EstimateFractions(ezbdo)
#'
#' # Estimate Kinetics
#' ezbdo <- EstimateKinetics(ezbdo)
#'
#' @export
EstimateKinetics <- function(obj,
                             strategy = c("standard", "tilac", "NSS",
                                          "shortfeed",
                                          "pulse-chase"),
                             features = NULL,
                             populations = NULL,
                             fraction_design = NULL,
                             repeatID = NULL,
                             exactMatch = TRUE,
                             grouping_factor = NULL,
                             reference_factor = NULL,
                             character_limit = 20,
                             feature_lengths = NULL,
                             exclude_pulse_estimates = TRUE,
                             scale_factors = NULL,
                             overwrite = TRUE){

  ### Check that input is valid

  # EZbakRData object on which EstimateFractions() has been run?
  if(!methods::is(obj, "EZbakRFractions")){

    if(methods::is(obj, "EZbakRData") || methods::is(obj, "EZbakRArrowData")){

      stop("obj is not an EZbakRFractions object! Run `obj <- EstimateFractions(obj, ...)`,
           where ... represents optional parameters, before running `EstimateKinetics`.")

    }else{

      stop("obj is not an EZbakRData object with fraction estimates!")

    }


  }

  # Analysis strategy
  strategy <- match.arg(strategy)

  # features
  if(!is.null(features)){

    if(!is.character(features)){

      stop("features is not a character vector!")

    }

  }




  ### "Method dispatch"

  if(strategy == "tilac"){

    obj <- tilac_ratio_estimation(obj,
                                  features = features,
                                  populations = populations,
                                  fraction_design = fraction_design,
                                  repeatID = repeatID,
                                  exactMatch = exactMatch,
                                  character_limit = character_limit,
                                  feature_lengths = feature_lengths,
                                  scale_factors = scale_factors,
                                  overwrite = overwrite)


  }else{

    obj <- Standard_kinetic_estimation(obj,
                                       strategy = strategy,
                                       features = features,
                                       populations = populations,
                                       fraction_design = fraction_design,
                                       reference_factor = reference_factor,
                                       grouping_factor = grouping_factor,
                                       repeatID = repeatID,
                                       exactMatch = exactMatch,
                                       character_limit = character_limit,
                                       feature_lengths = feature_lengths,
                                       scale_factors = scale_factors,
                                       exclude_pulse_estimates = exclude_pulse_estimates,
                                       overwrite = overwrite)

  }


  return(obj)

}


# kdeg = -log(1 - fn)/tl
# ksyn = (normalized read count)*kdeg
Standard_kinetic_estimation <- function(obj,
                                        strategy = c("standard", "NSS",
                                                     "shortfeed",
                                                     "pulse-chase"),
                                        features = NULL,
                                        populations = NULL,
                                        fraction_design = NULL,
                                        repeatID = NULL,
                                        exactMatch = TRUE,
                                        character_limit = 20,
                                        reference_factor = NULL,
                                        grouping_factor = NULL,
                                        feature_lengths = NULL,
                                        scale_factors = NULL,
                                        exclude_pulse_estimates = TRUE,
                                        overwrite = TRUE){

  ### Hack to deal with devtools::check() NOTEs
  tl <- kdeg <- log_kdeg <- se_log_kdeg <- ..cols_to_keep <- ..kinetics_cols_to_keep <- NULL
  ksyn <- normalized_reads <- log_ksyn <- se_log_ksyn <- scale_factor <- n <- nolabel_rpm <- NULL
  old_rpm <- new_rpm <- geom_mean <- rpm <- nolabel_n <- nolabel_reps <- NULL
  reference_rpm <- reference_n <- reference_reps <- tpulse <- tchase <- NULL
  fn_safe <- se_safe <- NULL


  rm(..cols_to_keep)
  rm(..kinetics_cols_to_keep)

  ci_lower_kdeg <- ci_upper_kdeg <- thalf <- ci_lower_thalf <- ci_upper_thalf <- NULL
  ci_lower_ksyn <- ci_upper_ksyn <- fn_lo_tmp <- fn_hi_tmp <- NULL

  `.` <- list


  strategy <- match.arg(strategy)

  ### Figure out which fraction new estimates to use

  # Function is in Helpers.R
  fractions_name <- EZget(obj,
                          features = features,
                          populations = populations,
                          fraction_design = fraction_design,
                          repeatID = repeatID,
                          exactMatch = exactMatch,
                          returnNameOnly = TRUE)

  if(is.null(fractions_name)){
    stop("No tables from your EZbakR analysis match your search criteria!")
  }

  # Get fractions
  kinetics <- obj[["fractions"]][[fractions_name]]

  features_to_analyze <- obj[["metadata"]][["fractions"]][[fractions_name]][["features"]]


  # Determine which column to use for kinetic parameter estimation
  fraction_cols <- colnames(kinetics)

  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]

  logit_of_interest <- paste0("logit_", fraction_of_interest)
  se_of_interest <- paste0("se_logit_", fraction_of_interest)
  ci_lo_col <- paste0("ci_lower_", fraction_of_interest)
  ci_hi_col <- paste0("ci_upper_", fraction_of_interest)
  has_prof_ci <- all(c(ci_lo_col, ci_hi_col) %in% fraction_cols)
  z95 <- stats::qnorm(0.975)

  extract_fn_ci <- function(df_in){
    if(has_prof_ci){
      lo <- pmin(pmax(as.numeric(df_in[[ci_lo_col]]), 0), 1)
      hi <- pmin(pmax(as.numeric(df_in[[ci_hi_col]]), 0), 1)
    }else{
      fn_v <- as.numeric(df_in[[fraction_of_interest]])
      lfn <- if(logit_of_interest %in% colnames(df_in)){
        as.numeric(df_in[[logit_of_interest]])
      }else{
        logit(fn_v)
      }
      se_v <- as.numeric(df_in[[se_of_interest]])
      se_v <- ifelse(is.finite(se_v), se_v, 1.5)
      lo <- ifelse(!is.finite(lfn) & lfn < 0, 0, inv_logit(lfn - z95 * se_v))
      hi <- ifelse(!is.finite(lfn) & lfn > 0, 1, inv_logit(lfn + z95 * se_v))
    }
    list(lo = lo, hi = hi)
  }

  ### Normalize read counts


  # TO-DO: ALLOW USERS TO JUST USE THE TPM FROM ISOFORM QUANTIFICATION
  reads_norm <- get_normalized_read_counts(obj = obj,
                                           features_to_analyze = features_to_analyze,
                                           fractions_name = fractions_name,
                                           feature_lengths = feature_lengths,
                                           scale_factors = scale_factors)



  if(strategy == "standard"){

    ### Estimate kdegs

    kinetics <- setDT(data.table::copy(kinetics))

    # Add label time info
    metadf <- obj$metadf

    metadf <- setDT(data.table::copy(metadf))
    setkey(metadf, sample)
    setkey(kinetics, sample)

    kinetics <- kinetics[metadf[,c("sample", "tl")], nomatch = NULL]

    # Make sure no -s4U controls made it through
    kinetics <- kinetics[tl > 0]


    kinetics[, kdeg := -log1p(-get(fraction_of_interest))/tl]
    kinetics[, log_kdeg := log(kdeg)]

    fn_ci <- extract_fn_ci(kinetics)
    fn_lo_v <- fn_ci$lo
    fn_hi_v <- fn_ci$hi
    kinetics[, ci_lower_kdeg := ifelse(fn_lo_v <= 0, 0, -log1p(-pmin(fn_lo_v, 1)) / tl)]
    kinetics[, ci_upper_kdeg := ifelse(fn_hi_v >= 1, Inf, -log1p(-pmax(fn_hi_v, 0)) / tl)]
    kinetics[, thalf := log(2) / kdeg]
    kinetics[, ci_lower_thalf := log(2) / ci_upper_kdeg]
    kinetics[, ci_upper_thalf := ifelse(ci_lower_kdeg <= 0, Inf, log(2) / ci_lower_kdeg)]


    ### Estimate uncertainty in log(kdeg)

    lkdeg_uncert <- function(fn, se_lfn){

      se_lfn <- ifelse(is.nan(se_lfn) | !is.finite(se_lfn),
             1.5,
             se_lfn)

      deriv <- (1/log1p(-fn)) * (1 / (1 - fn)) * (fn*(1 - fn))

      uncert <- abs(deriv)*se_lfn

      return(uncert)

    }

    se_delta <- lkdeg_uncert(fn = kinetics[[fraction_of_interest]],
                             se_lfn = kinetics[[se_of_interest]])

    if(has_prof_ci){
      log_ci_lo <- ifelse(fn_lo_v > 0, log(kinetics$ci_lower_kdeg), -Inf)
      log_ci_hi <- ifelse(fn_hi_v < 1, log(kinetics$ci_upper_kdeg), Inf)
      prior_sd_lk <- calibrate_logk_prior_sd(
        sample_vec = kinetics$sample,
        est = kinetics$log_kdeg,
        ci_lo = log_ci_lo,
        ci_hi = log_ci_hi,
        z = z95
      )
      # Concavity of h(eta) = log(-log(1 - inv_logit(eta))) compresses the upper
      # secant log(ci_upper_kdeg) - log(kdeg) logarithmically as fn_hi -> 1;
      # flooring the upper scale with the tangent projection h'(eta_hat) * (eta_hi - eta_hat)
      # agrees to second order in the interior and diverges smoothly as fn_hi -> 1.
      n_safe <- pmax(kinetics$n, 1)
      se_cr_lo <- 1 / sqrt(0.64761 * n_safe)
      fn_v <- pmin(pmax(as.numeric(kinetics[[fraction_of_interest]]), 1e-15), 1 - 1e-15)
      deriv_abs <- abs(fn_v / log1p(-fn_v))
      se_tan_hi <- ifelse(
        fn_hi_v < 1,
        pmax(se_cr_lo, deriv_abs * pmax(stats::qlogis(fn_hi_v) - stats::qlogis(fn_v), 0) / z95),
        Inf
      )

      se_delta <- posterior_se_from_ci(
        est = kinetics$log_kdeg,
        ci_lo = log_ci_lo,
        ci_hi = log_ci_hi,
        prior_sd = prior_sd_lk,
        se_min_lo = se_cr_lo,
        se_min_hi = se_tan_hi,
        z = z95
      )
    }

    kinetics[, se_log_kdeg := se_delta]


    ### Estimate ksyn

    # Merge with kinetics
    setkeyv(reads_norm, c("sample", features_to_analyze))
    setkeyv(kinetics, c("sample", features_to_analyze))

    cols_to_keep <- c("sample", features_to_analyze, "normalized_reads", "scale_factor")
    kinetics_cols_to_keep <- c(cols_to_keep, "kdeg", "log_kdeg", "se_log_kdeg",
                               "ci_lower_kdeg", "ci_upper_kdeg",
                               "thalf", "ci_lower_thalf", "ci_upper_thalf", "n")

    kinetics <- kinetics[reads_norm[,..cols_to_keep], ..kinetics_cols_to_keep, nomatch = NULL]

    # Estimate ksyn
    kinetics[, ksyn := normalized_reads*kdeg]
    kinetics[, log_ksyn := log(ksyn)]

    # Estimate uncertainty (assuming raw count n ~ Poisson(n))
    kinetics[, se_log_ksyn := sqrt( (1/n) + se_log_kdeg^2)]

    # MOVER 95% CIs for ksyn = normalized_reads * kdeg
    hw_r_std <- z95 / sqrt(pmax(kinetics$n, 1))
    kinetics[, ci_lower_ksyn := ifelse(
      ci_lower_kdeg <= 0,
      0,
      ksyn * exp(-sqrt(pmax(log_kdeg - log(ci_lower_kdeg), 0)^2 + hw_r_std^2))
    )]
    kinetics[, ci_upper_ksyn := ifelse(
      !is.finite(ci_upper_kdeg),
      Inf,
      ksyn * exp(sqrt(pmax(log(ci_upper_kdeg) - log_kdeg, 0)^2 + hw_r_std^2))
    )]



  }else if(strategy == "NSS"){

    ### IDEA
    # After working through the mathematical formalism, turns out
    # the ideas I had are identical to those presented inNarain et al., 2021.
    # Best you can do is fit a single exponential model using -s4U and +s4U data
    # to infer the "fraction of old reads that are still present by the end of
    # the label time". Turns out to be identical to estimating the average rate
    # constant for a non-homogeneous Poisson process during the label time. If
    # decay dynamics are non-Poissonian (time to decay is not ~ Exponential()),
    # then memorylessness is out the window and you are always going to be screwed.
    # I get the sense that exponential departure time is something of an approximation
    # or central limit theorem of almost any "standard" departure process. Pretty
    # much will be a good approximation as long as there is a predominantly rate-limiting
    # step in the RNA decay process.
    #
    # Analysis steps:
    # 1) Extract -s4U read counts from relevant source
    # Default is just from fractions object
    # Can also use a read count data frame stored in EZbakRData object
    # 2) Normalize read counts
    # Default is TMM
    # Can also provide a table of scale factors
    # Maybe can also choose not to normalize (technically would be same
    # as providing scale table of 1, but I can autmoate that for user)
    # 3) Estimate kdeg = -log(norm. old reads +s4U / norm. -s4U reads)/tl
    # 4) Estimate ksyn = (fn * norm. total reads +s4U * kdeg)/(1 - exp(-kdeg*tl))

    ### NOTE: a lot of steps 1 and 2 are duplicated from CorrectDropout, so
    ### really should refactor this in the future


    ##### STEP 1: GET -S4U READ COUNTS

    ### Which columns should -s4U samples be grouped by?

    metadf <- obj$metadf


    if(is.null(grouping_factor) | is.null(reference_factor)){

      stop("Grouping factors and reference factors must both be specified if
           using the NSS strategy!")

    }else if(length(grouping_factor) != 1 | length(reference_factor) != 1){

      stop("grouping_factor and reference_factor must each be a single
           column of your metadf!")

    }else if(!(grouping_factor %in% colnames(metadf)) | !(reference_factor %in% colnames(metadf))){

      stop("grouping_factor and reference_factor must both share a name
           with columns in your metadf!")

    }

    # Necessary generalizations:
    # 1) Metadf column used (e.g., pulse-chase)
    reference_data <- kinetics %>%
      dplyr::inner_join(metadf,
                        by = "sample") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(reference_rpm = n/(sum(n)/1000000),
                    reference_n = n) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(grouping_factor, features_to_analyze)))) %>%
      dplyr::summarise(reference_rpm = mean(reference_rpm),
                       reference_n = mean(reference_n),
                       reference_reps = dplyr::n()) %>%
      dplyr::select(!!grouping_factor, !!features_to_analyze,
                    reference_rpm, reference_n, reference_reps)


    ##### STEP 2: INTEGRATE WITH +S4U TO GET ADJUSTED FRACTION NEW
    kinetics <- kinetics %>%
      dplyr::inner_join(metadf %>% dplyr::filter(tl > 0) %>%
                          dplyr::select(sample, tl,
                                        !!reference_factor),
                        by = "sample") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(rpm = n/(sum(n)/1000000)) %>%
      dplyr::inner_join(reference_data %>%
                          dplyr::rename(!!reference_factor := !!dplyr::sym(grouping_factor)),
                        by = c(reference_factor, features_to_analyze)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(old_rpm = (1 - !!dplyr::sym(fraction_of_interest)) * rpm,
                    new_rpm = (!!dplyr::sym(fraction_of_interest)) * rpm) %>%
      dplyr::mutate(kdeg = dplyr::case_when(
        old_rpm >= reference_rpm ~ -log(1 - !!dplyr::sym(fraction_of_interest))/tl,
        .default = -log(old_rpm/reference_rpm)/tl
      ),
      ksyn = new_rpm * kdeg / (1 - exp(-kdeg*tl)),
      log_kdeg = log(kdeg),
      log_ksyn = log(ksyn))


    # Add normalized read counts
    kinetics <- kinetics %>%
      dplyr::inner_join(reads_norm %>%
                          dplyr::select(sample, !!features_to_analyze, "normalized_reads"),
                        by = c("sample", features_to_analyze))


    ##### STEP 3: ESTIMATE UNCERTAINTY
    ##### (higher than in standard due to -s4U read variance)

    ### TO-DO: Incorporate number of -s4U replicates used to infer nolabel_n

    ### Uncertainty approximated with delta method
    lkdeg_uncert_nss <- function(fn, se_lfn, Rs4U, Rctl, tl, kdeg, reference_reps = 1){

      se_lfn <- ifelse(is.nan(se_lfn), 1.5, se_lfn)

      part_1 <- (fn^2)*(se_lfn^2)/(tl^2)
      part_2 <- ((1/(Rs4U + 1))^2)*Rs4U/(tl^2)
      part_3 <- ((1/(Rctl + 1))^2)*Rctl/((tl^2)*reference_reps)
      kdeg_var <- part_1 + part_2 + part_3
      uncert <- sqrt(((1/kdeg)^2)*kdeg_var)

      return(uncert)

    }

    lksyn_uncert_nss <- function(fn, se_lfn, Rs4U, Rctl, tl, kdeg, reference_reps = 1){

      se_lfn <- ifelse(is.nan(se_lfn), 1.5, se_lfn)

      # Differentiate log(ksyn) = log(fn) + log(Rs4U) - log(tl) + log(u/(1 - exp(-u))) + const
      # directly wrt independent (logit(fn), log(Rs4U), log(Rctl)), where u = tl * kdeg.
      u <- tl * kdeg
      h_prime <- ifelse(u < 1e-5, 0.5 - u / 12, (1 / u) - exp(-u) / (1 - exp(-u)))

      part_1s <- ((1 - fn + fn * h_prime)^2) * (se_lfn^2)
      part_2s <- ((1 - h_prime)^2) * ((1/(Rs4U + 1))^2) * Rs4U
      part_3s <- (h_prime^2) * ((1/(Rctl + 1))^2) * Rctl / reference_reps
      lksyn_var <- part_1s + part_2s + part_3s

      return(sqrt(lksyn_var))


    }

    kinetics <- setDT(kinetics)
    kinetics[, se_log_kdeg := lkdeg_uncert_nss(fn = get(fraction_of_interest),
                                               se_lfn = get(se_of_interest),
                                               Rs4U = n,
                                               Rctl = reference_n,
                                               tl = tl,
                                               kdeg = kdeg,
                                               reference_reps = reference_reps)]

    # Estimate uncertainty (assuming normalized_reads ~ Poisson(normalized_reads)/scale_factor)
    kinetics[, se_log_ksyn := lksyn_uncert_nss(fn = get(fraction_of_interest),
                                               se_lfn = get(se_of_interest),
                                               Rs4U = n,
                                               Rctl = reference_n,
                                               tl = tl,
                                               kdeg = kdeg,
                                               reference_reps = reference_reps)]

    fn_ci_nss <- extract_fn_ci(kinetics)
    kinetics[, ci_lower_kdeg := ifelse(
      fn_ci_nss$lo <= 0 & old_rpm >= reference_rpm,
      0,
      exp(log_kdeg - z95 * se_log_kdeg)
    )]
    kinetics[, ci_upper_kdeg := ifelse(
      fn_ci_nss$hi >= 1,
      Inf,
      exp(log_kdeg + z95 * se_log_kdeg)
    )]
    kinetics[, thalf := log(2) / kdeg]
    kinetics[, ci_lower_thalf := log(2) / ci_upper_kdeg]
    kinetics[, ci_upper_thalf := ifelse(ci_lower_kdeg <= 0, Inf, log(2) / ci_lower_kdeg)]
    kinetics[, ci_lower_ksyn := ifelse(fn_ci_nss$lo <= 0, 0, exp(log_ksyn - z95 * se_log_ksyn))]
    kinetics[, ci_upper_ksyn := ifelse(fn_ci_nss$hi >= 1, Inf, exp(log_ksyn + z95 * se_log_ksyn))]



  }else if(strategy == "shortfeed"){

    ### IDEA
    # Take it to the limit where there is no degradation (so dR/dt = ks).
    # kdeg = fn / tl
    # ksyn = fn * (normalized read count) / tl; basically new RNA produced per unit time
    # Unclear how much this will matter, as it is literally just the limit of
    # the standard analysis, but worth implementing I think.


    ### NOTE: Currently a near perfect duplication of standard code, just
    ### refactor into a function for goodness sake!

    ### Estimate kdegs

    kinetics <- setDT(data.table::copy(kinetics))

    # Add label time info
    metadf <- obj$metadf

    metadf <- setDT(data.table::copy(metadf))
    setkey(metadf, sample)
    setkey(kinetics, sample)

    kinetics <- kinetics[metadf[,c("sample", "tl")], nomatch = NULL]

    # Make sure no -s4U controls made it through
    kinetics <- kinetics[tl > 0]


    kinetics[, kdeg := get(fraction_of_interest)/tl]
    kinetics[, log_kdeg := log(kdeg)]

    fn_ci_sf <- extract_fn_ci(kinetics)
    kinetics[, ci_lower_kdeg := fn_ci_sf$lo / tl]
    kinetics[, ci_upper_kdeg := fn_ci_sf$hi / tl]
    kinetics[, thalf := log(2) / kdeg]
    kinetics[, ci_lower_thalf := ifelse(ci_upper_kdeg <= 0, Inf, log(2) / ci_upper_kdeg)]
    kinetics[, ci_upper_thalf := ifelse(ci_lower_kdeg <= 0, Inf, log(2) / ci_lower_kdeg)]




    ### Estimate uncertainty in log(kdeg)

    # This needs to be a bit different in this case
    lkdeg_uncert_sf <- function(fn, se_lfn){

      deriv <- (1/fn)*(fn*(1-fn))

      uncert <- abs(deriv)*se_lfn

      return(uncert)

    }


    se_sf <- lkdeg_uncert_sf(fn = kinetics[[fraction_of_interest]],
                             se_lfn = kinetics[[se_of_interest]])
    if(has_prof_ci){
      log_ci_lo_sf <- ifelse(fn_ci_sf$lo > 0, log(kinetics$ci_lower_kdeg), -Inf)
      log_ci_hi_sf <- ifelse(fn_ci_sf$hi < 1 & kinetics$ci_upper_kdeg > 0, log(kinetics$ci_upper_kdeg), Inf)
      prior_sd_sf <- calibrate_logk_prior_sd(
        sample_vec = kinetics$sample,
        est = kinetics$log_kdeg,
        ci_lo = log_ci_lo_sf,
        ci_hi = log_ci_hi_sf,
        z = z95
      )
      se_sf <- posterior_se_from_ci(
        est = kinetics$log_kdeg,
        ci_lo = log_ci_lo_sf,
        ci_hi = log_ci_hi_sf,
        prior_sd = prior_sd_sf,
        se_min_lo = 0.5 / sqrt(pmax(kinetics$n, 1)),
        se_min_hi = ifelse(fn_ci_sf$hi < 1,
                           0.5 / sqrt(pmax(kinetics$n * (1 - fn_ci_sf$hi), 1e-30)),
                           Inf),
        z = z95
      )
    }
    kinetics[, se_log_kdeg := se_sf]


    ### Estimate ksyn

    # Merge with kinetics
    setkeyv(reads_norm, c("sample", features_to_analyze))
    setkeyv(kinetics, c("sample", features_to_analyze))

    cols_to_keep <- c("sample", features_to_analyze, "normalized_reads", "scale_factor")
    kinetics_cols_to_keep <- c(cols_to_keep, "kdeg", "log_kdeg", "se_log_kdeg",
                               "ci_lower_kdeg", "ci_upper_kdeg",
                               "thalf", "ci_lower_thalf", "ci_upper_thalf", "n")

    kinetics <- kinetics[reads_norm[,..cols_to_keep], ..kinetics_cols_to_keep, nomatch = NULL]

    # Estimate ksyn
    kinetics[, ksyn := normalized_reads*kdeg]
    kinetics[, log_ksyn := log(ksyn)]

    # Estimate uncertainty (assuming raw count n ~ Poisson(n))
    kinetics[, se_log_ksyn := sqrt( (1/n) + se_log_kdeg^2)]

    hw_r_sf <- z95 / sqrt(pmax(kinetics$n, 1))
    kinetics[, ci_lower_ksyn := ifelse(
      ci_lower_kdeg <= 0,
      0,
      ksyn * exp(-sqrt(pmax(log_kdeg - log(ci_lower_kdeg), 0)^2 + hw_r_sf^2))
    )]
    kinetics[, ci_upper_ksyn := ksyn * exp(sqrt(pmax(log(ci_upper_kdeg) - log_kdeg, 0)^2 + hw_r_sf^2))]


  }else if(strategy == "pulse-chase"){

    # Add label time info
    metadf <- obj$metadf


    group_cols <- c(features_to_analyze, grouping_factor)


    if(!all(c("tpulse", "tchase") %in% colnames(metadf))){

      stop("If strategy == pulse-chase, metadf needs tpulse and tchase columns!")

    }

    cols_to_keep <- c("sample", features_to_analyze, "normalized_reads", "scale_factor")
    kinetics_cols_to_keep <- c(cols_to_keep, "kdeg", "log_kdeg", "se_log_kdeg",
                               "fn_lo_tmp", "fn_hi_tmp", "n")

    fn_ci_pc <- extract_fn_ci(kinetics)
    kinetics$fn_lo_tmp <- fn_ci_pc$lo
    kinetics$fn_hi_tmp <- fn_ci_pc$hi

    kinetics <- kinetics %>%
      dplyr::inner_join(
        metadf %>%
          dplyr::select(
            sample, tpulse, tchase, !!grouping_factor
          ),
        by = "sample"
      ) %>%
      dplyr::filter(
        tpulse > 0
      ) %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(group_cols)
        )
      ) %>%
      dplyr::mutate(
        fn_safe = dplyr::case_when(
          !!dplyr::sym(fraction_of_interest) <= 0 ~ 1 / (n + 2),
          !!dplyr::sym(fraction_of_interest) >= 1 ~ 1 - 1 / (n + 2),
          .default = !!dplyr::sym(fraction_of_interest)
        ),
        se_safe = ifelse(is.finite(!!dplyr::sym(se_of_interest)),
                         !!dplyr::sym(se_of_interest),
                         1.5),
        kdeg = dplyr::case_when(
          tchase == 0 ~ -log(1 - fn_safe) / tpulse,
          .default = -log(fn_safe / mean(fn_safe[tchase == 0])) / tchase
        ),
        kdeg = ifelse(!is.finite(kdeg) | kdeg <= 0,
                      -log(1 - 1/(n + 2)) / ifelse(tchase == 0, tpulse, tchase), # theoretical dynamic range of experiment
                      kdeg),
        log_kdeg = log(kdeg),
        se_log_kdeg = dplyr::case_when(
          tchase == 0 ~ se_safe * abs(fn_safe / log(1 - fn_safe)),
          .default = sqrt(
            sum(
              (fn_safe[tchase == 0] *
                 (1 - fn_safe[tchase == 0]) *
                 se_safe[tchase == 0])^2
            ) / (sum(fn_safe[tchase == 0])^2) +
              ((1 - fn_safe) * se_safe)^2
          ) / (tchase * kdeg)
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::inner_join(
        reads_norm %>%
          dplyr::select(!!cols_to_keep),
        by = c("sample", features_to_analyze)
      ) %>%
      dplyr::select(!!c(kinetics_cols_to_keep, "tpulse", "tchase")) %>%
      dplyr::mutate(
        ci_lower_kdeg = ifelse(fn_lo_tmp <= 0 & tchase == 0, 0, exp(log_kdeg - z95 * se_log_kdeg)),
        ci_upper_kdeg = ifelse((fn_hi_tmp >= 1 & tchase == 0) | (fn_lo_tmp <= 0 & tchase > 0),
                               Inf, exp(log_kdeg + z95 * se_log_kdeg)),
        thalf = log(2) / kdeg,
        ci_lower_thalf = log(2) / ci_upper_kdeg,
        ci_upper_thalf = ifelse(ci_lower_kdeg <= 0, Inf, log(2) / ci_lower_kdeg),
        ksyn = normalized_reads * kdeg,
        log_ksyn = log(ksyn),
        se_log_ksyn = sqrt( (1/n) + se_log_kdeg^2),
        ci_lower_ksyn = ifelse(ci_lower_kdeg <= 0, 0, exp(log_ksyn - z95 * se_log_ksyn)),
        ci_upper_ksyn = ifelse(!is.finite(ci_upper_kdeg), Inf, exp(log_ksyn + z95 * se_log_ksyn))
      ) %>%
      dplyr::select(-fn_lo_tmp, -fn_hi_tmp)


    if(exclude_pulse_estimates){

      kinetics <- kinetics %>%
        dplyr::filter(tchase > 0)

    }


    # Remove metadf columns
    kinetics <- kinetics %>%
      dplyr::select(
        -tpulse, -tchase
      )


  }



  ##### PROCESS OUTPUT AND RETURN

  # What should output be named?
  kinetics_vect <- paste(gsub("_","",features_to_analyze), collapse = "_")

  if(nchar(kinetics_vect) > character_limit){

    num_kinetics <- length(obj[['kinetics']])
    kinetics_vect <- paste0("kinetics", num_kinetics + 1)

  }


  readcount_vect <- paste0(paste(gsub("_", "", features_to_analyze), collapse = "_"),
                           "_readcount_df")

  if(nchar(readcount_vect) > character_limit){

    num_reads <- sum(grepl("^readcount_df", names(obj[['readcounts']])))
    readcount_vect <- paste0("readcount_df", num_reads + 1)

  }


  # Are there any metadata objects
  if(length(obj[['metadata']][['kinetics']]) > 0){

    kinetics_vect <- decide_output(obj,
                                   proposed_name = kinetics_vect,
                                   type = "kinetics",
                                   features = features_to_analyze,
                                   kstrat = strategy,
                                   overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      repeatID <- 1

    }else{

      repeatID <- length(EZget(obj,
                               type = 'kinetics',
                               features = features_to_analyze,
                               kstrat = strategy,
                               returnNameOnly = TRUE,
                               exactMatch = TRUE,
                               alwaysCheck = TRUE)) + 1
    }

  }else{

    repeatID <- 1

  }


  if(length(obj[['metadata']][['readcounts']]) > 0){

    readcount_vect <- decide_output(obj,
                                    proposed_name = readcount_vect,
                                    type = "readcounts",
                                    features = features_to_analyze,
                                    counttype = "TMM_normalized",
                                    overwrite = overwrite)

    # How many identical tables already exist?
    if(overwrite){

      rc_repeatID <- 1

    }else{

      rc_repeatID <- length(EZget(obj,
                                  type = 'readcounts',
                                  features = features_to_analyze,
                                  counttype = "TMM_normalized",
                                  returnNameOnly = TRUE,
                                  exactMatch = TRUE,
                                  alwaysCheck = TRUE)) + 1
    }

  }else{

    rc_repeatID <- 1

  }

  obj[["kinetics"]][[kinetics_vect]] <- dplyr::as_tibble(kinetics) %>%
    dplyr::select(sample, !!features_to_analyze,
                  kdeg, log_kdeg, se_log_kdeg, ci_lower_kdeg, ci_upper_kdeg,
                  thalf, ci_lower_thalf, ci_upper_thalf,
                  ksyn, log_ksyn, se_log_ksyn, ci_lower_ksyn, ci_upper_ksyn,
                  normalized_reads, n)

  # Eventually want to add count matrix output
  obj[["readcounts"]][[readcount_vect]] <- dplyr::as_tibble(reads_norm) %>%
    dplyr::select(sample, !!features_to_analyze, n, normalized_reads, dplyr::any_of("geom_mean"), scale_factor)


  obj[["metadata"]][["kinetics"]][[kinetics_vect]] <- list(features = features_to_analyze,
                                                           kstrat = strategy,
                                                           repeatID = repeatID)
  obj[["metadata"]][["readcounts"]][[readcount_vect]] <- list(features = features_to_analyze,
                                                              counttype = "TMM_normalized",
                                                              repeatID = rc_repeatID)


  if(!methods::is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)


}


# Estimate TILAC ratio
tilac_ratio_estimation <- function(obj,
                                   features = NULL,
                                   populations = NULL,
                                   exactMatch = TRUE,
                                   fraction_design = NULL,
                                   grouping_factor = NULL,
                                   repeatID = NULL,
                                   character_limit = 20,
                                   feature_lengths = NULL,
                                   scale_factors = NULL,
                                   overwrite = TRUE){

  ### Hack to deal with devtools::check() NOTEs
  ..cols_to_keep <- TILAC_ratio <- log_TILAC_ratio <- se_log_TILAC_ratio <- ..cols_keep <- NULL
  ci_lower_TILAC_ratio <- ci_upper_TILAC_ratio <- NULL
  n <- normalized_reads <- scale_factor <- NULL

  rm(..cols_to_keep)
  rm(..cols_keep)

  `.` <- list

  ### Figure out which fraction new estimates to use

  fractions_name <- EZget(obj,
                          type = "fractions",
                          features = features,
                          populations = populations,
                          fraction_design = fraction_design,
                          repeatID = repeatID,
                          exactMatch = exactMatch,
                          returnNameOnly = TRUE)

  if(is.null(fractions_name)){
    stop("No tables from your EZbakR analysis match your search criteria!")
  }

  # Get fractions
  kinetics <- obj[["fractions"]][[fractions_name]]
  features_to_analyze <- obj[["metadata"]][["fractions"]][[fractions_name]][["features"]]


  ### Estimate TILAC ratio

  # Determine which columns represent the two single-label-high TILAC populations
  fraction_cols <- colnames(kinetics)
  all_fraction_cols <- fraction_cols[grepl("^fraction_", fraction_cols)]
  nhigh <- lengths(regmatches(all_fraction_cols, gregexpr("high", all_fraction_cols)))
  fraction_of_interest <- all_fraction_cols[nhigh == 1]

  # TODO: COULD ADD PARAMETER TO ALLOW USER TO DECIDE WHICH FRACTION COLUMNS TO USE
  if(length(fraction_of_interest) != 2){

    stop("There are not exactly two single-`high` mutational content fractions!
         Are you sure this is TILAC data?")

  }

  kinetics <- setDT(data.table::copy(kinetics))

  # Add label time info
  metadf <- obj$metadf
  metacols <- colnames(metadf)

  metadf <- setDT(data.table::copy(metadf))
  setkey(metadf, sample)
  setkey(kinetics, sample)

  tl_cols <- metacols[grepl("^tl($|_)", metacols)]


  cols_keep <- c("sample", tl_cols)
  kinetics <- kinetics[metadf[,..cols_keep], nomatch = NULL]

  # Make sure no -label controls made it through
  kinetics <- setDT(dplyr::as_tibble(kinetics) %>%
                      dplyr::filter(dplyr::if_all(dplyr::all_of(tl_cols), ~ . > 0)))


  kinetics[, TILAC_ratio := get(fraction_of_interest[1])/get(fraction_of_interest[2])]
  kinetics[, log_TILAC_ratio := log(TILAC_ratio)]

  # Estimate uncertainty in log_TILAC_ratio via exact CLR parallelogram identity
  # Var(z_1 - z_2) = 2*(Var(z_1) + Var(z_2)) - Var(z_3) when K <= 3, with
  # Bayesian multinomial fallback if se_invsoftmax_* columns are absent.
  pop_suffixes <- sub("^fraction_", "", all_fraction_cols)
  foi_suffixes <- sub("^fraction_", "", fraction_of_interest)
  se_ism_cols <- paste0("se_invsoftmax_", pop_suffixes)
  se_foi_cols <- paste0("se_invsoftmax_", foi_suffixes)

  fallback_var <- 1 / (kinetics[[fraction_of_interest[1]]] * kinetics[["n"]] + 1) +
    1 / (kinetics[[fraction_of_interest[2]]] * kinetics[["n"]] + 1)

  if(all(se_ism_cols %in% colnames(kinetics)) && length(all_fraction_cols) <= 3){
    other_se_cols <- setdiff(se_ism_cols, se_foi_cols)
    if(length(other_se_cols) == 1){
      clr_var <- 2 * (kinetics[[se_foi_cols[1]]]^2 + kinetics[[se_foi_cols[2]]]^2) -
        kinetics[[other_se_cols[1]]]^2
    }else{
      clr_var <- 2 * (kinetics[[se_foi_cols[1]]]^2 + kinetics[[se_foi_cols[2]]]^2)
    }
    var_log_ratio <- ifelse(!is.finite(clr_var) | clr_var <= 0, fallback_var, clr_var)
  }else{
    var_log_ratio <- fallback_var
  }

  kinetics[, se_log_TILAC_ratio := sqrt(var_log_ratio)]

  z95_tilac <- stats::qnorm(0.975)
  ci_lo_foi <- paste0("ci_lower_", foi_suffixes)
  ci_hi_foi <- paste0("ci_upper_", foi_suffixes)
  # Check both ci_lower_fraction_<pop> and ci_lower_<foi>
  ci_lo_cols_t <- paste0("ci_lower_", fraction_of_interest)
  ci_hi_cols_t <- paste0("ci_upper_", fraction_of_interest)
  if(all(c(ci_lo_cols_t, ci_hi_cols_t) %in% colnames(kinetics))){
    f1 <- as.numeric(kinetics[[fraction_of_interest[1L]]])
    f2 <- as.numeric(kinetics[[fraction_of_interest[2L]]])
    f1_lo <- as.numeric(kinetics[[ci_lo_cols_t[1L]]])
    f1_hi <- as.numeric(kinetics[[ci_hi_cols_t[1L]]])
    f2_lo <- as.numeric(kinetics[[ci_lo_cols_t[2L]]])
    f2_hi <- as.numeric(kinetics[[ci_hi_cols_t[2L]]])

    hw1_lo <- ifelse(f1_lo > 0, pmax(log(f1 / f1_lo), 0), Inf)
    hw1_hi <- ifelse(f1_hi > 0, pmax(log(f1_hi / f1), 0), 0)
    hw2_lo <- ifelse(f2_lo > 0, pmax(log(f2 / f2_lo), 0), Inf)
    hw2_hi <- ifelse(f2_hi > 0, pmax(log(f2_hi / f2), 0), 0)

    hw_tot_lo <- pmax(z95_tilac * kinetics$se_log_TILAC_ratio, sqrt(hw1_lo^2 + hw2_hi^2))
    hw_tot_hi <- pmax(z95_tilac * kinetics$se_log_TILAC_ratio, sqrt(hw1_hi^2 + hw2_lo^2))

    kinetics[, ci_lower_TILAC_ratio := ifelse(!is.finite(hw_tot_lo), 0, exp(log_TILAC_ratio - hw_tot_lo))]
    kinetics[, ci_upper_TILAC_ratio := ifelse(!is.finite(hw_tot_hi), Inf, exp(log_TILAC_ratio + hw_tot_hi))]

    log_ci_lo_t <- ifelse(is.finite(hw_tot_lo), kinetics$log_TILAC_ratio - hw_tot_lo, -Inf)
    log_ci_hi_t <- ifelse(is.finite(hw_tot_hi), kinetics$log_TILAC_ratio + hw_tot_hi, Inf)
    prior_sd_tilac <- calibrate_logk_prior_sd(
      sample_vec = kinetics$sample,
      est = kinetics$log_TILAC_ratio,
      ci_lo = log_ci_lo_t,
      ci_hi = log_ci_hi_t,
      z = z95_tilac
    )
    se_tilac_fl <- posterior_se_from_ci(
      est = kinetics$log_TILAC_ratio,
      ci_lo = log_ci_lo_t,
      ci_hi = log_ci_hi_t,
      prior_sd = prior_sd_tilac,
      se_min = 0.5 / sqrt(pmax(kinetics$n, 1)),
      z = z95_tilac
    )
    kinetics[, se_log_TILAC_ratio := se_tilac_fl]
  }else{
    kinetics[, ci_lower_TILAC_ratio := exp(log_TILAC_ratio - z95_tilac * se_log_TILAC_ratio)]
    kinetics[, ci_upper_TILAC_ratio := exp(log_TILAC_ratio + z95_tilac * se_log_TILAC_ratio)]
  }


  ### Get normalized read counts

  reads_norm <- get_normalized_read_counts(obj = obj,
                                           features_to_analyze = features_to_analyze,
                                           fractions_name = fractions_name,
                                           feature_lengths = feature_lengths,
                                           scale_factors = scale_factors)


  # Merge normalized_reads with kinetics
  setkeyv(reads_norm, c("sample", features_to_analyze))
  setkeyv(kinetics, c("sample", features_to_analyze))

  cols_to_keep <- c("sample", features_to_analyze, "normalized_reads")

  kinetics <- kinetics[reads_norm[,..cols_to_keep], nomatch = NULL]


  # Naming & metadata identical to Standard_kinetic_estimation
  kinetics_vect <- paste(gsub("_", "", features_to_analyze), collapse = "_")

  if(nchar(kinetics_vect) > character_limit){
    num_kinetics <- length(obj[['kinetics']])
    kinetics_vect <- paste0("kinetics", num_kinetics + 1)
  }

  readcount_vect <- paste0(paste(gsub("_", "", features_to_analyze), collapse = "_"),
                           "_readcount_df")

  if(nchar(readcount_vect) > character_limit){
    num_reads <- sum(grepl("^readcount_df", names(obj[['readcounts']])))
    readcount_vect <- paste0("readcount_df", num_reads + 1)
  }

  if(length(obj[['metadata']][['kinetics']]) > 0){
    kinetics_vect <- decide_output(obj,
                                   proposed_name = kinetics_vect,
                                   type = "kinetics",
                                   features = features_to_analyze,
                                   kstrat = "tilac",
                                   overwrite = overwrite)
    if(overwrite){
      repeatID_out <- 1
    }else{
      repeatID_out <- length(EZget(obj,
                                   type = 'kinetics',
                                   features = features_to_analyze,
                                   kstrat = "tilac",
                                   returnNameOnly = TRUE,
                                   exactMatch = TRUE,
                                   alwaysCheck = TRUE)) + 1
    }
  }else{
    repeatID_out <- 1
  }

  if(length(obj[['metadata']][['readcounts']]) > 0){
    readcount_vect <- decide_output(obj,
                                    proposed_name = readcount_vect,
                                    type = "readcounts",
                                    features = features_to_analyze,
                                    counttype = "TMM_normalized",
                                    overwrite = overwrite)
    if(overwrite){
      rc_repeatID <- 1
    }else{
      rc_repeatID <- length(EZget(obj,
                                  type = 'readcounts',
                                  features = features_to_analyze,
                                  counttype = "TMM_normalized",
                                  returnNameOnly = TRUE,
                                  exactMatch = TRUE,
                                  alwaysCheck = TRUE)) + 1
    }
  }else{
    rc_repeatID <- 1
  }

  obj[['kinetics']][[kinetics_vect]] <- dplyr::as_tibble(kinetics) %>%
    dplyr::select(sample, !!features_to_analyze, TILAC_ratio, log_TILAC_ratio,
                  se_log_TILAC_ratio, ci_lower_TILAC_ratio, ci_upper_TILAC_ratio,
                  normalized_reads, n)

  obj[['readcounts']][[readcount_vect]] <- dplyr::as_tibble(reads_norm) %>%
    dplyr::select(sample, !!features_to_analyze, n, normalized_reads, dplyr::any_of("geom_mean"), scale_factor)

  obj[['metadata']][['kinetics']][[kinetics_vect]] <- list(features = features_to_analyze,
                                                           kstrat = "tilac",
                                                           repeatID = repeatID_out)
  obj[['metadata']][['readcounts']][[readcount_vect]] <- list(features = features_to_analyze,
                                                              counttype = "TMM_normalized",
                                                              repeatID = rc_repeatID)

  if(!methods::is(obj, "EZbakRKinetics")){

    class(obj) <- c( "EZbakRKinetics", class(obj))

  }

  return(obj)



}
