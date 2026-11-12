# Helper to construct clean EZbakRFractions object with realistic finite point estimates & profile CIs
make_gof_fractions_obj <- function(rows_df, metadf){
  z95 <- stats::qnorm(0.975)
  fn_raw <- pmin(pmax(rows_df$fraction_highTC, 0), 1)
  n_r <- pmax(rows_df$n, 1)

  # Match EstimateFractions: finite [-9, 9] logit point estimates + exact boundary profile CIs [0, U] / [L, 1]
  a <- fn_raw * n_r + 1
  b <- (1 - fn_raw) * n_r + 1
  f_post <- a / (a + b)
  se_logit <- sqrt(a * b / (((a + b)^2) * (a + b + 1))) / pmax(f_post * (1 - f_post), 1e-6)

  logit_fn <- ifelse(fn_raw <= 0, pmax(stats::qlogis(f_post), -9),
                     ifelse(fn_raw >= 1, pmin(stats::qlogis(f_post), 9),
                            pmin(pmax(stats::qlogis(fn_raw), -9), 9)))
  fn_est <- stats::plogis(logit_fn)

  ci_lo <- ifelse(fn_raw <= 0, 0,
                  ifelse(fn_raw >= 1, 0.05^(1 / n_r),
                         stats::plogis(stats::qlogis(pmin(pmax(fn_raw, 1e-6), 1 - 1e-6)) - z95 * se_logit)))
  ci_hi <- ifelse(fn_raw >= 1, 1,
                  ifelse(fn_raw <= 0, 1 - 0.05^(1 / n_r),
                         stats::plogis(stats::qlogis(pmin(pmax(fn_raw, 1e-6), 1 - 1e-6)) + z95 * se_logit)))

  rows_df$fraction_highTC <- fn_est
  rows_df$logit_fraction_highTC <- logit_fn
  rows_df$se_logit_fraction_highTC <- se_logit
  rows_df$ci_lower_fraction_highTC <- ci_lo
  rows_df$ci_upper_fraction_highTC <- ci_hi

  EZbakRFractions(rows_df, metadf)
}


test_that("AssessGoF is calibrated under multi-tl steady-state null, handles 0/1 boundaries, and disambiguates log_ksyn", {
  set.seed(101)
  tls <- c(1, 2, 4, 8)
  nreps <- 2
  metadf <- dplyr::tibble(
    sample = c(paste0("s_tl", rep(tls, each = nreps), "_r", rep(seq_len(nreps), times = length(tls))), "ctl_1"),
    tl = c(rep(tls, each = nreps), 0),
    condition = "WT"
  )

  ngenes <- 36
  genes <- paste0("Gene", seq_len(ngenes))
  # Span slow (boundary ci_lower == 0 at tl=1) to fast (boundary ci_upper == 1 at tl=8)
  true_kdeg <- exp(seq(log(0.012), log(0.65), length.out = ngenes))

  cB_list <- vector("list", nrow(metadf))
  for(s in seq_len(nrow(metadf))){
    tl_s <- metadf$tl[s]
    fn_s <- if(tl_s == 0) rep(0, ngenes) else (1 - exp(-true_kdeg * tl_s))
    sim_s <- SimulateOneRep(
      nfeatures = ngenes,
      read_vect = rep(250, ngenes),
      fn_vect = fn_s,
      sample_name = metadf$sample[s],
      pnew = 0.08,
      pold = 0.002
    )
    cB_list[[s]] <- sim_s$cB
  }
  cB <- dplyr::bind_rows(cB_list)

  ezbdo <- EZbakRData(cB, metadf)
  ezbdo <- suppressMessages(EstimateFractions(ezbdo, pold_from_nolabel = TRUE))
  ezbdo <- suppressMessages(EstimateKinetics(ezbdo))
  # Use inverse-variance weighted averaging across multi-tl samples so saturated/zero-mutation
  # timepoints with large se_log_kdeg do not bias the kinetic trajectory
  ezbdo <- suppressMessages(AverageAndRegularize(ezbdo, parameter = "log_kdeg", formula_mean = ~ 1, force_optim = TRUE))
  ezbdo <- suppressMessages(AverageAndRegularize(ezbdo, parameter = "log_ksyn", formula_mean = ~ 1))

  # Verify boundary profile intervals (ci_lower == 0 and/or ci_upper == 1) actually occurred in EstimateFractions
  fxn_tbl <- EZget(ezbdo, type = "fractions")
  fxn_lab <- fxn_tbl[fxn_tbl$sample != "ctl_1", ]
  expect_true(any(fxn_lab$ci_lower_fraction_highTC == 0) || any(fxn_lab$ci_upper_fraction_highTC == 1))

  # With parameter = NULL, AssessGoF prioritizes log_kdeg over co-existing log_ksyn cleanly
  ezbdo <- AssessGoF(ezbdo)
  expect_s3_class(ezbdo, "EZbakRGoF")

  # Check print.EZbakRData reports `gof: 1 analysis` without listing companion `gof_residuals:`
  prt <- paste(utils::capture.output(print(ezbdo)), collapse = "\n")
  expect_match(prt, "gof: 1 analysis")
  expect_no_match(prt, "gof_residuals:")

  gof_tbl <- EZget(ezbdo, type = "gof")
  res_tbl <- EZget(ezbdo, type = "gof_residuals")

  expect_equal(nrow(gof_tbl), ngenes)
  expect_equal(nrow(res_tbl), ngenes * (nrow(metadf) - 1L))
  expect_true(all(is.finite(res_tbl$std_resid)))
  expect_true(all(is.finite(gof_tbl$chisq_lof)))
  expect_true(all(is.finite(gof_tbl$wrmse_fn)))

  # Exact degrees-of-freedom decomposition: 8 obs = 4 timepoint groups (df_rep = 4) + (4 - 1 = 3 df_lof)
  expect_true(all(gof_tbl$n_obs == 8L))
  expect_true(all(gof_tbl$n_groups == 4L))
  expect_true(all(gof_tbl$df_rep == 4L))
  expect_true(all(gof_tbl$df_lof == 3L))
  expect_true(all(gof_tbl$df_total == 7L))

  # Null calibration: well-specified mono-exponential genes (including 0/1 boundaries) do not explode
  expect_lte(sum(gof_tbl$padj_lof < 0.05, na.rm = TRUE), 2L)
  expect_lt(stats::median(gof_tbl$wrmse_fn), 0.05)

  # Asking for log_ksyn explicitly throws an informative error
  expect_error(AssessGoF(ezbdo, parameter = "log_ksyn"), "does not directly predict")
})


test_that("AssessGoF detects biphasic decay, processing delays, and non-additive interactions", {
  set.seed(202)
  tls <- c(0.5, 1, 2, 4, 8)
  nreps <- 2
  metadf <- dplyr::tibble(
    sample = paste0("s_", rep(tls, each = nreps), "_", rep(seq_len(nreps), times = length(tls))),
    tl = rep(tls, each = nreps),
    condition = "WT"
  )

  n_good <- 30L
  n_biphasic <- 5L
  n_delay <- 5L
  all_genes <- c(
    paste0("Mono_", seq_len(n_good)),
    paste0("Biphasic_", seq_len(n_biphasic)),
    paste0("Delay_", seq_len(n_delay))
  )
  k_mono <- exp(seq(log(0.08), log(0.9), length.out = n_good))

  rows_list <- list()
  for(s in seq_len(nrow(metadf))){
    t_s <- metadf$tl[s]
    fn_true <- c(
      1 - exp(-k_mono * t_s),
      # Biphasic: 55% rapidly decaying (k=3.5) + 45% ultra-stable (k=0.005)
      rep(0.55 * (1 - exp(-3.5 * t_s)) + 0.45 * (1 - exp(-0.005 * t_s)), n_biphasic),
      # Strong nuclear/processing lag: virtually 0 labeled mature RNA until t > 1.2h, then steep rise
      rep(ifelse(t_s <= 1, 0.008 * t_s, 1 - exp(-0.85 * (t_s - 1))), n_delay)
    )
    # Sample realistic binomial noise at coverage n = 400
    n_vec <- rep(400, length(all_genes))
    fn_obs <- stats::rbinom(length(all_genes), size = n_vec, prob = pmin(pmax(fn_true, 1e-4), 1 - 1e-4)) / n_vec
    rows_list[[s]] <- dplyr::tibble(
      sample = metadf$sample[s],
      XF = all_genes,
      fraction_highTC = fn_obs,
      n = n_vec
    )
  }

  ezbfo <- make_gof_fractions_obj(dplyr::bind_rows(rows_list), metadf)
  ezbfo <- suppressMessages(EstimateKinetics(ezbfo))
  ezbfo <- suppressMessages(AverageAndRegularize(ezbfo, formula_mean = ~ 1, force_optim = TRUE))
  ezbfo <- AssessGoF(ezbfo)

  gof_tbl <- EZget(ezbfo, type = "gof")
  mono_rows <- gof_tbl[grepl("^Mono_", gof_tbl$XF), ]
  bad_rows <- gof_tbl[!grepl("^Mono_", gof_tbl$XF), ]

  expect_equal(sum(bad_rows$padj_lof < 0.01), n_biphasic + n_delay)
  expect_gt(min(bad_rows$wrmse_fn), max(mono_rows$wrmse_fn))
  expect_lte(sum(mono_rows$padj_lof < 0.05), 2L)

  # Also test direct logit_fraction_highTC additive formula_mean (~ geno + treat) catching epistasis/interaction
  meta_2x2 <- dplyr::tibble(
    sample = paste0("s", seq_len(8)),
    tl = rep(2, 8),
    geno = rep(c("WT", "WT", "KO", "KO"), each = 2),
    treat = rep(c("DMSO", "Drug", "DMSO", "Drug"), each = 2)
  )
  g_2x2 <- paste0("G", seq_len(25))
  r_2x2 <- lapply(seq_len(nrow(meta_2x2)), function(i){
    is_ko <- meta_2x2$geno[i] == "KO"
    is_dr <- meta_2x2$treat[i] == "Drug"
    eta_base <- seq(-1.5, 1.0, length.out = 25) + (if(is_ko) 0.4 else 0) + (if(is_dr) -0.3 else 0)
    # Inject massive synthetic-lethal interaction in G21..G25
    if(is_ko && is_dr) eta_base[21:25] <- eta_base[21:25] + 2.8
    fn_i <- stats::rbinom(25, size = 500, prob = stats::plogis(eta_base)) / 500
    dplyr::tibble(sample = meta_2x2$sample[i], XF = g_2x2, fraction_highTC = fn_i, n = 500)
  })
  ez_2x2 <- make_gof_fractions_obj(dplyr::bind_rows(r_2x2), meta_2x2)
  ez_2x2 <- suppressMessages(AverageAndRegularize(
    ez_2x2, type = "fractions", parameter = "logit_fraction_highTC", formula_mean = ~ geno + treat
  ))
  ez_2x2 <- AssessGoF(ez_2x2)
  gof_2x2 <- EZget(ez_2x2, type = "gof")
  expect_true(all(gof_2x2$df_lof == 1L))
  is_epi <- gof_2x2$XF %in% paste0("G", 21:25)
  expect_equal(sum(gof_2x2$padj_lof[is_epi] < 0.01), 5L)
  expect_lte(sum(gof_2x2$padj_lof[!is_epi] < 0.05), 1L)
})


test_that("AssessGoF supports shortfeed, pulse-chase, and single-tl replicate dispersion diagnostics", {
  set.seed(303)
  # 1. Pulse-chase strategy
  meta_pc <- dplyr::tibble(
    sample = paste0("pc_", seq_len(6)),
    tpulse = rep(2, 6),
    tchase = rep(c(0, 2, 6), each = 2),
    condition = "WT"
  )
  genes_pc <- paste0("PCGene", seq_len(20))
  k_pc <- exp(seq(log(0.1), log(0.8), length.out = 20))
  rows_pc <- lapply(seq_len(nrow(meta_pc)), function(i){
    tp <- meta_pc$tpulse[i]
    tc <- meta_pc$tchase[i]
    fn_true <- (1 - exp(-k_pc * tp)) * exp(-k_pc * tc)
    # Make PCGene20 chase-resistant (stabilized during chase!)
    if(tc > 0) fn_true[20] <- (1 - exp(-k_pc[20] * tp)) * 0.98
    fn_obs <- stats::rbinom(20, size = 400, prob = fn_true) / 400
    dplyr::tibble(sample = meta_pc$sample[i], XF = genes_pc, fraction_highTC = fn_obs, n = 400)
  })
  ez_pc <- make_gof_fractions_obj(dplyr::bind_rows(rows_pc), meta_pc)
  ez_pc <- suppressMessages(EstimateKinetics(ez_pc, strategy = "pulse-chase"))
  ez_pc <- suppressMessages(AverageAndRegularize(ez_pc, formula_mean = ~ 1))
  ez_pc <- AssessGoF(ez_pc)
  gof_pc <- EZget(ez_pc, type = "gof")
  expect_lt(gof_pc$padj_lof[gof_pc$XF == "PCGene20"], 0.01)
  expect_s3_class(EZFitPlot(ez_pc, features_to_plot = "PCGene20"), "ggplot")

  # 2. Shortfeed strategy
  meta_sf <- dplyr::tibble(
    sample = paste0("sf_", seq_len(6)),
    tl = rep(c(0.1, 0.25, 0.5), each = 2)
  )
  rows_sf <- lapply(seq_len(nrow(meta_sf)), function(i){
    t_i <- meta_sf$tl[i]
    fn_true <- pmin(0.4 * t_i, 0.95)
    fn_obs <- stats::rbinom(15, size = 400, prob = rep(fn_true, 15)) / 400
    dplyr::tibble(sample = meta_sf$sample[i], XF = paste0("SF", seq_len(15)), fraction_highTC = fn_obs, n = 400)
  })
  ez_sf <- make_gof_fractions_obj(dplyr::bind_rows(rows_sf), meta_sf)
  ez_sf <- suppressMessages(EstimateKinetics(ez_sf, strategy = "shortfeed"))
  ez_sf <- suppressMessages(AverageAndRegularize(ez_sf, formula_mean = ~ 1))
  ez_sf <- AssessGoF(ez_sf)
  expect_equal(nrow(EZget(ez_sf, type = "gof")), 15L)

  # 3. Single-tl design: df_lof == 0, df_rep > 0 catches hyper-variable replicate outliers
  meta_1tl <- dplyr::tibble(
    sample = paste0("r_", seq_len(4)),
    tl = rep(2, 4),
    condition = rep("WT", 4)
  )
  g_1tl <- paste0("RepG", seq_len(25))
  rows_1tl <- lapply(seq_len(4), function(i){
    eta <- rep(0, 25)
    # Inject massive replicate-to-replicate swing in RepG24 and RepG25
    eta[24:25] <- if(i %% 2 == 0) 2.2 else -2.2
    fn_obs <- stats::rbinom(25, size = 400, prob = stats::plogis(eta)) / 400
    dplyr::tibble(sample = meta_1tl$sample[i], XF = g_1tl, fraction_highTC = fn_obs, n = 400)
  })
  ez_1tl <- make_gof_fractions_obj(dplyr::bind_rows(rows_1tl), meta_1tl)
  ez_1tl <- suppressMessages(EstimateKinetics(ez_1tl))
  ez_1tl <- suppressMessages(AverageAndRegularize(ez_1tl, formula_mean = ~ 1))
  ez_1tl <- AssessGoF(ez_1tl)
  gof_1tl <- EZget(ez_1tl, type = "gof")
  expect_true(all(gof_1tl$df_lof == 0L))
  expect_true(all(is.na(gof_1tl$pval_lof)))
  expect_true(all(gof_1tl$df_rep == 3L))
  is_hyper <- gof_1tl$XF %in% c("RepG24", "RepG25")
  expect_true(all(gof_1tl$padj_rep[is_hyper] < 1e-6))

  # Volcano auto-falls back to replicate dispersion when df_lof == 0
  p_vol_rep <- EZGoFPlot(ez_1tl, plot_type = "volcano")
  expect_s3_class(p_vol_rep, "ggplot")
})


test_that("AssessGoF evaluates EZDynamics models, active boundary violations, and premature-decay coverage GoF", {
  set.seed(404)
  # Reaction graph: 0 -> P -> M -> 0 (preRNA model)
  graph_pre <- matrix(
    c(0, 1, 0,
      0, 0, 2,
      3, 0, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("0", "P", "M"), c("0", "P", "M"))
  )
  # Note: XF == "__no_feature" maps to measured_specie = "GF" (pre-mRNA P),
  # and XF != "__no_feature" maps to measured_specie = "XF" (mature mRNA M).
  mtom_pre <- list(GF ~ P, XF ~ M)

  tls <- c(1, 3)
  nreps <- 2
  metadf <- dplyr::tibble(
    sample = paste0("dyn_tl", rep(tls, each = nreps), "_r", rep(seq_len(nreps), times = length(tls))),
    tl = rep(tls, each = nreps)
  )

  n_ok <- 16L
  all_gf <- c(paste0("PreOK_", seq_len(n_ok)), "ReversedPM", "PrematureDecay")
  n_tot <- length(all_gf)

  ksyn_vec <- rep(800, n_tot)
  kp_vec <- c(exp(seq(log(0.4), log(1.5), length.out = n_ok)), 0.8, 0.25)
  kd_vec <- c(exp(seq(log(0.15), log(0.6), length.out = n_ok)), 0.3, 0.50)

  rows_dyn <- list()
  for(s in seq_len(nrow(metadf))){
    t_s <- metadf$tl[s]
    for(g in seq_len(n_tot)){
      gf <- all_gf[g]
      kp <- kp_vec[g]
      kd <- kd_vec[g]
      ksyn <- ksyn_vec[g]

      if(gf == "ReversedPM"){
        # Impossible under 0 -> P -> M -> 0: mature RNA is MUCH more labeled than pre-mRNA!
        fn_P <- 0.06
        fn_M <- 0.82
        ss_P <- 400
        ss_M <- 600
      }else if(gf == "PrematureDecay"){
        # True system has 0 -> P -> M -> 0 PLUS massive nuclear decay P -> 0 (k_prem = 2.75):
        # Total P exit rate k_tot = kp + k_prem = 3.0, while only kp/(kp + k_prem) = 1/12 reaches M!
        # Fractions follow a valid (k_tot=3.0, kd=0.5) trajectory, so fractions alone look fine,
        # BUT steady-state P/M abundance ratio is 12x higher than kd / k_tot = 0.5 / 3.0!
        k_tot <- 3.0
        fn_P <- 1 - exp(-k_tot * t_s)
        # Exact mature fraction new when P exits at k_tot and matures at kp: same ODE shape!
        fn_M <- 1 - (k_tot * exp(-kd * t_s) - kd * exp(-k_tot * t_s)) / (k_tot - kd)
        ss_P <- 2400 # ksyn_big / k_tot
        ss_M <- 400  # (kp * ss_P) / kd = (0.25 * 2400) / 0.5 -> wait: 12x excess of P vs M!
        # Under 0->P->M->0 without premature decay, ss_P / ss_M MUST equal kd / k_tot = 1/6,
        # whereas here ss_P / ss_M = 6 (36-fold abundance contradiction!)
      }else{
        fn_P <- 1 - exp(-kp * t_s)
        fn_M <- 1 - (kp * exp(-kd * t_s) - kd * exp(-kp * t_s)) / (kp - kd)
        ss_P <- ksyn / kp
        ss_M <- ksyn / kd
      }

      n_P <- max(stats::rpois(1, ss_P), 30L)
      n_M <- max(stats::rpois(1, ss_M), 30L)
      obs_P <- stats::rbinom(1, n_P, pmin(pmax(fn_P, 1e-4), 1 - 1e-4)) / n_P
      obs_M <- stats::rbinom(1, n_M, pmin(pmax(fn_M, 1e-4), 1 - 1e-4)) / n_M

      rows_dyn[[length(rows_dyn) + 1L]] <- dplyr::tibble(
        sample = metadf$sample[s],
        GF = c(gf, gf),
        XF = c("__no_feature", gf),
        fraction_highTC = c(obs_P, obs_M),
        n = c(n_P, n_M)
      )
    }
  }

  ez_dyn <- make_gof_fractions_obj(dplyr::bind_rows(rows_dyn), metadf)
  ez_dyn <- suppressMessages(AverageAndRegularize(
    ez_dyn,
    type = "fractions",
    parameter = "logit_fraction_highTC",
    formula_mean = ~ tl
  ))

  # 1. Fit EZDynamics with use_coverage = TRUE
  ez_dyn_cov <- suppressMessages(EZDynamics(
    ez_dyn,
    graph = graph_pre,
    sub_features = c("GF", "XF"),
    grouping_features = "GF",
    modeled_to_measured = mtom_pre,
    use_coverage = TRUE
  ))
  ez_dyn_cov <- AssessGoF(ez_dyn_cov, type = "dynamics")
  gof_cov <- EZget(ez_dyn_cov, type = "gof", source_type = "dynamics")

  rev_row <- gof_cov[gof_cov$GF == "ReversedPM", ]
  prem_row <- gof_cov[gof_cov$GF == "PrematureDecay", ]
  ok_rows <- gof_cov[grepl("^PreOK_", gof_cov$GF), ]

  # ReversedPM fails fraction GoF massively; PrematureDecay fails joint (fraction + coverage) GoF!
  expect_lt(rev_row$padj_lof, 0.01)
  expect_gt(rev_row$wrmse_fn, 0.25)
  expect_lt(prem_row$padj_lof_joint, 0.01)
  expect_gt(prem_row$wrmse_log10_cov, max(ok_rows$wrmse_log10_cov))
  expect_lte(sum(ok_rows$padj_lof < 0.05, na.rm = TRUE), 2L)

  # 2. Fit EZDynamics with use_coverage = FALSE (tests our ODEs.R parameter_names bugfix too!)
  ez_dyn_nocov <- suppressMessages(EZDynamics(
    ez_dyn_cov,
    graph = graph_pre,
    sub_features = c("GF", "XF"),
    grouping_features = "GF",
    modeled_to_measured = mtom_pre,
    use_coverage = FALSE,
    overwrite = FALSE
  ))
  ez_dyn_nocov <- AssessGoF(ez_dyn_nocov, type = "dynamics", repeatID = 2, overwrite = FALSE)
  gof_nocov <- EZget(ez_dyn_nocov, type = "gof", repeatID = 2)
  expect_false("padj_lof_joint" %in% colnames(gof_nocov))
  expect_lt(gof_nocov$padj_lof[gof_nocov$GF == "ReversedPM"], 0.01)

  # 3a. Single-sample EZDynamics(type = "fractions"): 4 obs (2 fxn + 2 cov) vs 3 params -> df_lof_joint = 1
  meta_1s <- metadf[1, , drop = FALSE]
  fxn_1s <- EZget(ez_dyn, type = "fractions")
  fxn_1s <- fxn_1s[fxn_1s$sample == meta_1s$sample[1], , drop = FALSE]
  ez_1s <- EZbakRFractions(
    fxn_1s[, c("sample", "GF", "XF", "fraction_highTC", "logit_fraction_highTC",
               "se_logit_fraction_highTC", "ci_lower_fraction_highTC", "ci_upper_fraction_highTC", "n")],
    meta_1s
  )
  ez_1s <- suppressMessages(EZDynamics(
    ez_1s,
    type = "fractions",
    graph = graph_pre,
    sub_features = c("GF", "XF"),
    grouping_features = "GF",
    modeled_to_measured = mtom_pre
  ))
  ez_1s <- AssessGoF(ez_1s, type = "dynamics")
  gof_1s <- EZget(ez_1s, type = "gof")
  rev_1s <- gof_1s[gof_1s$GF == "ReversedPM", ]
  expect_equal(rev_1s$df_lof_joint, 1L)
  expect_lt(rev_1s$pval_lof_joint, 0.01)
  expect_gt(rev_1s$wrmse_fn, 0.25)

  # 3b. Single-sample EZDynamics(type = "fractions", use_coverage = FALSE) with weak priors:
  # ReversedPM hits active box bound on logk3 (kd -> upper_bounds)!
  ez_1s_nocov <- suppressMessages(EZDynamics(
    ez_1s,
    type = "fractions",
    graph = graph_pre,
    sub_features = c("GF", "XF"),
    grouping_features = "GF",
    modeled_to_measured = mtom_pre,
    use_coverage = FALSE,
    prior_sds = c(100, 100, 100),
    overwrite = FALSE
  ))
  ez_1s_nocov <- AssessGoF(ez_1s_nocov, type = "dynamics", repeatID = 2)
  gof_1s_nocov <- EZget(ez_1s_nocov, type = "gof", repeatID = 2)
  rev_1s_nocov <- gof_1s_nocov[gof_1s_nocov$GF == "ReversedPM", ]
  expect_gt(rev_1s_nocov$n_at_bounds, 0L)
  expect_lt(rev_1s_nocov$pval_lof, 0.01)

  # Visualizations on EZDynamics GoF build cleanly
  expect_s3_class(ggplot2::ggplot_build(EZGoFPlot(ez_dyn_cov, source_type = "dynamics", plot_type = "volcano")), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(EZGoFPlot(ez_dyn_cov, source_type = "dynamics", plot_type = "volcano", metric = "joint")), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(EZGoFPlot(ez_dyn_cov, source_type = "dynamics", plot_type = "obs_vs_pred")), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(EZGoFPlot(ez_dyn_cov, source_type = "dynamics", plot_type = "residual_trend")), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(EZFitPlot(ez_dyn_cov, source_type = "dynamics", features_to_plot = c("ReversedPM", "PreOK_1"))), "ggplot_built")
})


test_that("EZGoFPlot, EZFitPlot, print, EZget, store_residuals = FALSE, and EZconcat work end-to-end", {
  set.seed(505)
  metadf <- dplyr::tibble(
    sample = paste0("s", seq_len(8)),
    tl = rep(c(1, 4), each = 4),
    condition = rep(c("A", "A", "B", "B"), times = 2)
  )
  genes <- paste0("Gene", seq_len(20))
  rows <- lapply(seq_len(nrow(metadf)), function(i){
    k <- seq(0.1, 0.8, length.out = 20) * (if(metadf$condition[i] == "B") 1.5 else 1.0)
    fn_true <- 1 - exp(-k * metadf$tl[i])
    # Bias sample s8 upward globally for sample_qc check
    if(metadf$sample[i] == "s8") fn_true <- pmin(fn_true + 0.12, 0.98)
    fn_obs <- stats::rbinom(20, size = 300, prob = fn_true) / 300
    dplyr::tibble(sample = metadf$sample[i], XF = genes, fraction_highTC = fn_obs, n = 300)
  })
  ezbfo <- make_gof_fractions_obj(dplyr::bind_rows(rows), metadf)
  ezbfo <- suppressMessages(EstimateKinetics(ezbfo))
  ezbfo <- suppressMessages(AverageAndRegularize(ezbfo, formula_mean = ~ condition))

  # On-the-fly fallback when AssessGoF hasn't been called yet
  expect_message(p_otf <- EZGoFPlot(ezbfo, plot_type = "volcano"), "running `AssessGoF\\(\\)`")
  expect_s3_class(ggplot2::ggplot_build(p_otf), "ggplot_built")
  expect_message(f_otf <- EZFitPlot(ezbfo, top_n = 3), "running `AssessGoF\\(\\)`")
  expect_s3_class(ggplot2::ggplot_build(f_otf), "ggplot_built")

  # Explicit AssessGoF + all plot modes & scales
  ezbfo <- AssessGoF(ezbfo)
  for(pt in c("volcano", "sample_qc", "obs_vs_pred", "residual_trend")){
    g_pt <- EZGoFPlot(ezbfo, plot_type = pt)
    expect_s3_class(ggplot2::ggplot_build(g_pt), "ggplot_built")
  }
  expect_s3_class(ggplot2::ggplot_build(
    EZGoFPlot(ezbfo, plot_type = "volcano", signed = TRUE, features_to_highlight = c("Gene1", "Gene5"))
  ), "ggplot_built")
  for(sc in c("fraction", "standardized", "logit")){
    expect_s3_class(ggplot2::ggplot_build(
      EZGoFPlot(ezbfo, plot_type = "sample_qc", residual_scale = sc)
    ), "ggplot_built")
  }
  for(rb in c("worst_lof", "worst_wrmse", "worst_rep", "best_fit")){
    expect_s3_class(ggplot2::ggplot_build(EZFitPlot(ezbfo, rank_by = rb, top_n = 2)), "ggplot_built")
  }
  expect_s3_class(ggplot2::ggplot_build(EZFitPlot(ezbfo, x_axis = "condition", top_n = 2)), "ggplot_built")

  # Check store_residuals = FALSE behavior
  ez_nores <- AssessGoF(ezbfo, store_residuals = FALSE, overwrite = TRUE)
  expect_equal(length(ez_nores$gof_residuals), 0L)
  expect_error(EZFitPlot(ez_nores), "store_residuals = TRUE")
})
