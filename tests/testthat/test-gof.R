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


test_that("AssessGoF(type = 'fractions') evaluates exact conditional Binomial 2-comp GoF without touching obj$fractions", {
  set.seed(606)
  metadf <- dplyr::tibble(
    sample = c("s_nl", "s_lb1", "s_lb2"),
    tl = c(0, 2, 2),
    condition = c("ctl", "A", "A")
  )
  pnew_true <- 0.065
  pold_true <- 0.001

  clean_genes <- paste0("Clean_", seq_len(28))
  all_genes <- c(clean_genes, "LowDepthGene", "HetSNPGene", "ParalogTailGene")

  sim_reads_1d <- function(samp, tl_s){
    rows_g <- vector("list", length(all_genes))
    for(gi in seq_along(all_genes)){
      gn <- all_genes[gi]
      N <- if(gn == "LowDepthGene" && samp == "s_lb2") 7L else 300L
      nT <- pmax(stats::rbinom(N, size = 60, prob = 0.30), 6L)

      if(tl_s == 0){
        tc <- stats::rbinom(N, size = nT, prob = pold_true)
        if(gn == "HetSNPGene"){
          has_snp <- stats::rbinom(N, size = 1, prob = 0.55) == 1L
          tc[has_snp] <- pmin(tc[has_snp] + 1L, nT[has_snp])
        }
      }else{
        fn_g <- stats::runif(1, 0.25, 0.70)
        is_new <- stats::rbinom(N, size = 1, prob = fn_g) == 1L
        tc <- ifelse(
          is_new,
          stats::rbinom(N, size = nT, prob = pnew_true),
          stats::rbinom(N, size = nT, prob = pold_true)
        )
        if(gn == "HetSNPGene"){
          # Uncalled SNP adds +2 mismatched positions on 50% of reads
          has_snp <- stats::rbinom(N, size = 1, prob = 0.50) == 1L
          tc[has_snp] <- pmin(tc[has_snp] + 2L, nT[has_snp])
        }else if(gn == "ParalogTailGene"){
          # Cross-mapping paralog reads carry 6 T->C mismatches
          is_par <- stats::rbinom(N, size = 1, prob = 0.25) == 1L
          nT[is_par] <- pmax(nT[is_par], 10L)
          tc[is_par] <- 6L
        }
      }

      dt_g <- data.table::data.table(sample = samp, XF = gn, TC = tc, nT = nT)[
        , .(n = .N), by = .(sample, XF, TC, nT)
      ]
      rows_g[[gi]] <- dt_g
    }
    data.table::rbindlist(rows_g)
  }

  cB_sim <- dplyr::bind_rows(
    sim_reads_1d("s_nl", 0),
    sim_reads_1d("s_lb1", 2),
    sim_reads_1d("s_lb2", 2)
  )

  ez <- EZbakRData(cB_sim, metadf)
  # Fit with Poisson = TRUE to verify AssessGoF still uses exact conditional Binomial from cB!
  ez <- suppressMessages(EstimateFractions(ez, Poisson = TRUE))
  fxn_snapshot <- ez$fractions[[1]]

  # Auto-detection chooses "fractions" when neither averages nor dynamics exists
  ez <- AssessGoF(ez)
  expect_identical(ez$fractions[[1]], fxn_snapshot)

  gof_def <- EZget(ez, type = "gof", source_type = "fractions")
  # LowDepthGene has only 7 reads in s_lb2 (<= min_reads = 10), so filtered out of default GoF table
  expect_false("LowDepthGene" %in% gof_def$XF)
  # Yet LowDepthGene remains untouched in obj$fractions!
  expect_true("LowDepthGene" %in% EZget(ez, type = "fractions")$XF)

  snp_row <- gof_def[gof_def$XF == "HetSNPGene", ]
  par_row <- gof_def[gof_def$XF == "ParalogTailGene", ]
  ok_rows <- gof_def[gof_def$XF %in% clean_genes, ]

  expect_lt(snp_row$padj_gof, 1e-6)
  expect_gt(snp_row$mean_tvd, 0.10)
  expect_lt(snp_row$nolabel_padj, 1e-5)

  expect_lt(par_row$padj_gof, 1e-6)
  expect_gt(par_row$mean_tvd, 0.10)
  expect_equal(par_row$dominant_failure_mode, "high_mut_tail")

  # Clean Binomial genes remain well calibrated even though EstimateFractions used Poisson = TRUE
  expect_lte(sum(ok_rows$padj_gof < 0.05 & ok_rows$mean_tvd >= 0.05, na.rm = TRUE), 2L)

  # Test min_samples = 1 and feature_sample_counts (filtering only GoF output, never obj$fractions)
  ez_ms1 <- AssessGoF(ez, type = "fractions", min_samples = 1, overwrite = FALSE)
  expect_identical(ez_ms1$fractions[[1]], fxn_snapshot)
  gof_ms1 <- EZget(ez_ms1, type = "gof", repeatID = 2)
  expect_true("LowDepthGene" %in% gof_ms1$XF)

  fsc_tbl <- dplyr::tibble(XF = c("LowDepthGene", "HetSNPGene"), nsamps = c(1L, 2L))
  ez_fsc <- AssessGoF(ez, type = "fractions", feature_sample_counts = fsc_tbl, overwrite = FALSE)
  expect_identical(ez_fsc$fractions[[1]], fxn_snapshot)
  gof_fsc <- EZget(ez_fsc, type = "gof", repeatID = 2)
  expect_setequal(gof_fsc$XF, c("LowDepthGene", "HetSNPGene"))
})


test_that("Empirical sample overdispersion calibration protects ultra-deep genes while catching real artifacts", {
  set.seed(707)
  metadf <- dplyr::tibble(sample = c("deep1", "deep2"), tl = c(2, 2))
  pold_s <- 0.001
  pnew_center <- 0.055

  mod_genes <- paste0("Mod_", seq_len(36))
  all_g <- c(mod_genes, "Housekeeping20k", "RealParalog")

  rows_list <- list()
  for(s_i in seq_len(nrow(metadf))){
    s_nm <- metadf$sample[s_i]
    for(gn in all_g){
      N <- if(gn == "Housekeeping20k") 20000L else 350L
      nT <- pmax(stats::rbinom(N, size = 55, prob = 0.30), 8L)
      # Realistic slight gene-to-gene pnew variation around 0.055
      pnew_g <- if(gn == "Housekeeping20k"){
        0.061
      }else{
        stats::plogis(stats::qlogis(pnew_center) + stats::rnorm(1, 0, 0.14))
      }
      is_new <- stats::rbinom(N, size = 1, prob = 0.45) == 1L
      tc <- ifelse(
        is_new,
        stats::rbinom(N, size = nT, prob = pnew_g),
        stats::rbinom(N, size = nT, prob = pold_s)
      )
      if(gn == "RealParalog"){
        sp <- stats::rbinom(N, size = 1, prob = 0.30) == 1L
        tc[sp] <- 5L
      }
      dt_g <- data.table::data.table(sample = s_nm, XF = gn, TC = tc, nT = nT)[
        , .(n = .N), by = .(sample, XF, TC, nT)
      ]
      rows_list[[length(rows_list) + 1L]] <- dt_g
    }
  }

  ez_od <- EZbakRData(data.table::rbindlist(rows_list), metadf)
  ez_od <- suppressMessages(EstimateFractions(ez_od, Poisson = FALSE))

  # Nominal mode rejects Housekeeping20k due to N = 20,000 magnifying tiny pnew drift
  ez_nom <- AssessGoF(ez_od, type = "fractions", dispersion_mode = "nominal")
  gof_nom <- EZget(ez_nom, type = "gof")
  hk_nom <- gof_nom[gof_nom$XF == "Housekeeping20k", ]
  expect_lt(hk_nom$pval_gof, 0.01)

  # Default empirical mode calibrates sample-level phi_s > 0, rescuing Housekeeping20k
  # while still catching RealParalog with p << 1e-6!
  ez_emp <- AssessGoF(ez_od, type = "fractions", dispersion_mode = "empirical")
  gof_emp <- EZget(ez_emp, type = "gof")
  disp_tbl <- ez_emp$metadata$gof[[1]]$sample_dispersion
  expect_true(all(disp_tbl$phi_dispersion > 0))

  hk_emp <- gof_emp[gof_emp$XF == "Housekeeping20k", ]
  par_emp <- gof_emp[gof_emp$XF == "RealParalog", ]
  expect_equal(hk_emp$dominant_failure_mode, "ok")
  expect_lt(hk_emp$mean_tvd, 0.04)
  expect_gt(hk_emp$pval_gof, hk_nom$pval_gof)
  expect_lt(par_emp$padj_gof, 1e-6)
  expect_equal(par_emp$dominant_failure_mode, "high_mut_tail")
})


test_that("Hierarchical 2-component GoF rescues feature-specific pnew shifts while catching paralog spikes", {
  set.seed(808)
  metadf <- dplyr::tibble(sample = c("h1", "h2"), tl = c(2, 2))
  norm_genes <- paste0("HNorm_", seq_len(28))
  all_g <- c(norm_genes, "ShiftedPnew", "HierParalog")

  norm_pnews <- stats::plogis(stats::qlogis(0.070) + seq(-0.18, 0.18, length.out = length(norm_genes)))
  names(norm_pnews) <- norm_genes

  rows_h <- list()
  for(s_nm in metadf$sample){
    for(gn in all_g){
      N <- if(gn == "ShiftedPnew") 1500L else 650L
      nT <- rep(28L, N)
      pnew_g <- if(gn == "ShiftedPnew") 0.036 else if(gn == "HierParalog") 0.070 else unname(norm_pnews[gn])
      fn_g <- if(gn == "ShiftedPnew") 0.80 else 0.50
      is_new <- stats::rbinom(N, size = 1, prob = fn_g) == 1L
      tc <- ifelse(
        is_new,
        stats::rbinom(N, size = nT, prob = pnew_g),
        stats::rbinom(N, size = nT, prob = 0.001)
      )
      if(gn == "HierParalog"){
        sp <- stats::rbinom(N, size = 1, prob = 0.28) == 1L
        tc[sp] <- 6L
      }
      rows_h[[length(rows_h) + 1L]] <- data.table::data.table(
        sample = s_nm, XF = gn, TC = tc, nT = nT
      )[, .(n = .N), by = .(sample, XF, TC, nT)]
    }
  }

  ez_h <- EZbakRData(data.table::rbindlist(rows_h), metadf)
  ez_std <- suppressMessages(EstimateFractions(ez_h, strategy = "standard", Poisson = FALSE))
  ez_std <- AssessGoF(ez_std, type = "fractions")
  gof_std <- EZget(ez_std, type = "gof")

  ez_hier <- suppressMessages(EstimateFractions(ez_h, strategy = "hierarchical", Poisson = FALSE))
  ez_hier <- AssessGoF(ez_hier, type = "fractions")
  gof_hier <- EZget(ez_hier, type = "gof")

  sh_std <- gof_std[gof_std$XF == "ShiftedPnew", ]
  sh_hier <- gof_hier[gof_hier$XF == "ShiftedPnew", ]
  par_hier <- gof_hier[gof_hier$XF == "HierParalog", ]

  # ShiftedPnew fails standard 2-component GoF strongly, and is rescued by hierarchical GoF!
  expect_lt(sh_std$padj_gof, 1e-4)
  expect_gt(sh_std$mean_tvd, 0.07)
  expect_lt(sh_hier$mean_tvd, 0.04)
  expect_equal(sh_hier$dominant_failure_mode, "ok")

  # True paralog spike still fails hierarchical GoF!
  expect_lt(par_hier$padj_gof, 1e-6)
})


test_that("Multi-label (>2-component) TILAC mixture GoF catches forbidden double-labeled populations", {
  set.seed(909)
  metadf <- dplyr::tibble(
    sample = c("tilac_nl", "tilac_1", "tilac_2"),
    tl_TC = c(0, 2, 2),
    tl_GA = c(0, 2, 2)
  )
  # TILAC design: (highTC, highGA) is absent (K = 3 populations across J = 2 channels)
  fdesign_tilac <- dplyr::tibble(
    TC = c(TRUE, TRUE, FALSE, FALSE),
    GA = c(TRUE, FALSE, TRUE, FALSE),
    present = c(FALSE, TRUE, TRUE, TRUE)
  )

  clean_t <- paste0("TilacOK_", seq_len(22))
  all_t <- c(clean_t, "DoubleLabeledArtifact")

  rows_t <- list()
  for(i in seq_len(nrow(metadf))){
    s_nm <- metadf$sample[i]
    is_nl <- metadf$tl_TC[i] == 0
    for(gn in all_t){
      N <- 400L
      nT <- rep(22L, N)
      nG <- rep(22L, N)
      if(is_nl){
        tc <- stats::rbinom(N, nT, 0.001)
        ga <- stats::rbinom(N, nG, 0.001)
      }else{
        comp <- sample.int(3L, size = N, replace = TRUE, prob = c(0.35, 0.35, 0.30))
        tc <- ifelse(comp == 1L, stats::rbinom(N, nT, 0.075), stats::rbinom(N, nT, 0.001))
        ga <- ifelse(comp == 2L, stats::rbinom(N, nG, 0.075), stats::rbinom(N, nG, 0.001))
        if(gn == "DoubleLabeledArtifact"){
          # Inject forbidden (highTC & highGA) chimeric doublets
          dbl <- stats::rbinom(N, 1, 0.40) == 1L
          tc[dbl] <- stats::rbinom(sum(dbl), nT[dbl], 0.12)
          ga[dbl] <- stats::rbinom(sum(dbl), nG[dbl], 0.12)
        }
      }
      rows_t[[length(rows_t) + 1L]] <- data.table::data.table(
        sample = s_nm, XF = gn, TC = tc, GA = ga, nT = nT, nG = nG
      )[, .(n = .N), by = .(sample, XF, TC, GA, nT, nG)]
    }
  }

  ez_tilac <- EZbakRData(data.table::rbindlist(rows_t), metadf)
  ez_tilac <- suppressMessages(EstimateFractions(
    ez_tilac,
    mutrate_populations = c("TC", "GA"),
    fraction_design = fdesign_tilac,
    Poisson = FALSE
  ))
  ez_tilac <- AssessGoF(ez_tilac, type = "fractions")
  gof_tilac <- EZget(ez_tilac, type = "gof")

  art_t <- gof_tilac[gof_tilac$XF == "DoubleLabeledArtifact", ]
  ok_t <- gof_tilac[gof_tilac$XF %in% clean_t, ]

  expect_lt(art_t$padj_gof, 1e-6)
  expect_gt(art_t$mean_tvd, 0.08)
  expect_equal(art_t$dominant_failure_mode, "multilabel_mismatch")
  expect_lte(sum(ok_t$dominant_failure_mode != "ok"), 1L)

  # Multi-channel EZFitPlot builds cleanly with [TC] and [GA] facets
  expect_s3_class(
    ggplot2::ggplot_build(EZFitPlot(ez_tilac, source_type = "fractions", top_n = 2)),
    "ggplot_built"
  )
})


test_that("EZbakRArrowData backend and coexisting fractions + averages GoF plots work identically", {
  skip_if_not_installed("arrow")
  set.seed(1010)
  sim <- EZSimulate(25, nreps = 2)
  # Add a -s4U control sample and an artifact feature to test all plot modes
  fcol_sim <- setdiff(colnames(sim$cB), c("sample", "TC", "nT", "n"))[1L]
  nl_cB <- sim$cB[sim$cB$sample == sim$metadf$sample[1], ]
  nl_cB$sample <- "sample_nl"
  nl_cB$TC <- 0L
  # Inject spike into first feature across labeled samples
  feat1 <- unique(sim$cB[[fcol_sim]])[1L]
  cB_mod <- dplyr::bind_rows(sim$cB, nl_cB)
  is_f1_lb <- (cB_mod[[fcol_sim]] == feat1) & (cB_mod$sample != "sample_nl") & (cB_mod$nT >= 6)
  cB_mod$TC[is_f1_lb] <- pmin(cB_mod$nT[is_f1_lb], 6L)

  meta_mod <- dplyr::bind_rows(
    sim$metadf,
    dplyr::tibble(sample = "sample_nl", tl = 0, treatment = sim$metadf$treatment[1])
  )

  # 1. In-memory EZbakRData
  ez_mem <- EZbakRData(cB_mod, meta_mod)
  ez_mem <- suppressMessages(EstimateFractions(ez_mem, Poisson = FALSE))
  ez_mem <- AssessGoF(ez_mem, type = "fractions")

  # 2. On-disk EZbakRArrowData
  pq_dir <- tempfile("ez_arrow_gof_")
  dir.create(pq_dir, recursive = TRUE)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  arrow::write_dataset(cB_mod, path = pq_dir, format = "parquet", partitioning = "sample")

  ez_arr <- EZbakRArrowData(pq_dir, meta_mod)
  ez_arr <- suppressMessages(EstimateFractions(ez_arr, Poisson = FALSE))
  ez_arr <- AssessGoF(ez_arr, type = "fractions")

  gof_mem <- EZget(ez_mem, type = "gof", source_type = "fractions")
  gof_arr <- EZget(ez_arr, type = "gof", source_type = "fractions")
  gof_mem <- gof_mem[order(gof_mem[[fcol_sim]]), ]
  gof_arr <- gof_arr[order(gof_arr[[fcol_sim]]), ]

  expect_equal(gof_arr$mean_tvd, gof_mem$mean_tvd, tolerance = 1e-10)
  expect_equal(gof_arr$pval_gof, gof_mem$pval_gof, tolerance = 1e-10)
  expect_lt(gof_arr$padj_gof[gof_arr[[fcol_sim]] == feat1], 1e-5)

  # Coexist with kinetic averages GoF in the same object + test all EZGoFPlot & EZFitPlot modes
  ez_arr <- suppressMessages(EstimateKinetics(ez_arr))
  ez_arr <- suppressMessages(AverageAndRegularize(ez_arr))
  ez_arr <- AssessGoF(ez_arr, type = "averages")

  expect_equal(length(ez_arr$gof), 2L)
  for(pt in c("volcano", "sample_qc", "obs_vs_pred", "residual_trend")){
    expect_s3_class(
      ggplot2::ggplot_build(EZGoFPlot(ez_arr, source_type = "fractions", plot_type = pt)),
      "ggplot_built"
    )
  }
  expect_s3_class(
    ggplot2::ggplot_build(EZGoFPlot(ez_arr, source_type = "fractions", plot_type = "volcano",
                                    metric = "nolabel", features_to_highlight = feat1)),
    "ggplot_built"
  )
  for(sc in c("fraction", "standardized", "logit")){
    expect_s3_class(
      ggplot2::ggplot_build(EZGoFPlot(ez_arr, source_type = "fractions", plot_type = "sample_qc",
                                      residual_scale = sc)),
      "ggplot_built"
    )
  }
  for(rb in c("worst_gof", "worst_tvd", "worst_rep", "best_fit")){
    expect_s3_class(
      ggplot2::ggplot_build(EZFitPlot(ez_arr, source_type = "fractions", rank_by = rb, top_n = 2)),
      "ggplot_built"
    )
  }
  expect_s3_class(
    ggplot2::ggplot_build(EZFitPlot(ez_arr, source_type = "fractions",
                                    features_to_plot = feat1, samples = "all")),
    "ggplot_built"
  )
  p_top3 <- EZFitPlot(ez_arr, rank_by = "worst_lof", top_n = 3, source_type = "fractions")
  strip_lvls <- levels(p_top3$layers[[1]]$data$panel_label)
  max_line_len <- max(nchar(unlist(strsplit(strip_lvls, "\n", fixed = TRUE))))
  expect_lte(max_line_len, 20L)
})

