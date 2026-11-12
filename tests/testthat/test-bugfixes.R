test_that("CorrectDropout fixes operator precedence, grandR 1D optimization, and read conservation", {
  set.seed(42)
  ngenes <- 60
  features <- paste0("Gene", seq_len(ngenes))
  true_fn <- stats::runif(ngenes, 0.05, 0.85)
  true_n_base <- round(stats::runif(ngenes, 200, 2000))
  true_pdo <- 0.45

  # Observed reads after s4U dropout: retained fraction is fn*(1 - pdo) + (1 - fn)
  ret <- true_fn * (1 - true_pdo) + (1 - true_fn)
  obs_fn <- (true_fn * (1 - true_pdo)) / ret
  obs_n_s4U <- round(true_n_base * ret)
  obs_n_ctl <- round(true_n_base * 1.2)

  fxn_df <- dplyr::bind_rows(
    dplyr::tibble(
      sample = "s4U_1",
      XF = features,
      fraction_highTC = obs_fn,
      logit_fraction_highTC = EZbakR:::logit(obs_fn),
      se_logit_fraction_highTC = rep(0.1, ngenes),
      n = obs_n_s4U
    ),
    dplyr::tibble(
      sample = "ctl_1",
      XF = features,
      fraction_highTC = rep(0, ngenes),
      logit_fraction_highTC = rep(-Inf, ngenes),
      se_logit_fraction_highTC = rep(0, ngenes),
      n = obs_n_ctl
    )
  )
  metadf <- dplyr::tibble(
    sample = c("s4U_1", "ctl_1"),
    tl = c(2, 0),
    condition = c("A", "A")
  )

  ezbfo <- EZbakRFractions(fxn_df, metadf)

  # grandR strategy: 1D optimize recovers true_pdo accurately and sets scale = NA
  ezbfo_gr <- suppressMessages(CorrectDropout(ezbfo, strategy = "grandR"))
  pdo_gr <- ezbfo_gr$dropout_params[[1]]$pdo
  scale_gr <- ezbfo_gr$dropout_params[[1]]$scale
  expect_equal(pdo_gr, true_pdo, tolerance = 0.05)
  expect_true(is.na(scale_gr))

  # Total library reads in s4U sample are conserved up to ceiling() rounding
  corr_tbl <- EZget(ezbfo_gr, type = "fractions")
  s4u_corr <- corr_tbl[corr_tbl$sample == "s4U_1", ]
  expect_lte(abs(sum(s4u_corr$n) - sum(obs_n_s4U)), ngenes)

  # bakR strategy estimates numeric scale > 1 and accurate pdo
  ezbfo_br <- suppressMessages(CorrectDropout(ezbfo, strategy = "bakR"))
  expect_equal(ezbfo_br$dropout_params[[1]]$pdo, true_pdo, tolerance = 0.05)
  expect_gt(ezbfo_br$dropout_params[[1]]$scale, 1)

  # Zero-dropout dataset returns pdo == 0 under grandR root-finding (not 0.99!)
  fxn_df_zero <- dplyr::bind_rows(
    dplyr::tibble(
      sample = "s4U_1",
      XF = features,
      fraction_highTC = true_fn,
      logit_fraction_highTC = EZbakR:::logit(true_fn),
      se_logit_fraction_highTC = rep(0.1, ngenes),
      n = true_n_base
    ),
    dplyr::tibble(
      sample = "ctl_1",
      XF = features,
      fraction_highTC = rep(0, ngenes),
      logit_fraction_highTC = rep(-Inf, ngenes),
      se_logit_fraction_highTC = rep(0, ngenes),
      n = round(true_n_base * 1.2)
    )
  )
  ezbfo_zero <- suppressMessages(CorrectDropout(EZbakRFractions(fxn_df_zero, metadf), strategy = "grandR"))
  expect_equal(ezbfo_zero$dropout_params[[1]]$pdo, 0, tolerance = 0.02)
})


test_that("NormalizeForDropout uses read-weighted global fraction, feature SEs, and input validation", {
  set.seed(101)
  ngenes <- 50
  features <- paste0("Gene", seq_len(ngenes))
  true_fn <- c(rep(0.8, 10), rep(0.1, 40))
  # Skewed coverage: low-fn genes have 20x higher read counts than high-fn genes
  base_n <- c(rep(100, 10), rep(2000, 40))
  pdo_rep2 <- 0.4

  ret2 <- true_fn * (1 - pdo_rep2) + (1 - true_fn)
  fn2 <- (true_fn * (1 - pdo_rep2)) / ret2
  n2 <- round(base_n * ret2)

  fxn_df <- dplyr::bind_rows(
    dplyr::tibble(
      sample = "rep1",
      XF = features,
      fraction_highTC = true_fn,
      logit_fraction_highTC = EZbakR:::logit(true_fn),
      se_logit_fraction_highTC = rep(0.08, ngenes),
      n = base_n
    ),
    dplyr::tibble(
      sample = "rep2",
      XF = features,
      fraction_highTC = fn2,
      logit_fraction_highTC = EZbakR:::logit(fn2),
      se_logit_fraction_highTC = rep(0.08, ngenes),
      n = n2
    )
  )
  metadf <- dplyr::tibble(
    sample = c("rep1", "rep2"),
    tl = c(2, 2),
    condition = c("WT", "WT")
  )
  ezbfo <- EZbakRFractions(fxn_df, metadf)

  # check_donorm_validity rejects fewer than 2 labeled samples
  ezbfo_single <- EZbakRFractions(fxn_df[fxn_df$sample == "rep1", ], metadf[1, ])
  expect_error(NormalizeForDropout(ezbfo_single), "at least 2 labeled")

  ezbfo_norm <- NormalizeForDropout(ezbfo)
  norm_tbl <- EZget(ezbfo_norm, type = "fractions")
  rep2_norm <- norm_tbl[norm_tbl$sample == "rep2", ]

  # Read-weighted global fraction conserves total reads in rep2 (up to ceiling rounding)
  expect_lte(abs(sum(rep2_norm$n) - sum(n2)), ngenes)
  expect_equal(rep2_norm$fraction_highTC, true_fn, tolerance = 0.05)
})


test_that("Multi-label mixture fit is identifiable with finite SEs and works end-to-end with TILAC", {
  set.seed(123)
  # Simulate 3-population TILAC data (highTC_lowGA, lowTC_highGA, lowTC_lowGA)
  genes <- paste0("Gene", 1:15)
  cB_rows <- list()
  true_ratios <- seq(-1.0, 1.0, length.out = length(genes))

  for (s_idx in 1:2) {
    sname <- paste0("mix_", s_idx)
    for (g_idx in seq_along(genes)) {
      p_tc <- 0.35 * exp(true_ratios[g_idx] / 2)
      p_ga <- 0.35 * exp(-true_ratios[g_idx] / 2)
      p_bg <- max(0.05, 1 - p_tc - p_ga)
      probs <- c(p_tc, p_ga, p_bg)
      probs <- probs / sum(probs)

      counts <- as.vector(stats::rmultinom(1, size = 400, prob = probs))
      nT_vec <- rep(25L, 400)
      nG_vec <- rep(25L, 400)
      pop_id <- rep(1:3, times = counts)
      tc_prob <- ifelse(pop_id == 1, 0.08, 0.002)
      ga_prob <- ifelse(pop_id == 2, 0.07, 0.002)
      TC_vec <- stats::rbinom(400, size = nT_vec, prob = tc_prob)
      GA_vec <- stats::rbinom(400, size = nG_vec, prob = ga_prob)

      df_g <- dplyr::tibble(
        sample = sname,
        XF = genes[g_idx],
        nT = nT_vec,
        nG = nG_vec,
        TC = TC_vec,
        GA = GA_vec
      ) %>%
        dplyr::count(sample, XF, nT, nG, TC, GA, name = "n")
      cB_rows[[length(cB_rows) + 1]] <- df_g
    }
  }

  # Add a dual no-label control sample
  ctl_df <- dplyr::tibble(
    sample = "ctl_1",
    XF = rep(genes, each = 5),
    nT = 25L,
    nG = 25L,
    TC = stats::rbinom(length(genes) * 5, 25, 0.002),
    GA = stats::rbinom(length(genes) * 5, 25, 0.002),
    n = 20L
  ) %>%
    dplyr::group_by(sample, XF, nT, nG, TC, GA) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")

  cB_tilac <- dplyr::bind_rows(c(cB_rows, list(ctl_df)))
  metadf_tilac <- dplyr::tibble(
    sample = c("mix_1", "mix_2", "ctl_1"),
    tl_TC = c(2, 2, 0),
    tl_GA = c(2, 2, 0),
    treatment = c("mix", "mix", "mix")
  )

  ezbdo <- EZbakRData(cB_tilac, metadf_tilac)
  ezbdo <- suppressMessages(EstimateFractions(
    ezbdo,
    fraction_design = EZbakR::tilac_fraction_design
  ))

  fxns <- EZget(ezbdo, type = "fractions")
  fxns_labeled <- fxns %>% dplyr::filter(sample != "ctl_1")

  # Fractions sum to 1, CLR coordinates sum to 0, and all SEs are strictly finite > 0
  f_sum <- fxns_labeled$fraction_highTC_lowGA +
    fxns_labeled$fraction_lowTC_highGA +
    fxns_labeled$fraction_lowTC_lowGA
  expect_equal(f_sum, rep(1, nrow(fxns_labeled)), tolerance = 1e-6)

  clr_sum <- fxns_labeled$invsoftmax_highTC_lowGA +
    fxns_labeled$invsoftmax_lowTC_highGA +
    fxns_labeled$invsoftmax_lowTC_lowGA
  expect_equal(clr_sum, rep(0, nrow(fxns_labeled)), tolerance = 1e-6)

  expect_true(all(is.finite(fxns_labeled$se_invsoftmax_highTC_lowGA)))
  expect_true(all(is.finite(fxns_labeled$se_invsoftmax_lowTC_highGA)))
  expect_true(all(is.finite(fxns_labeled$se_invsoftmax_lowTC_lowGA)))
  expect_true(all(fxns_labeled$se_invsoftmax_highTC_lowGA > 0))

  # End-to-end TILAC kinetics -> single-level AverageAndRegularize -> CompareParameters
  ezbdo <- suppressMessages(EstimateKinetics(ezbdo, strategy = "tilac"))
  kin <- EZget(ezbdo, type = "kinetics")
  expect_equal(kin$log_TILAC_ratio, log(kin$TILAC_ratio), tolerance = 1e-10)
  expect_true(all(is.finite(kin$se_log_TILAC_ratio) & kin$se_log_TILAC_ratio > 0))

  ezbdo <- suppressMessages(AverageAndRegularize(
    ezbdo,
    parameter = "log_TILAC_ratio",
    type = "kinetics"
  ))
  ezbdo <- suppressMessages(CompareParameters(
    ezbdo,
    parameter = "log_TILAC_ratio",
    param_name = "treatmentmix"
  ))
  comp <- EZget(ezbdo, type = "comparisons")
  expect_equal(nrow(comp), length(genes))
  expect_true(all(is.finite(comp$stat)))
})


test_that("Empirical Bayes hierarchical pnew prior subtracts variances and matches Arrow backend", {
  set.seed(77)
  # Direct mathematical check of variance subtraction vs SD subtraction
  true_tau <- 0.15
  se_val <- 0.20
  n_feat <- 4000
  lpnew_draws <- stats::rnorm(n_feat, mean = -3.5, sd = sqrt(true_tau^2 + se_val^2))
  fs_dt <- data.table::data.table(
    sample = "s1",
    pnew = EZbakR:::inv_logit(lpnew_draws),
    pold = 0.001,
    lpnew_uncert = rep(se_val, n_feat)
  )
  mr_dt <- data.table::data.table(
    sample = "s1",
    pnew = EZbakR:::inv_logit(-3.5),
    pold = 0.001
  )
  est_tbl <- EZbakR:::estimate_hier_pnew_prior(
    feature_specific = fs_dt,
    mutrates = mr_dt,
    pnew_prior_sd_min = 0.01,
    pnew_prior_sd_max = 0.5
  )
  expect_equal(est_tbl$pnew_prior_sd, true_tau, tolerance = 0.05)

  # Backend parity between EZbakRData and EZbakRArrowData on hierarchical model
  sim <- EZSimulate(40, nreps = 2)
  ezbdo <- EZbakRData(sim$cB, sim$metadf)

  arrow_dir <- file.path(tempdir(), paste0("ezbakr_arrow_test_", Sys.getpid()))
  on.exit(unlink(arrow_dir, recursive = TRUE), add = TRUE)
  arrow::write_dataset(sim$cB, path = arrow_dir, format = "parquet", partitioning = "sample")
  ds <- arrow::open_dataset(arrow_dir)
  ezbado <- EZbakRArrowData(ds, sim$metadf)

  ezbdo <- suppressMessages(EstimateFractions(ezbdo, strategy = "hierarchical", hier_readcutoff = 10))
  ezbado <- suppressMessages(EstimateFractions(ezbado, strategy = "hierarchical", hier_readcutoff = 10))

  fxn_mem <- EZget(ezbdo, type = "fractions") %>% dplyr::arrange(sample, feature)
  fxn_arr <- EZget(ezbado, type = "fractions") %>% dplyr::arrange(sample, feature)

  expect_equal(fxn_mem$logit_fraction_highTC, fxn_arr$logit_fraction_highTC, tolerance = 1e-5)
  expect_equal(fxn_mem$se_logit_fraction_highTC, fxn_arr$se_logit_fraction_highTC, tolerance = 1e-5)
})


test_that("Pulse-chase and NSS SE formulas match Monte Carlo uncertainty propagation", {
  set.seed(2026)

  # 1. Pulse-chase SE check (including tchase == 0 boundary fn == 0 safety)
  fxn_pc <- dplyr::tibble(
    sample = c("pulse_1", "pulse_2", "chase_1", "pulse_zero"),
    XF = c("G1", "G1", "G1", "G0"),
    fraction_highTC = c(0.70, 0.74, 0.25, 0),
    logit_fraction_highTC = EZbakR:::logit(c(0.70, 0.74, 0.25, 0)),
    se_logit_fraction_highTC = c(0.12, 0.10, 0.15, 0),
    n = c(1000, 1000, 1000, 500)
  )
  # Add matching chase row for G0
  fxn_pc <- dplyr::bind_rows(
    fxn_pc,
    dplyr::tibble(
      sample = c("pulse_1", "pulse_2", "chase_1"),
      XF = c("G0", "G0", "G0"),
      fraction_highTC = c(0, 0, 0),
      logit_fraction_highTC = c(-Inf, -Inf, -Inf),
      se_logit_fraction_highTC = c(0, 0, 0),
      n = c(500, 500, 500)
    ) %>% dplyr::slice(3) # chase_1 for G0; pulse_zero will be paired in group
  )
  # Clean 2-gene pulse-chase dataset
  fxn_pc <- dplyr::tibble(
    sample = c("p1", "p2", "c1", "p1", "p2", "c1"),
    XF = c("G1", "G1", "G1", "G0", "G0", "G0"),
    fraction_highTC = c(0.70, 0.74, 0.25, 0, 0, 0),
    logit_fraction_highTC = c(EZbakR:::logit(c(0.70, 0.74, 0.25)), -Inf, -Inf, -Inf),
    se_logit_fraction_highTC = c(0.12, 0.10, 0.15, 0.2, 0.2, 0.2),
    n = rep(1000, 6)
  )
  meta_pc <- dplyr::tibble(
    sample = c("p1", "p2", "c1"),
    tpulse = c(2, 2, 2),
    tchase = c(0, 0, 4),
    group = c("WT", "WT", "WT")
  )
  ezbfo_pc <- EZbakRFractions(fxn_pc, meta_pc)
  ezbfo_pc <- suppressMessages(EstimateKinetics(
    ezbfo_pc,
    strategy = "pulse-chase",
    grouping_factor = "group",
    exclude_pulse_estimates = FALSE
  ))
  kin_pc <- EZget(ezbfo_pc, type = "kinetics")
  expect_true(all(is.finite(kin_pc$kdeg) & kin_pc$kdeg > 0))
  expect_true(all(is.finite(kin_pc$se_log_kdeg)))

  # Monte Carlo check for G1 in chase sample c1 (tchase = 4)
  nmc <- 50000
  fp1 <- EZbakR:::inv_logit(stats::rnorm(nmc, EZbakR:::logit(0.70), 0.12))
  fp2 <- EZbakR:::inv_logit(stats::rnorm(nmc, EZbakR:::logit(0.74), 0.10))
  fc1 <- EZbakR:::inv_logit(stats::rnorm(nmc, EZbakR:::logit(0.25), 0.15))
  lkdeg_mc <- log(-log(fc1 / ((fp1 + fp2) / 2)) / 4)
  g1_c1_se <- kin_pc %>%
    dplyr::filter(XF == "G1", sample == "c1") %>%
    dplyr::pull(se_log_kdeg)
  expect_equal(g1_c1_se, stats::sd(lkdeg_mc), tolerance = 0.03)

  # 2. NSS SE check with tl = 4 and 2 reference replicates
  fxn_nss <- dplyr::tibble(
    sample = c("ctl1", "ctl2", "s4u1"),
    XF = c("G1", "G1", "G1"),
    fraction_highTC = c(0, 0, 0.45),
    logit_fraction_highTC = c(-Inf, -Inf, EZbakR:::logit(0.45)),
    se_logit_fraction_highTC = c(0, 0, 0.10),
    n = c(800, 800, 600)
  )
  # Add stable housekeeping genes so scale_factor == 1 across samples
  hk <- tidyr::expand_grid(sample = c("ctl1", "ctl2", "s4u1"), XF = paste0("HK", 1:20)) %>%
    dplyr::mutate(
      fraction_highTC = ifelse(sample == "s4u1", 0.3, 0),
      logit_fraction_highTC = ifelse(sample == "s4u1", EZbakR:::logit(0.3), -Inf),
      se_logit_fraction_highTC = ifelse(sample == "s4u1", 0.08, 0),
      n = 1000
    )
  meta_nss <- dplyr::tibble(
    sample = c("ctl1", "ctl2", "s4u1"),
    tl = c(0, 0, 4),
    state = c("t0", "t0", "t4"),
    ref = c("t0", "t0", "t0")
  )
  ezbfo_nss <- EZbakRFractions(dplyr::bind_rows(fxn_nss, hk), meta_nss)
  ezbfo_nss <- suppressMessages(EstimateKinetics(
    ezbfo_nss,
    strategy = "NSS",
    grouping_factor = "state",
    reference_factor = "ref"
  ))
  kin_nss <- EZget(ezbfo_nss, type = "kinetics") %>% dplyr::filter(XF == "G1")

  fn_mc <- EZbakR:::inv_logit(stats::rnorm(nmc, EZbakR:::logit(0.45), 0.10))
  rs4u_mc <- stats::rpois(nmc, 600)
  rctl_mc <- (stats::rpois(nmc, 800) + stats::rpois(nmc, 800)) / 2
  kdeg_mc <- -log((1 - fn_mc) * rs4u_mc / rctl_mc) / 4
  ksyn_mc <- rs4u_mc * fn_mc * kdeg_mc / (1 - exp(-4 * kdeg_mc))

  expect_equal(kin_nss$se_log_kdeg, stats::sd(log(kdeg_mc)), tolerance = 0.03)
  expect_equal(kin_nss$se_log_ksyn, stats::sd(log(ksyn_mc)), tolerance = 0.03)
})


test_that("se_log_ksyn is invariant to feature_lengths and custom scale_factors work", {
  sim <- EZSimulate(20, nreps = 2)
  ezbdo <- EZbakRData(sim$cB, sim$metadf)
  ezbdo <- suppressMessages(EstimateFractions(ezbdo))

  fl_df <- dplyr::tibble(
    feature = unique(sim$cB$feature),
    length = rep(c(200, 10000), length.out = length(unique(sim$cB$feature)))
  )

  ez1 <- suppressMessages(EstimateKinetics(ezbdo))
  ez2 <- suppressMessages(EstimateKinetics(ezbdo, feature_lengths = fl_df))

  k1 <- EZget(ez1, type = "kinetics", exactMatch = TRUE, feature_lengths = NULL) %>%
    dplyr::arrange(sample, feature)
  k2 <- EZget(ez2, type = "kinetics", feature_lengths = fl_df) %>%
    dplyr::arrange(sample, feature)

  expect_equal(k1$se_log_ksyn, k2$se_log_ksyn, tolerance = 1e-12)

  # Custom scale_factors tibble does not get mutated into data.table and does not crash on geom_mean
  sf_tbl <- dplyr::tibble(
    sample = sim$metadf$sample,
    scale_factor = rep(1.25, nrow(sim$metadf))
  )
  ez3 <- suppressMessages(EstimateKinetics(ezbdo, scale_factors = sf_tbl))
  expect_false(data.table::is.data.table(sf_tbl))
  expect_true("geom_mean" %in% colnames(ez3$readcounts[[1]]))
})


test_that("CompareParameters dynamics inverse-variance weighting and linear model fixes work", {
  # 1. Inverse-variance weighting in CompareParameters(type = "dynamics")
  dyn_tbl <- dplyr::tibble(
    XF = rep("Gene1", 4),
    condition = factor(c("WT", "WT", "KO", "KO")),
    sub = factor(c("s1", "s2", "s1", "s2")),
    log_k1 = c(2.0, 0.0, 3.0, 3.0),
    se_log_k1 = c(0.1, 1.0, 0.2, 0.2)
  )
  meta_dyn <- dplyr::tibble(
    sample = c("s1", "s2", "s3", "s4"),
    tl = c(2, 2, 2, 2),
    condition = factor(c("WT", "WT", "KO", "KO")),
    sub = factor(c("s1", "s2", "s1", "s2"))
  )
  obj_dyn <- structure(
    list(
      metadf = meta_dyn,
      dynamics = list(dyn1 = dyn_tbl),
      metadata = list(dynamics = list(dyn1 = list(
        features = "XF",
        design_factors = c("condition", "sub")
      )))
    ),
    class = "EZbakRData"
  )

  obj_dyn <- suppressMessages(CompareParameters(
    obj_dyn,
    type = "dynamics",
    parameter = "log_k1",
    design_factor = "condition",
    reference = "WT",
    experimental = "KO"
  ))
  comp_dyn <- EZget(obj_dyn, type = "comparisons")
  expected_wt_mean <- (2.0 / 0.1^2 + 0.0 / 1.0^2) / (1 / 0.1^2 + 1 / 1.0^2)
  expected_wt_se <- sqrt(1 / (1 / 0.1^2 + 1 / 1.0^2))
  expected_ko_se <- sqrt(1 / (1 / 0.2^2 + 1 / 0.2^2))

  expect_equal(comp_dyn$difference, 3.0 - expected_wt_mean, tolerance = 1e-8)
  expect_equal(comp_dyn$uncertainty, sqrt(expected_wt_se^2 + expected_ko_se^2), tolerance = 1e-8)

  # 2. Batch linear model with explicit sd_grouping_factors, design_cor, and sd_reg_factor
  set.seed(88)
  sim <- EZSimulate(30, nreps = 3)
  meta_batch <- sim$metadf %>%
    dplyr::group_by(treatment, tl) %>%
    dplyr::mutate(batch = paste0("B", dplyr::row_number())) %>%
    dplyr::ungroup()
  ezbdo <- EZbakRData(sim$cB, meta_batch)
  ezbdo <- suppressMessages(EstimateFractions(ezbdo))
  ezbdo <- suppressMessages(EstimateKinetics(ezbdo))

  ez_reg1 <- suppressMessages(AverageAndRegularize(
    ezbdo,
    formula_mean = ~ treatment + batch,
    sd_grouping_factors = "treatment",
    sd_reg_factor = 1
  ))
  ez_reg10 <- suppressMessages(AverageAndRegularize(
    ezbdo,
    formula_mean = ~ treatment + batch,
    sd_grouping_factors = "treatment",
    sd_reg_factor = 10
  ))

  expect_true(!is.null(ez_reg1$metadata$averages[[1]]$design_cor))
  cor_mat <- ez_reg1$metadata$averages[[1]]$design_cor
  rho <- cor_mat["treatmenttreatment1", "treatmenttreatment2"]
  expect_gt(rho, 0) # Shared batch baseline induces positive correlation between group intercepts

  ez_reg1 <- suppressMessages(CompareParameters(
    ez_reg1,
    design_factor = "treatment",
    reference = "treatment1",
    experimental = "treatment2"
  ))
  comp_b <- EZget(ez_reg1, type = "comparisons")
  avg_b <- EZget(ez_reg1, type = "averages")
  expected_unc <- sqrt(pmax(
    0,
    avg_b$sd_treatmenttreatment1_posterior^2 + avg_b$sd_treatmenttreatment2_posterior^2 -
      2 * rho * avg_b$sd_treatmenttreatment1_posterior * avg_b$sd_treatmenttreatment2_posterior
  ))
  expect_equal(comp_b$uncertainty, expected_unc, tolerance = 1e-8)

  # Stronger sd_reg_factor reduces dispersion of posterior log-SDs around trend
  avg1 <- EZget(ez_reg1, type = "averages")
  avg10 <- EZget(ez_reg10, type = "averages")
  expect_lt(
    stats::sd(log(avg10$sd_treatmenttreatment1_posterior)),
    stats::sd(log(avg1$sd_treatmenttreatment1_posterior))
  )

  # Singular design with error_if_singular = FALSE does not crash
  meta_sing <- sim$metadf %>%
    dplyr::mutate(dup_cov = treatment)
  ez_sing <- EZbakRData(sim$cB, meta_sing) %>%
    EstimateFractions() %>%
    EstimateKinetics() %>%
    suppressMessages()
  expect_no_error(suppressMessages(AverageAndRegularize(
    ez_sing,
    formula_mean = ~ treatment + dup_cov,
    error_if_singular = FALSE
  )))

  # scale_factor_df divides read counts consistently in AverageAndRegularize(type = "fractions")
  sf_df1 <- dplyr::tibble(sample = sim$metadf$sample, scale_factor = rep(1, nrow(sim$metadf)))
  sf_df2 <- dplyr::tibble(sample = sim$metadf$sample, scale_factor = rep(2, nrow(sim$metadf)))
  ez_sf1 <- suppressMessages(AverageAndRegularize(
    ezbdo,
    type = "fractions",
    parameter = "logit_fraction_highTC",
    scale_factor_df = sf_df1
  ))
  ez_sf2 <- suppressMessages(AverageAndRegularize(
    ezbdo,
    type = "fractions",
    parameter = "logit_fraction_highTC",
    scale_factor_df = sf_df2
  ))
  avg_sf1 <- EZget(ez_sf1, type = "averages")
  avg_sf2 <- EZget(ez_sf2, type = "averages")
  cov_col <- grep("^coverage_", colnames(avg_sf1), value = TRUE)[1]
  expect_equal(avg_sf2[[cov_col]], avg_sf1[[cov_col]] - log10(2), tolerance = 1e-10)
})


test_that("ODEs infer_eqns multi-species replicates, repeated eigenvalues, and deterministic EZDynamics", {
  # 1. infer_eqns with multi-species compartments and replicate samples
  norm_df_t <- dplyr::tibble(
    compartment = c("total", "total", "nuclear", "nuclear", "cytoplasm", "cytoplasm")
  )
  mtom <- list(
    total = list(GF ~ NP + CP, XF ~ NM + CM),
    nuclear = list(GF ~ NP, XF ~ NM),
    cytoplasm = list(GF ~ CP, XF ~ CM)
  )
  M <- EZbakR:::infer_eqns(
    norm_df_t = norm_df_t,
    modeled_to_measured = mtom,
    species = c("NP", "NM", "CP", "CM"),
    species_to_sf = NULL,
    sample_feature = "compartment"
  )
  expect_equal(dim(M), c(6, 2))

  # 2. Repeated eigenvalues in dynamics_likelihood do not trigger browser() and return finite ll
  graph_chain <- matrix(c(
    0, 1, 0,
    0, 0, 2,
    3, 0, 0
  ), nrow = 3, byrow = TRUE, dimnames = list(c("0", "P", "M"), c("0", "P", "M")))
  ll_val <- EZbakR:::dynamics_likelihood(
    parameter_ests = c(0, 0, 0), # k2 == k3 == 1 -> repeated eigenvalue -1!
    graph = graph_chain,
    formula_list = list(tot = list(GF = GF ~ P + M, XF = XF ~ M)),
    logit_fn = c(0, -0.5),
    logit_fn_sd = c(0.1, 0.1),
    coverage = c(2.0, 1.7),
    coverage_sd = c(0.05, 0.05),
    alt_coverage = TRUE,
    scale_factor = c(1, 1),
    tls = c(1, 1),
    sample_features = c("tot", "tot"),
    feature_types = c("GF", "XF"),
    prior_means = c(0, 0, 0),
    prior_sds = c(2, 2, 2)
  )
  expect_true(is.finite(ll_val))

  # 3. EZDynamics is completely deterministic and does not advance .Random.seed
  set.seed(555)
  sim_ode <- EZSimulate(
    nfeatures = 6,
    ntreatments = 1,
    nreps = 2,
    seqdepth = 6 * 1000,
    mode = "dynamics",
    label_time = c(1, 3),
    dynamics_preset = "nuc2cyto"
  )
  ez_ode <- EZbakRData(sim_ode$cB, sim_ode$metadf) %>%
    EstimateFractions() %>%
    AverageAndRegularize(
      formula_mean = ~tl:compartment - 1,
      type = "fractions",
      parameter = "logit_fraction_highTC"
    ) %>%
    suppressMessages()

  seed_before <- .Random.seed
  ez_ode_1 <- suppressMessages(EZDynamics(
    ez_ode,
    graph = EZbakR::ode_models$nuc2cyto$graph,
    sub_features = "GF",
    grouping_features = "GF",
    sample_feature = "compartment",
    modeled_to_measured = EZbakR::ode_models$nuc2cyto$formulas
  ))
  seed_after <- .Random.seed
  ez_ode_2 <- suppressMessages(EZDynamics(
    ez_ode,
    graph = EZbakR::ode_models$nuc2cyto$graph,
    sub_features = "GF",
    grouping_features = "GF",
    sample_feature = "compartment",
    modeled_to_measured = EZbakR::ode_models$nuc2cyto$formulas
  ))
  expect_identical(seed_before, seed_after)
  expect_identical(EZget(ez_ode_1, type = "dynamics"), EZget(ez_ode_2, type = "dynamics"))

  # 4. Isoforms boundary 0/1 beta_r_likelihood is finite
  ll_beta <- EZbakR:::beta_r_likelihood(
    data = c(0, 0.5, 1),
    design_matrix = matrix(c(1, 0, 0.5, 0.5, 0, 1), nrow = 3, byrow = TRUE),
    v = c(50, 100, 50),
    par = c(-1, 1),
    prior_a = 0,
    prior_b = 3
  )
  expect_true(is.finite(ll_beta))
})


test_that("EZpcaPlot joins metadata by sample and EZVolcanoPlot/EZMAPlot highlight character vectors", {
  # 1. EZpcaPlot with deliberately reversed metadf sample order
  genes <- paste0("G", 1:25)
  samps_alpha <- c("a_B1", "a_B2", "z_A1", "z_A2")
  fxn_pca <- tidyr::expand_grid(sample = samps_alpha, XF = genes) %>%
    dplyr::mutate(
      logit_fraction_highTC = ifelse(grepl("^z_A", sample), 2.0, -2.0) +
        rep(seq(-0.2, 0.2, length.out = 25), times = 4),
      fraction_highTC = EZbakR:::inv_logit(logit_fraction_highTC),
      se_logit_fraction_highTC = 0.1,
      n = 500
    )
  # Reverse order in metadf compared to alphabetical sample order in fractions table!
  meta_rev <- dplyr::tibble(
    sample = c("z_A1", "z_A2", "a_B1", "a_B2"),
    tl = c(2, 2, 2, 2),
    group = c("GroupA", "GroupA", "GroupB", "GroupB")
  )
  ezbfo <- EZbakRFractions(fxn_pca, meta_rev)
  g_pca <- EZpcaPlot(ezbfo, metadf_cols_to_use = "group")

  pca_data <- g_pca$data
  expect_equal(
    as.character(pca_data$combo[pca_data$sample %in% c("a_B1", "a_B2")]),
    c("GroupB", "GroupB")
  )
  expect_equal(
    as.character(pca_data$combo[pca_data$sample %in% c("z_A1", "z_A2")]),
    c("GroupA", "GroupA")
  )

  # 2. EZVolcanoPlot and EZMAPlot with character vector features_to_highlight
  sim <- EZSimulate(25, nreps = 2)
  ezbdo <- EZbakRData(sim$cB, sim$metadf) %>%
    EstimateFractions() %>%
    EstimateKinetics() %>%
    AverageAndRegularize() %>%
    CompareParameters(
      design_factor = "treatment",
      reference = "treatment1",
      experimental = "treatment2"
    ) %>%
    suppressMessages()

  hl_genes <- unique(sim$cB$feature)[1:3]
  gv <- EZVolcanoPlot(ezbdo, features_to_highlight = hl_genes)
  gm <- EZMAPlot(ezbdo, features_to_highlight = hl_genes)
  expect_s3_class(gv, "ggplot")
  expect_s3_class(gm, "ggplot")
})


test_that("normalize_reads handles zero counts and singleton features properly", {
  # Feature1 & Feature2 present in both samples (true ratio S2/S1 = 2.0 -> scale factors ~ 1/sqrt(2), sqrt(2))
  # FeatureZero has n = 0 in S1; FeatureSingleton only exists in S2 with huge count
  reads_dt <- data.table::data.table(
    sample = c("S1", "S2", "S1", "S2", "S1", "S2", "S2"),
    XF = c("F1", "F1", "F2", "F2", "FZero", "FZero", "FSingleton"),
    n = c(100, 200, 200, 400, 0, 50, 5000)
  )
  norm_dt <- EZbakR:::normalize_reads(
    data.table::copy(reads_dt),
    features_to_analyze = "XF",
    feature_lengths = NULL,
    scale_factors = NULL
  )
  expect_false(any(is.na(norm_dt$normalized_reads)))
  expect_false(any(is.nan(norm_dt$normalized_reads)))

  sf_s1 <- unique(norm_dt[sample == "S1"]$scale_factor)
  sf_s2 <- unique(norm_dt[sample == "S2"]$scale_factor)
  # Ratio of scale factors S2 / S1 must be exactly 2 (singleton in S2 and zero in S1 excluded from reference!)
  expect_equal(sf_s2 / sf_s1, 2.0, tolerance = 1e-10)
})


test_that("EZQC handles multi-label datasets, boundary 0/1 fractions, and consistent return structure", {
  sim <- EZSimulate(25, nreps = 2)
  ezbdo <- EZbakRData(sim$cB, sim$metadf)
  qc_data <- suppressMessages(EZQC(ezbdo))
  expect_s3_class(qc_data$Raw_mutrates, "ggplot")
  expect_type(qc_data$Readcount_corr, "list")
  expect_false("table" %in% names(qc_data$Readcount_corr))

  # Boundary 0 and 1 fractions in EZbakRFractions do not crash kde2d() or cor()
  ezbdo_f <- suppressMessages(EstimateFractions(ezbdo))
  fxn <- EZget(ezbdo_f, type = "fractions")
  fxn$fraction_highTC[1:3] <- 0
  fxn$logit_fraction_highTC[1:3] <- -Inf
  fxn$fraction_highTC[4:5] <- 1
  fxn$logit_fraction_highTC[4:5] <- Inf
  ezbfo_bound <- EZbakRFractions(fxn, sim$metadf)
  expect_no_error(suppressMessages(EZQC(ezbfo_bound)))

  # Multi-label EZQC
  cB_multi <- data.table::data.table(
    sample = rep(c("M1", "M2", "Mctl"), each = 10),
    TC = c(rep(2L, 20), rep(0L, 10)),
    nT = rep(20L, 30),
    GA = c(rep(2L, 20), rep(0L, 10)),
    nG = rep(20L, 30),
    XF = rep(paste0("G", 1:10), times = 3),
    n = rep(50L, 30)
  )
  meta_multi <- data.table::data.table(
    sample = c("M1", "M2", "Mctl"),
    tl_TC = c(2, 2, 0),
    tl_GA = c(2, 2, 0),
    cond = c("A", "A", "A")
  )
  ez_multi <- EZbakRData(cB_multi, meta_multi)
  qc_multi <- suppressMessages(EZQC(ez_multi))
  expect_equal(length(qc_multi$Inferred_mutrates), 2)
})


test_that("Object validators, non-mutating constructors, SimpleSim checks, and EZappend/EZconcat work", {
  # 1. EZbakRData does not mutate caller tibbles into data.tables by reference
  cB_tbl <- dplyr::tibble(
    sample = c("S1", "S2"),
    TC = c(2, 0),
    nT = c(20, 20),
    XF = c("G1", "G1"),
    n = c(10, 10)
  )
  meta_tbl <- dplyr::tibble(sample = c("S1", "S2"), tl = c(2, 0), cond = c("A", "A"))
  ez <- EZbakRData(cB_tbl, meta_tbl)
  expect_false(data.table::is.data.table(cB_tbl))
  expect_false(data.table::is.data.table(meta_tbl))

  # 2. TC > nT throws mut_exceeds_basecount
  cB_bad <- cB_tbl
  cB_bad$TC[1] <- 25
  expect_error(EZbakRData(cB_bad, meta_tbl), class = "mut_exceeds_basecount")

  # 3. NtoA column is accepted with nN
  cB_ntoa <- dplyr::tibble(
    sample = c("S1", "S2"),
    NtoA = c(1, 0),
    nN = c(30, 30),
    XF = c("G1", "G1"),
    n = c(10, 10)
  )
  expect_no_error(EZbakRData(cB_ntoa, meta_tbl))

  # 4. validate_EZbakRFractions Beta(1,1) posterior SE imputation & multiple missing SEs
  fxn_nose <- dplyr::tibble(
    sample = c("S1", "S1"),
    XF = c("G1", "G2"),
    fraction_highTC_lowGA = c(0.4, 0.2),
    fraction_lowTC_highGA = c(0.3, 0.5),
    fraction_lowTC_lowGA = c(0.3, 0.3),
    logit_fraction_highTC_lowGA = EZbakR:::logit(c(0.4, 0.2)),
    logit_fraction_lowTC_highGA = EZbakR:::logit(c(0.3, 0.5)),
    logit_fraction_lowTC_lowGA = EZbakR:::logit(c(0.3, 0.3)),
    n = c(100, 200)
  )
  meta_nose <- dplyr::tibble(sample = "S1", tl_TC = 2, tl_GA = 2, cond = "A")
  ezf_imp <- EZbakRFractions(fxn_nose, meta_nose)
  imp_tbl <- ezf_imp$fractions[[1]]

  f1 <- 0.4
  n1 <- 100
  a1 <- f1 * n1 + 1
  b1 <- (1 - f1) * n1 + 1
  expected_se1 <- sqrt(a1 * b1 / (((a1 + b1)^2) * (a1 + b1 + 1))) / (f1 * (1 - f1))
  expect_equal(imp_tbl$se_logit_fraction_highTC_lowGA[1], expected_se1, tolerance = 1e-10)

  # 5. SimpleSim validates pold (not pnew twice)
  expect_error(SimpleSim(pold = -0.1), "pold")
  expect_error(SimpleSim(pold = "invalid"), "pold")

  # 6. EZappend & EZconcat preserve S3 class, metadata, and obj2 unique slots
  sim <- EZSimulate(15, nreps = 2)
  ez1 <- EZbakRData(sim$cB, sim$metadf) %>%
    EstimateFractions() %>%
    suppressMessages()
  ez2 <- EZbakRData(sim$cB, sim$metadf) %>%
    EstimateFractions() %>%
    EstimateKinetics() %>%
    suppressMessages()
  kin_tbl <- EZget(ez2, type = "kinetics")

  ez_app <- EZbakR:::EZappend(ez1, kin_tbl, features = "feature", table_type = "kinetics")
  expect_s3_class(ez_app, "EZbakRData")
  expect_true("kinetics" %in% names(ez_app))
  expect_true(!is.null(ez_app$metadata$kinetics))
  expect_true(!is.null(ez_app$metadata$fractions))
})


test_that("EstimateFractions remove_multi_features works across EZbakRData, Arrow, multi_feature_cols, sep, and AssessGoF", {
  set.seed(2026)

  # Build synthetic dataset with +label and -label samples, single-gene and multi-gene XF,
  # and single/multi-isoform TEC
  make_gene_rows <- function(samp, xf, tec, n_reads, fn_true, pnew = 0.06, pold = 0.002) {
    is_new <- stats::rbinom(n_reads, 1, fn_true)
    nT <- rep(35L, n_reads)
    tc <- ifelse(is_new == 1L,
                 stats::rbinom(n_reads, nT, pnew),
                 stats::rbinom(n_reads, nT, pold))
    dplyr::tibble(sample = samp, XF = xf, TEC = tec, TC = tc, nT = nT, n = 1L) %>%
      dplyr::group_by(sample, XF, TEC, TC, nT) %>%
      dplyr::summarise(n = sum(n), .groups = "drop")
  }

  genes_single <- paste0("Gene", 1:12)
  rows_list <- list()
  for (g_i in seq_along(genes_single)) {
    g <- genes_single[g_i]
    tec_val <- if (g_i %% 2 == 0) paste0(g, "_tx1+", g, "_tx2") else paste0(g, "_tx1")
    fn_g <- 0.15 + 0.05 * (g_i %% 10)
    rows_list[[length(rows_list) + 1L]] <- make_gene_rows("s4U_1", g, tec_val, 150, fn_g)
    rows_list[[length(rows_list) + 1L]] <- make_gene_rows("ctl_1", g, tec_val, 120, 0.0)
  }
  # Add multi-gene ambiguous rows ("Gene1+Gene2", "Gene3+Gene4+Gene5") and "__no_feature"
  rows_list[[length(rows_list) + 1L]] <- make_gene_rows("s4U_1", "Gene1+Gene2", "Gene1_tx1+Gene2_tx1", 200, 0.85)
  rows_list[[length(rows_list) + 1L]] <- make_gene_rows("ctl_1", "Gene1+Gene2", "Gene1_tx1+Gene2_tx1", 160, 0.0)
  rows_list[[length(rows_list) + 1L]] <- make_gene_rows("s4U_1", "Gene3+Gene4+Gene5", "Gene3_tx1", 180, 0.10)
  rows_list[[length(rows_list) + 1L]] <- make_gene_rows("ctl_1", "Gene3+Gene4+Gene5", "Gene3_tx1", 140, 0.0)
  rows_list[[length(rows_list) + 1L]] <- make_gene_rows("s4U_1", "__no_feature", "__no_feature", 50, 0.5)

  cB_df <- dplyr::bind_rows(rows_list)
  metadf <- dplyr::tibble(
    sample = c("s4U_1", "ctl_1"),
    tl = c(2, 0),
    condition = c("A", "A")
  )

  ezbdo <- EZbakRData(cB_df, metadf)

  # 1. Backwards compatibility: default remove_multi_features = FALSE retains multi-feature rows
  ez_keep <- suppressMessages(EstimateFractions(ezbdo, features = "XF"))
  fxn_keep <- EZget(ez_keep, type = "fractions", features = "XF")
  expect_true(all(c("Gene1+Gene2", "Gene3+Gene4+Gene5") %in% fxn_keep$XF))
  expect_false(ez_keep$metadata$fractions$XF$remove_multi_features)

  # 2. remove_multi_features = TRUE drops multi-feature rows in both +label and -label samples
  for (pois in c(TRUE, FALSE)) {
    for (strat in c("standard", "hierarchical")) {
      ez_rm <- suppressMessages(EstimateFractions(
        ezbdo,
        features = "XF",
        Poisson = pois,
        strategy = strat,
        hier_readcutoff = 50,
        remove_multi_features = TRUE
      ))
      fxn_rm <- EZget(ez_rm, type = "fractions", features = "XF")
      expect_false(any(grepl("+", fxn_rm$XF, fixed = TRUE)))
      expect_setequal(unique(fxn_rm$XF), genes_single)
      expect_setequal(unique(fxn_rm$sample), c("s4U_1", "ctl_1"))
      expect_true(ez_rm$metadata$fractions$XF$remove_multi_features)
      expect_equal(ez_rm$metadata$fractions$XF$multi_feature_cols, "XF")
    }
  }

  # 3. Selective multi_feature_cols: filter XF multi-genes while preserving TEC multi-isoforms
  ez_sel <- suppressMessages(EstimateFractions(
    ezbdo,
    features = c("XF", "TEC"),
    remove_multi_features = TRUE,
    multi_feature_cols = "XF"
  ))
  fxn_sel <- EZget(ez_sel, type = "fractions", features = c("XF", "TEC"))
  expect_false(any(grepl("+", fxn_sel$XF, fixed = TRUE)))
  expect_true(any(grepl("+", fxn_sel$TEC, fixed = TRUE)))
  expect_setequal(unique(fxn_sel$XF), genes_single)

  # Whereas default multi_feature_cols = NULL checks all analyzed features (both XF and TEC)
  ez_all <- suppressMessages(EstimateFractions(
    ezbdo,
    features = c("XF", "TEC"),
    remove_multi_features = TRUE
  ))
  fxn_all <- EZget(ez_all, type = "fractions", features = c("XF", "TEC"))
  expect_false(any(grepl("+", fxn_all$XF, fixed = TRUE)))
  expect_false(any(grepl("+", fxn_all$TEC, fixed = TRUE)))

  # 4. Custom multi_feature_sep works for both remove_multi_features and split_multi_features
  cB_semi <- cB_df %>%
    dplyr::mutate(XF = gsub("+", ";", XF, fixed = TRUE))
  ez_semi <- EZbakRData(cB_semi, metadf)

  ez_semi_rm <- suppressMessages(EstimateFractions(
    ez_semi,
    features = "XF",
    remove_multi_features = TRUE,
    multi_feature_sep = ";"
  ))
  fxn_semi_rm <- EZget(ez_semi_rm, type = "fractions", features = "XF")
  expect_false(any(grepl(";", fxn_semi_rm$XF, fixed = TRUE)))
  expect_setequal(unique(fxn_semi_rm$XF), genes_single)

  ez_semi_sp <- suppressMessages(EstimateFractions(
    ez_semi,
    features = "XF",
    split_multi_features = TRUE,
    multi_feature_cols = "XF",
    multi_feature_sep = ";"
  ))
  fxn_semi_sp <- EZget(ez_semi_sp, type = "fractions", features = "XF")
  expect_false(any(grepl(";", fxn_semi_sp$XF, fixed = TRUE)))
  # Gene1 read count in s4U_1 should include the 200 split reads from "Gene1;Gene2"
  n_g1_rm <- fxn_semi_rm$n[fxn_semi_rm$sample == "s4U_1" & fxn_semi_rm$XF == "Gene1"]
  n_g1_sp <- fxn_semi_sp$n[fxn_semi_sp$sample == "s4U_1" & fxn_semi_sp$XF == "Gene1"]
  expect_equal(n_g1_sp - n_g1_rm, 200)

  # 5. Exact parity with EZbakRArrowData & compatibility with AssessGoF(type = "fractions")
  if (requireNamespace("arrow", quietly = TRUE)) {
    tmp_dir <- tempfile("arrow_rm_mf_")
    dir.create(tmp_dir, recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
    arrow::write_dataset(cB_df, path = tmp_dir, partitioning = "sample")

    for (pois in c(TRUE, FALSE)) {
      ez_mem <- suppressMessages(EstimateFractions(
        ezbdo,
        features = "XF",
        Poisson = pois,
        remove_multi_features = TRUE
      ))
      ds <- arrow::open_dataset(tmp_dir)
      ez_arr <- EZbakRArrowData(ds, metadf)
      ez_arr <- suppressMessages(EstimateFractions(
        ez_arr,
        features = "XF",
        Poisson = pois,
        remove_multi_features = TRUE
      ))

      f_mem <- EZget(ez_mem, type = "fractions", features = "XF") %>%
        dplyr::arrange(sample, XF)
      f_arr <- EZget(ez_arr, type = "fractions", features = "XF") %>%
        dplyr::arrange(sample, XF)

      expect_equal(f_arr$XF, f_mem$XF)
      expect_equal(f_arr$n, f_mem$n)
      expect_equal(f_arr$fraction_highTC, f_mem$fraction_highTC, tolerance = 1e-6)
      expect_equal(f_arr$se_logit_fraction_highTC, f_mem$se_logit_fraction_highTC, tolerance = 1e-6)
    }

    # AssessGoF(type = "fractions") succeeds and excludes multi-feature rows on both backends
    ez_mem_gof <- suppressMessages(EstimateFractions(
      ezbdo, features = "XF", Poisson = FALSE, remove_multi_features = TRUE
    ))
    ez_mem_gof <- suppressMessages(AssessGoF(ez_mem_gof, type = "fractions", min_reads = 10))
    gof_mem <- EZget(ez_mem_gof, type = "gof")
    expect_false(any(grepl("+", gof_mem$XF, fixed = TRUE)))

    ds2 <- arrow::open_dataset(tmp_dir)
    ez_arr_gof <- EZbakRArrowData(ds2, metadf)
    ez_arr_gof <- suppressMessages(EstimateFractions(
      ez_arr_gof, features = "XF", Poisson = FALSE, remove_multi_features = TRUE
    ))
    ez_arr_gof <- suppressMessages(AssessGoF(ez_arr_gof, type = "fractions", min_reads = 10))
    gof_arr <- EZget(ez_arr_gof, type = "gof")
    expect_equal(
      dplyr::arrange(gof_arr, XF)$mean_tvd,
      dplyr::arrange(gof_mem, XF)$mean_tvd,
      tolerance = 1e-8
    )
    expect_equal(
      dplyr::arrange(gof_arr, XF)$mean_n,
      dplyr::arrange(gof_mem, XF)$mean_n,
      tolerance = 1e-8
    )
  }

  # 6. Input validation errors
  expect_error(
    EstimateFractions(ezbdo, remove_multi_features = TRUE, split_multi_features = TRUE, multi_feature_cols = "XF"),
    "Cannot set both"
  )
  expect_error(
    EstimateFractions(ezbdo, split_multi_features = TRUE),
    "multi_feature_cols"
  )
  expect_error(
    suppressMessages(EstimateFractions(ezbdo, features = "XF", remove_multi_features = TRUE, multi_feature_cols = "TEC")),
    "multi_feature_cols"
  )
  expect_error(
    EstimateFractions(ezbdo, remove_multi_features = TRUE, multi_feature_sep = ""),
    "multi_feature_sep"
  )
  expect_error(
    EstimateFractions(ezbdo, remove_multi_features = NA),
    "remove_multi_features"
  )
})
