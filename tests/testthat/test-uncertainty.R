### Comprehensive tests for Empirical Bayes prior calibration, 95% profile-likelihood
### confidence intervals, kinetics propagation, and spline + measurement-floor regularization

test_that("Empirical Bayes logit(fn) prior calibration reduces skewed-label-time shrinkage bias and matches Arrow", {
  set.seed(1101)
  G <- 120
  # Short label time tl = 0.25h where median thalf ~ 4h -> true mean logit(fn) ~ -3.14
  true_log_kdeg <- rnorm(G, mean = log(log(2) / 4), sd = 0.6)
  tl <- 0.25
  true_fn <- 1 - exp(-exp(true_log_kdeg) * tl)
  true_logit_fn <- qlogis(true_fn)
  n_reads <- sample(c(rep(25L, G / 2), rep(400L, G / 2)))
  pnew <- 0.06
  pold <- 0.002

  cB_list <- lapply(seq_len(G), function(g) {
    nr <- n_reads[g]
    is_new <- rbinom(nr, 1, true_fn[g])
    nT <- sample(22:32, nr, replace = TRUE)
    tc <- rbinom(nr, nT, ifelse(is_new == 1L, pnew, pold))
    dplyr::tibble(sample = "short_s1", feature = paste0("G", g), TC = tc, nT = nT, n = 1L)
  })
  cB_short <- dplyr::bind_rows(cB_list)
  meta_short <- dplyr::tibble(sample = "short_s1", tl = tl, group = "A")

  ez_eb <- suppressMessages(
    EZbakRData(cB_short, meta_short) %>%
      EstimateFractions(pold_est = pold, EBayes = TRUE)
  )
  ez_legacy <- suppressMessages(
    EZbakRData(cB_short, meta_short) %>%
      EstimateFractions(pold_est = pold, EBayes = FALSE)
  )

  fn_eb <- EZget(ez_eb, type = "fractions") %>%
    dplyr::mutate(idx = as.integer(sub("^G", "", feature))) %>%
    dplyr::arrange(idx)
  fn_leg <- EZget(ez_legacy, type = "fractions") %>%
    dplyr::mutate(idx = as.integer(sub("^G", "", feature))) %>%
    dplyr::arrange(idx)

  low_cov <- which(n_reads == 25L)
  rmse_eb <- sqrt(mean((fn_eb$logit_fraction_highTC[low_cov] - true_logit_fn[low_cov])^2))
  rmse_leg <- sqrt(mean((fn_leg$logit_fraction_highTC[low_cov] - true_logit_fn[low_cov])^2))
  # Legacy Normal(0, 1.7) pulls low-coverage short-pulse genes strongly upward toward fn = 0.5;
  # EBayes calibrates prior mean near -3.1 and substantially lowers upward bias & RMSE
  bias_eb <- abs(mean(fn_eb$logit_fraction_highTC[low_cov] - true_logit_fn[low_cov]))
  bias_leg <- abs(mean(fn_leg$logit_fraction_highTC[low_cov] - true_logit_fn[low_cov]))
  expect_lt(bias_eb, bias_leg)
  expect_lt(rmse_eb, rmse_leg)

  # EZbakRArrowData matches EZbakRData with EBayes = TRUE (including CI columns)
  if (requireNamespace("arrow", quietly = TRUE)) {
    pq_dir <- tempfile("ezbakr_eb_arrow_")
    dir.create(pq_dir)
    on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
    arrow::write_dataset(cB_short, pq_dir, partitioning = "sample")
    ds <- arrow::open_dataset(pq_dir)

    ez_arr <- suppressMessages(
      EZbakRArrowData(ds, meta_short) %>%
        EstimateFractions(pold_est = pold, EBayes = TRUE)
    )
    fn_arr <- EZget(ez_arr, type = "fractions") %>%
      dplyr::mutate(idx = as.integer(sub("^G", "", feature))) %>%
      dplyr::arrange(idx)

    expect_equal(fn_arr$logit_fraction_highTC, fn_eb$logit_fraction_highTC, tolerance = 1e-6)
    expect_equal(fn_arr$se_logit_fraction_highTC, fn_eb$se_logit_fraction_highTC, tolerance = 1e-6)
    expect_equal(fn_arr$ci_lower_fraction_highTC, fn_eb$ci_lower_fraction_highTC, tolerance = 1e-6)
    expect_equal(fn_arr$ci_upper_fraction_highTC, fn_eb$ci_upper_fraction_highTC, tolerance = 1e-6)
  }
})


test_that("95% Wilks profile-likelihood CIs match exact log-likelihood drops and contain Pseudo-Huber MAP", {
  set.seed(1202)
  sim <- EZSimulate(nfeatures = 40, ntreatments = 2, nreps = 1, nctlreps = 1, seqdepth = 40 * 300)
  # Add extreme high-coverage TC=0 and TC=nT edge features on a sample with high average fn
  cB_ext <- dplyr::bind_rows(
    sim$cB,
    dplyr::tibble(sample = "sample1", feature = "zero_deep", TC = 0L, nT = 30L, n = 10000L),
    dplyr::tibble(sample = "sample1", feature = "one_deep", TC = 30L, nT = 30L, n = 5000L)
  )

  crit <- 0.5 * stats::qchisq(0.95, df = 1)

  for (pois_mode in c(FALSE, TRUE)) {
    ez <- suppressMessages(
      EZbakRData(cB_ext, sim$metadf) %>%
        EstimateFractions(Poisson = pois_mode, EBayes = TRUE)
    )
    fns <- EZget(ez, type = "fractions")

    # Universal ordering invariant: 0 <= ci_lower <= fn_hat <= ci_upper <= 1
    expect_true(all(fns$ci_lower_fraction_highTC >= 0))
    expect_true(all(fns$ci_upper_fraction_highTC <= 1))
    expect_true(all(fns$ci_lower_fraction_highTC <= fns$fraction_highTC + 1e-12))
    expect_true(all(fns$fraction_highTC <= fns$ci_upper_fraction_highTC + 1e-12))

    # -s4U control samples have exact zeros
    ctl_rows <- fns[fns$sample %in% sim$metadf$sample[sim$metadf$tl == 0], ]
    expect_true(all(ctl_rows$ci_lower_fraction_highTC == 0))
    expect_true(all(ctl_rows$ci_upper_fraction_highTC == 0))

    # Verify exact Wilks log-likelihood drop against brute-force uniroot on labeled sample
    pn <- ez$mutation_rates[[1]]$pnew[ez$mutation_rates[[1]]$sample == "sample1"]
    po <- ez$mutation_rates[[1]]$pold[ez$mutation_rates[[1]]$sample == "sample1"]
    cB_s1 <- data.table::as.data.table(cB_ext)[sample == "sample1"]
    if (pois_mode) {
      cB_s1 <- cB_s1[, .(nT = sum(nT * n) / sum(n), n = sum(n)), by = .(feature, TC)]
    } else {
      cB_s1 <- cB_s1[, .(n = sum(n)), by = .(feature, TC, nT)]
    }

    ll_pure <- function(fn_val, sub) {
      if (pois_mode) {
        l1 <- stats::dpois(sub$TC, lambda = sub$nT * pn)
        l0 <- stats::dpois(sub$TC, lambda = sub$nT * po)
      } else {
        l1 <- stats::dbinom(sub$TC, size = sub$nT, prob = pn)
        l0 <- stats::dbinom(sub$TC, size = sub$nT, prob = po)
      }
      sum(sub$n * log(fn_val * l1 + (1 - fn_val) * l0))
    }

    lab_rows <- fns[fns$sample == "sample1", ]
    for (i in seq_len(nrow(lab_rows))) {
      g_i <- lab_rows$feature[i]
      sub <- cB_s1[feature == g_i]
      opt_mle <- stats::optimize(ll_pure, interval = c(0, 1), sub = sub, maximum = TRUE, tol = 1e-11)
      ll_max <- max(opt_mle$objective, ll_pure(0, sub), ll_pure(1, sub))

      ci_lo <- lab_rows$ci_lower_fraction_highTC[i]
      ci_hi <- lab_rows$ci_upper_fraction_highTC[i]

      if (ci_lo == 0) {
        expect_lte(ll_max - ll_pure(0, sub), crit + 1e-5)
      } else {
        expect_equal(ll_max - ll_pure(ci_lo, sub), crit, tolerance = 1e-4)
      }

      if (ci_hi == 1) {
        expect_lte(ll_max - ll_pure(1, sub), crit + 1e-5)
      } else {
        expect_equal(ll_max - ll_pure(ci_hi, sub), crit, tolerance = 1e-4)
      }
    }

    expect_equal(lab_rows$ci_lower_fraction_highTC[lab_rows$feature == "zero_deep"], 0)
    expect_equal(lab_rows$ci_upper_fraction_highTC[lab_rows$feature == "one_deep"], 1)
  }
})


test_that("Profile CIs propagate through hierarchical, TILAC, Dropout, Isoforms, and EstimateKinetics", {
  set.seed(1303)
  sim <- EZSimulate(nfeatures = 35, ntreatments = 2, nreps = 2, nctlreps = 1, seqdepth = 35 * 400)
  # Inject boundary features across all labeled samples so EstimateKinetics sees TC=0 and TC=nT
  lab_samps <- sim$metadf$sample[sim$metadf$tl > 0]
  cB_aug <- dplyr::bind_rows(
    sim$cB,
    dplyr::tibble(sample = lab_samps, feature = "ultra_stable", TC = 0L, nT = 30L, n = 300L),
    dplyr::tibble(sample = lab_samps, feature = "ultra_fast", TC = 30L, nT = 30L, n = 300L)
  )

  # 1. Standard EstimateFractions -> EstimateKinetics (standard, NSS, shortfeed)
  ez <- suppressMessages(
    EZbakRData(cB_aug, sim$metadf) %>%
      EstimateFractions(EBayes = TRUE) %>%
      EstimateKinetics(strategy = "standard")
  )
  kin <- EZget(ez, type = "kinetics")
  expect_true(all(c(
    "ci_lower_kdeg", "ci_upper_kdeg",
    "thalf", "ci_lower_thalf", "ci_upper_thalf",
    "ci_lower_ksyn", "ci_upper_ksyn"
  ) %in% colnames(kin)))

  expect_true(all(kin$ci_lower_kdeg <= kin$kdeg + 1e-12))
  expect_true(all(kin$kdeg <= kin$ci_upper_kdeg + 1e-12))
  expect_true(all(kin$ci_lower_thalf <= kin$thalf + 1e-12))
  expect_true(all(kin$thalf <= kin$ci_upper_thalf + 1e-12))
  expect_true(all(kin$ci_lower_ksyn <= kin$ksyn + 1e-12))
  expect_true(all(kin$ksyn <= kin$ci_upper_ksyn + 1e-12))

  # Boundary invariants: TC=0 -> ci_lower_kdeg == 0 & ci_upper_thalf == Inf;
  # TC=nT -> ci_upper_kdeg == Inf & ci_lower_thalf == 0 & se_log_kdeg floored away from 0
  kin_stable <- kin[kin$feature == "ultra_stable", ]
  expect_true(all(kin_stable$ci_lower_kdeg == 0))
  expect_true(all(is.infinite(kin_stable$ci_upper_thalf)))

  kin_fast <- kin[kin$feature == "ultra_fast", ]
  expect_true(all(is.infinite(kin_fast$ci_upper_kdeg)))
  expect_true(all(kin_fast$ci_lower_thalf == 0))
  expect_true(all(kin_fast$se_log_kdeg >= 0.5))

  # Shortfeed strategy also produces valid CIs
  ez_sf <- suppressMessages(EstimateKinetics(ez, strategy = "shortfeed", overwrite = FALSE))
  kin_sf <- ez_sf$kinetics[[2]]
  expect_true(all(kin_sf$ci_lower_thalf <= kin_sf$thalf + 1e-12))
  expect_true(all(kin_sf$thalf <= kin_sf$ci_upper_thalf + 1e-12))

  # 2. Hierarchical EstimateFractions 2D profile CIs
  ez_hier <- suppressMessages(
    EZbakRData(sim$cB, sim$metadf) %>%
      EstimateFractions(strategy = "hierarchical", hier_readcutoff = 80, overwrite = TRUE)
  )
  fn_hier <- EZget(ez_hier, type = "fractions")
  expect_true(all(fn_hier$ci_lower_fraction_highTC <= fn_hier$fraction_highTC + 1e-12))
  expect_true(all(fn_hier$fraction_highTC <= fn_hier$ci_upper_fraction_highTC + 1e-12))

  # 3. Dropout correction & normalization preserve CI ordering and shift upward monotonically
  ez_do <- suppressMessages(CorrectDropout(ez, strategy = "grandR"))
  fn_before <- ez$fractions[[1]] %>%
    dplyr::inner_join(sim$metadf, by = "sample") %>%
    dplyr::filter(tl > 0) %>%
    dplyr::arrange(sample, feature)
  fn_after <- ez_do$fractions[[1]] %>%
    dplyr::inner_join(sim$metadf, by = "sample") %>%
    dplyr::filter(tl > 0) %>%
    dplyr::arrange(sample, feature)
  expect_true(all(fn_after$ci_lower_fraction_highTC <= fn_after$fraction_highTC + 1e-12))
  expect_true(all(fn_after$fraction_highTC <= fn_after$ci_upper_fraction_highTC + 1e-12))
  expect_true(all(fn_after$ci_lower_fraction_highTC + 1e-12 >= fn_before$ci_lower_fraction_highTC))

  # 4. TILAC 3-component simplex profile CIs + TILAC ratio kinetics CIs
  G_tilac <- 20
  cB_tilac <- dplyr::bind_rows(lapply(seq_len(G_tilac), function(g) {
    pop <- sample(1:3, 150, replace = TRUE, prob = c(0.3, 0.35, 0.35))
    nT <- sample(22:30, 150, replace = TRUE)
    nG <- sample(22:30, 150, replace = TRUE)
    dplyr::tibble(
      sample = "tilac1",
      XF = paste0("TG", g),
      TC = rbinom(150, nT, ifelse(pop == 1L, 0.06, 0.002)),
      nT = nT,
      GA = rbinom(150, nG, ifelse(pop == 2L, 0.055, 0.002)),
      nG = nG,
      n = 1L
    )
  }))
  meta_tilac <- dplyr::tibble(sample = "tilac1", tl_TC = 2, tl_GA = 2, group = "A")
  ez_tilac <- suppressMessages(
    EZbakRData(cB_tilac, meta_tilac) %>%
      EstimateFractions(fraction_design = EZbakR::tilac_fraction_design) %>%
      EstimateKinetics(strategy = "tilac")
  )
  fn_tilac <- EZget(ez_tilac, type = "fractions")
  expect_true(all(c(
    "ci_lower_fraction_highTC_lowGA", "ci_upper_fraction_highTC_lowGA",
    "ci_lower_fraction_lowTC_highGA", "ci_upper_fraction_lowTC_highGA",
    "ci_lower_fraction_lowTC_lowGA", "ci_upper_fraction_lowTC_lowGA"
  ) %in% colnames(fn_tilac)))
  expect_true(all(fn_tilac$ci_lower_fraction_highTC_lowGA <= fn_tilac$fraction_highTC_lowGA + 1e-12))
  expect_true(all(fn_tilac$fraction_highTC_lowGA <= fn_tilac$ci_upper_fraction_highTC_lowGA + 1e-12))

  kin_tilac <- EZget(ez_tilac, type = "kinetics")
  expect_true(all(c("ci_lower_TILAC_ratio", "ci_upper_TILAC_ratio") %in% colnames(kin_tilac)))
  expect_true(all(kin_tilac$ci_lower_TILAC_ratio <= kin_tilac$TILAC_ratio + 1e-12))
  expect_true(all(kin_tilac$TILAC_ratio <= kin_tilac$ci_upper_TILAC_ratio + 1e-12))

  # 5. Isoform disambiguation CIs + S3 constructor & EZappend round-trip
  iso_sim <- SimulateIsoforms(nfeatures = 20)
  meta_iso <- dplyr::tibble(sample = "sampleA", tl = 4, condition = "A")
  ez_iso <- suppressMessages(
    EZbakRData(iso_sim$cB, meta_iso) %>%
      EstimateFractions()
  )
  ez_iso[["readcounts"]][["simulated_isoform_quant"]] <- iso_sim$ground_truth %>%
    dplyr::select(transcript_id, true_count, true_TPM) %>%
    dplyr::mutate(sample = "sampleA", effective_length = 10000) %>%
    dplyr::rename(expected_count = true_count, TPM = true_TPM)
  ez_iso <- suppressMessages(EstimateIsoformFractions(ez_iso))
  fn_iso <- EZget(ez_iso, type = "fractions", features = "transcript_id")
  expect_true(all(c("ci_lower_fraction_highTC", "ci_upper_fraction_highTC") %in% colnames(fn_iso)))
  expect_true(all(fn_iso$ci_lower_fraction_highTC <= fn_iso$fraction_highTC + 1e-12))
  expect_true(all(fn_iso$fraction_highTC <= fn_iso$ci_upper_fraction_highTC + 1e-12))

  # Validate standalone constructors and EZappend with new CI & thalf columns
  expect_silent(f_obj <- EZbakRFractions(ez$fractions[[1]], sim$metadf))
  k_obj <- suppressMessages(EZbakRKinetics(ez$kinetics[[1]], sim$metadf, features = "feature"))
  expect_false(any(grepl("^se_ci_|^se_thalf$", colnames(k_obj$kinetics[[1]]))))
  ez_rebuilt <- EZbakRData(cB_aug, sim$metadf)
  ez_rebuilt <- EZbakR:::EZappend(ez_rebuilt, ez$fractions[[1]], features = "feature", table_type = "fractions")
  ez_rebuilt <- EZbakR:::EZappend(ez_rebuilt, ez$kinetics[[1]], features = "feature", table_type = "kinetics")
  expect_true(length(ez_rebuilt$fractions) == 1L)
  expect_true(length(ez_rebuilt$kinetics) == 1L)
})


test_that("AverageAndRegularize enforces measurement SE floor across all LM paths and suppresses out-of-range FPs", {
  set.seed(1404)
  # Simulate 3v3 experiment at tl = 2h with:
  # - 80 ultra-stable nulls (thalf = 18h in both WT and KO, with 3x coverage difference to stress prior/dynamic range)
  # - 40 ultra-fast saturated nulls (thalf = 0.08h in both WT and KO)
  # - 60 interior nulls (thalf = 3.5h in both WT and KO)
  # - 25 true destabilized boundary crossings (WT thalf = 18h -> KO thalf = 1.4h)
  # - 25 true interior changes (WT thalf = 4.5h -> KO thalf = 1.5h)
  tl <- 2.0
  pnew <- 0.06
  pold <- 0.0015

  thalf_wt <- c(rep(18, 80), rep(0.08, 40), rep(3.5, 60), rep(18, 25), rep(4.5, 25))
  thalf_ko <- c(rep(18, 80), rep(0.08, 40), rep(3.5, 60), rep(1.4, 25), rep(1.5, 25))
  G <- length(thalf_wt)
  gene_type <- c(
    rep("stable_null", 80),
    rep("fast_null", 40),
    rep("mid_null", 60),
    rep("boundary_de", 25),
    rep("mid_de", 25)
  )

  samps <- paste0("s", 1:6)
  conds <- c(rep("WT", 3), rep("KO", 3))
  batches <- c("B1", "B2", "B1", "B1", "B2", "B2")
  meta <- dplyr::tibble(
    sample = samps,
    tl = rep(tl, 6),
    condition = factor(conds, levels = c("WT", "KO")),
    batch = factor(batches, levels = c("B1", "B2"))
  )

  base_cov <- round(exp(runif(G, log(40), log(300))))
  cB_sim <- dplyr::bind_rows(lapply(seq_along(samps), function(s_idx) {
    is_ko <- conds[s_idx] == "KO"
    th_base <- if (is_ko) thalf_ko else thalf_wt
    # Add modest biological replicate variability (sd = 0.08 on log_kdeg)
    kdeg_rep <- (log(2) / th_base) * exp(rnorm(G, 0, 0.08))
    fn_rep <- 1 - exp(-kdeg_rep * tl)

    dplyr::bind_rows(lapply(seq_len(G), function(g) {
      # Asymmetric depth on stable_null (WT 2.5x deeper than KO) to test prior/coverage confounding!
      nr <- if (gene_type[g] == "stable_null" && !is_ko) base_cov[g] * 2L else base_cov[g]
      is_new <- rbinom(nr, 1, fn_rep[g])
      nT <- sample(24:32, nr, replace = TRUE)
      tc <- rbinom(nr, nT, ifelse(is_new == 1L, pnew, pold))
      dplyr::tibble(sample = samps[s_idx], feature = paste0("G", g), TC = tc, nT = nT, n = 1L)
    }))
  }))

  ez <- suppressMessages(
    EZbakRData(cB_sim, meta) %>%
      EstimateFractions(pold_est = pold, EBayes = TRUE) %>%
      EstimateKinetics() %>%
      AverageAndRegularize(formula_mean = ~ condition - 1, include_all_parameters = TRUE) %>%
      CompareParameters(
        parameter = "log_kdeg",
        design_factor = "condition",
        reference = "WT",
        experimental = "KO"
      )
  )

  avg_full <- ez$averages[[grep("^fullfit_", names(ez$averages))[1]]]
  comp <- EZget(ez, type = "comparisons") %>%
    dplyr::mutate(
      idx = as.integer(sub("^G", "", feature)),
      gtype = gene_type[idx]
    )

  # 1. Measurement floor invariant holds strictly on can_simply_average path
  expect_true(all(c("meas_logse_conditionWT", "meas_logse_conditionKO") %in% colnames(avg_full)))
  expect_true(all(avg_full$sd_conditionWT_posterior + 1e-12 >= exp(avg_full$meas_logse_conditionWT)))
  expect_true(all(avg_full$sd_conditionKO_posterior + 1e-12 >= exp(avg_full$meas_logse_conditionKO)))

  # 2. Dynamic-range false-positive control & true-positive power
  fp_stable <- mean(comp$padj[comp$gtype == "stable_null"] < 0.05, na.rm = TRUE)
  fp_fast <- mean(comp$padj[comp$gtype == "fast_null"] < 0.05, na.rm = TRUE)
  tp_boundary <- mean(comp$padj[comp$gtype == "boundary_de"] < 0.05, na.rm = TRUE)
  tp_mid <- mean(comp$padj[comp$gtype == "mid_de"] < 0.05, na.rm = TRUE)

  expect_lte(fp_stable, 0.10)
  expect_lte(fp_fast, 0.10)
  expect_gte(tp_boundary, 0.80)
  expect_gte(tp_mid, 0.80)

  # 3. Also verify measurement floor on single_level_mean (~ 1) and full batch OLS (~ condition + batch - 1)
  ez_single <- suppressMessages(
    AverageAndRegularize(ez, formula_mean = ~ 1, include_all_parameters = TRUE, overwrite = FALSE)
  )
  avg_single <- ez_single$averages[[length(ez_single$averages)]]
  expect_true("meas_logse_Intercept" %in% colnames(avg_single))
  expect_true(all(avg_single$sd_Intercept_posterior + 1e-12 >= exp(avg_single$meas_logse_Intercept)))

  ez_batch <- suppressMessages(
    AverageAndRegularize(ez, formula_mean = ~ condition + batch - 1, include_all_parameters = TRUE, overwrite = FALSE)
  )
  avg_batch <- ez_batch$averages[[length(ez_batch$averages)]]
  expect_true(all(c("meas_logse_conditionWT", "meas_logse_conditionKO", "meas_logse_batchB2") %in% colnames(avg_batch)))
  expect_true(all(avg_batch$sd_conditionWT_posterior + 1e-12 >= exp(avg_batch$meas_logse_conditionWT)))
  expect_true(all(avg_batch$sd_conditionKO_posterior + 1e-12 >= exp(avg_batch$meas_logse_conditionKO)))
  expect_true(all(avg_batch$sd_batchB2_posterior + 1e-12 >= exp(avg_batch$meas_logse_batchB2)))
})


test_that("Posterior precision SE reconciliation eliminates 0.8 boundary shelf, near-boundary spikes, and scales smoothly with depth", {
  set.seed(1505)
  tl <- 2.0
  pnew <- 0.06
  pold <- 0.002

  # Realistic background transcriptome (~400 interior genes with SD(log_kdeg) ~ 1.1)
  G_bg <- 400
  true_lk <- rnorm(G_bg, mean = -1.2, sd = 1.1)
  true_fn <- 1 - exp(-exp(true_lk) * tl)
  cB_bg <- dplyr::bind_rows(lapply(seq_len(G_bg), function(g) {
    nr <- 200L
    is_new <- rbinom(nr, 1, true_fn[g])
    nT <- rep(28L, nr)
    tc <- rbinom(nr, nT, ifelse(is_new == 1L, pnew, pold))
    dplyr::tibble(sample = "s1", feature = paste0("BG_", g), TC = tc, nT = nT, n = 1L)
  }))

  # 1. Boundary genes at increasing read depths (TC = 0 and TC = nT)
  depths <- c(15L, 60L, 250L, 1000L)
  cB_depths <- dplyr::bind_rows(
    lapply(seq_along(depths), function(i) {
      dplyr::tibble(
        sample = "s1",
        feature = c(paste0("zero_d", i), paste0("one_d", i)),
        TC = c(0L, 28L),
        nT = c(28L, 28L),
        n = c(depths[i], depths[i])
      )
    })
  )

  # 2. Single-read sweeps right across the 95% Wilks boundary cutoffs (crit = 1.920729):
  # Lower sweep: 3 reads of (TC = 2, nT = 28) + M reads of (TC = 0, nT = 28) (crosses around M = 163)
  # Upper sweep: 25 reads of (TC = 0, nT = 28) + M reads of (TC = 4, nT = 28) (crosses around M = 64)
  M_grid_lo <- 150:175
  cB_sweep_lo <- dplyr::bind_rows(lapply(seq_along(M_grid_lo), function(i) {
    dplyr::tibble(
      sample = "s1",
      feature = paste0("sweep_lo_", sprintf("%03d", i)),
      TC = c(2L, 0L),
      nT = c(28L, 28L),
      n = c(3L, M_grid_lo[i])
    )
  }))
  M_grid_hi <- 52:77
  cB_sweep_hi <- dplyr::bind_rows(lapply(seq_along(M_grid_hi), function(i) {
    dplyr::tibble(
      sample = "s1",
      feature = paste0("sweep_hi_", sprintf("%03d", i)),
      TC = c(0L, 4L),
      nT = c(28L, 28L),
      n = c(25L, M_grid_hi[i])
    )
  }))

  meta <- dplyr::tibble(sample = "s1", tl = tl, condition = "A")
  ez <- suppressMessages(
    EZbakRData(dplyr::bind_rows(cB_bg, cB_depths, cB_sweep_lo, cB_sweep_hi), meta) %>%
      EstimateFractions(
        pold_est = pold,
        pnew_prior_mean = stats::qlogis(pnew),
        pnew_prior_sd = 1e-4,
        EBayes = TRUE
      ) %>%
      EstimateKinetics()
  )

  fn <- EZget(ez, type = "fractions")
  kin <- EZget(ez, type = "kinetics")

  # Check 1: Boundary CIs remain exact (0 and Inf), no flat 0.8 shelf, and SE strictly decreases with depth
  kin_z <- kin[match(paste0("zero_d", seq_along(depths)), kin$feature), ]
  kin_o <- kin[match(paste0("one_d", seq_along(depths)), kin$feature), ]

  expect_true(all(kin_z$ci_lower_kdeg == 0))
  expect_true(all(is.infinite(kin_o$ci_upper_kdeg)))
  expect_false(any(abs(kin_z$se_log_kdeg - 0.8) < 1e-12))
  expect_false(any(abs(kin_o$se_log_kdeg - 0.8) < 1e-12))
  expect_true(all(diff(kin_z$se_log_kdeg) < 0))
  expect_true(all(diff(kin_o$se_log_kdeg) < 0))
  expect_true(all(kin_z$se_log_kdeg >= 0.5))
  expect_true(all(kin_o$se_log_kdeg >= 0.5))

  # Check 2: Zero extreme near-boundary spikes across entire dataset
  expect_true(all(is.finite(fn$se_logit_fraction_highTC)))
  expect_true(all(is.finite(kin$se_log_kdeg)))
  expect_lte(max(fn$se_logit_fraction_highTC), 2.5)
  expect_lte(max(kin$se_log_kdeg), 2.0)

  # Check 3: Exact mathematical continuity of posterior_se_from_ci across +/-Inf
  ci_lo_cont <- c(seq(-3.5, -28, length.out = 1000), -Inf)
  se_lo_cont <- EZbakR:::posterior_se_from_ci(
    est = rep(-3.0, length(ci_lo_cont)),
    ci_lo = ci_lo_cont,
    ci_hi = rep(-2.2, length(ci_lo_cont)),
    prior_sd = 1.2
  )
  expect_lt(max(abs(diff(se_lo_cont))), 0.01)

  ci_hi_cont <- c(seq(1.5, 28, length.out = 1000), Inf)
  se_hi_cont <- EZbakR:::posterior_se_from_ci(
    est = rep(1.0, length(ci_hi_cont)),
    ci_lo = rep(0.4, length(ci_hi_cont)),
    ci_hi = ci_hi_cont,
    prior_sd = 1.2
  )
  expect_lt(max(abs(diff(se_hi_cont))), 0.01)

  # Check 4: End-to-end single-read sweeps across Wilks cutoffs glide monotonically
  # into the boundary without spikes or overshoot-and-snap reversals
  sw_lo_kin <- kin[grep("^sweep_lo_", kin$feature), ]
  sw_lo_kin <- sw_lo_kin[order(sw_lo_kin$feature), ]
  expect_true(any(sw_lo_kin$ci_lower_kdeg > 0) && any(sw_lo_kin$ci_lower_kdeg == 0))
  idx_cross_lo <- which(sw_lo_kin$ci_lower_kdeg == 0)[1]
  expect_true(all(diff(sw_lo_kin$se_log_kdeg[seq_len(idx_cross_lo)]) > 0))
  expect_lt(max(abs(diff(sw_lo_kin$se_log_kdeg))), 0.06)

  sw_hi_kin <- kin[grep("^sweep_hi_", kin$feature), ]
  sw_hi_kin <- sw_hi_kin[order(sw_hi_kin$feature), ]
  expect_true(any(is.finite(sw_hi_kin$ci_upper_kdeg)) && any(is.infinite(sw_hi_kin$ci_upper_kdeg)))
  idx_cross_hi <- which(is.infinite(sw_hi_kin$ci_upper_kdeg))[1]
  expect_true(all(diff(sw_hi_kin$se_log_kdeg[seq_len(idx_cross_hi)]) > 0))
  expect_lt(max(abs(diff(sw_hi_kin$se_log_kdeg))), 0.20)
})
