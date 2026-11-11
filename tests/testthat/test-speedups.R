### Numerical parity & edge-case tests for vectorized speedups

test_that("Vectorized standard EstimateFractions matches per-feature optim (Binomial, Poisson, bounds, -s4U)", {
  set.seed(101)
  sim <- EZSimulate(nfeatures = 30, ntreatments = 2, nreps = 1, nctlreps = 1, seqdepth = 30 * 400)
  # Add boundary edge cases (all zero TC, all mutated TC == nT)
  edge_cB <- dplyr::bind_rows(
    sim$cB,
    dplyr::tibble(sample = "sample1", feature = "edge_zero", TC = 0L, nT = 30L, n = 200L),
    dplyr::tibble(sample = "sample1", feature = "edge_one", TC = 30L, nT = 30L, n = 200L)
  )
  edge_meta <- sim$metadf
  ctl_samps <- edge_meta$sample[edge_meta$tl == 0]

  for (pois_mode in c(FALSE, TRUE)) {
    ez <- suppressMessages(
      EZbakRData(edge_cB, edge_meta) %>%
        EstimateFractions(Poisson = pois_mode)
    )
    fns <- EZget(ez, type = "fractions")

    # Reference optim(tcmml) check on labeled rows
    mutrates <- ez$mutation_rates[[1]]
    cB_dt <- data.table::as.data.table(edge_cB)
    if (pois_mode) {
      cB_sum <- cB_dt[, .(nT = sum(nT * n) / sum(n), n = sum(n)), by = .(sample, feature, TC)]
    } else {
      cB_sum <- cB_dt[, .(n = sum(n)), by = .(sample, feature, TC, nT)]
    }

    labeled_fns <- fns[!(fns$sample %in% ctl_samps), ]
    ref_logit <- numeric(nrow(labeled_fns))
    ref_se <- numeric(nrow(labeled_fns))
    for (i in seq_len(nrow(labeled_fns))) {
      s_i <- labeled_fns$sample[i]
      g_i <- labeled_fns$feature[i]
      sub <- cB_sum[sample == s_i & feature == g_i]
      pn <- mutrates$pnew[mutrates$sample == s_i]
      po <- mutrates$pold[mutrates$sample == s_i]
      fit <- stats::optim(
        0, fn = EZbakR:::tcmml,
        muts = sub$TC, nucs = sub$nT, n = sub$n,
        pnew = pn, pold = po, Poisson = pois_mode,
        lower = -9, upper = 9, method = "L-BFGS-B", hessian = TRUE
      )
      ref_logit[i] <- fit$par[1]
      ref_se[i] <- sqrt(solve(fit$hessian)[1])
    }

    expect_equal(labeled_fns$logit_fraction_highTC, ref_logit, tolerance = 5e-4)
    expect_equal(labeled_fns$se_logit_fraction_highTC, ref_se, tolerance = 5e-4)

    ctl_fns <- fns[fns$sample %in% ctl_samps, ]
    expect_true(all(ctl_fns$logit_fraction_highTC == -Inf))
    expect_true(all(ctl_fns$se_logit_fraction_highTC == 0))
    expect_true(all(ctl_fns$fraction_highTC == 0))
  }
})


test_that("Vectorized hierarchical EstimateFractions matches reference tcmml 2D optim", {
  set.seed(202)
  sim <- EZSimulate(
    nfeatures = 30, ntreatments = 2, nreps = 1, nctlreps = 1,
    seqdepth = 30 * 600
  )
  ctl_samps <- sim$metadf$sample[sim$metadf$tl == 0]

  for (pois_mode in c(FALSE, TRUE)) {
    ez <- suppressMessages(
      EZbakRData(sim$cB, sim$metadf) %>%
        EstimateFractions(strategy = "hierarchical", Poisson = pois_mode, hier_readcutoff = 100)
    )
    fns <- EZget(ez, type = "fractions")
    labeled_fns <- fns[!(fns$sample %in% ctl_samps), ]
    expect_equal(nrow(labeled_fns), 60)
    expect_true(all(is.finite(labeled_fns$logit_fraction_highTC)))
    expect_true(all(is.finite(labeled_fns$se_logit_fraction_highTC) & labeled_fns$se_logit_fraction_highTC > 0))

    # Direct check of solve_tcmm_hier_dt vs reference tcmml optim on random subset
    cB_dt <- data.table::as.data.table(sim$cB)[!(sample %in% ctl_samps), .(n = sum(n)), by = .(sample, feature, TC, nT)]
    cB_dt[, `:=`(pold = 0.002, pnew_prior = -2.9, pnew_prior_sd = 0.25)]
    hier_vec <- EZbakR:::solve_tcmm_hier_dt(
      dt = cB_dt,
      features_to_analyze = "feature",
      pop_col = "TC",
      base_col = "nT",
      samples_with_no_label = character(0),
      Poisson = FALSE,
      pnew_prior_mean_col = "pnew_prior",
      pnew_prior_sd_arg = "pnew_prior_sd",
      fraction_prior_sd = 1.7
    )

    for (k in 1:10) {
      sub <- cB_dt[sample == hier_vec$sample[k] & feature == hier_vec$feature[k]]
      ref_fit <- stats::optim(
        c(0, -2), fn = EZbakR:::tcmml,
        muts = sub$TC, nucs = sub$nT, n = sub$n,
        pold = 0.002, pnew_prior_mean = -2.9, pnew_prior_sd = 0.25,
        fraction_prior_sd = 1.7, Poisson = FALSE,
        hessian = TRUE, method = "L-BFGS-B", lower = c(-9, -9), upper = c(9, 0)
      )
      ref_u <- sqrt(abs(diag(solve(ref_fit$hessian))))
      expect_equal(hier_vec$p1[k], ref_fit$par[1], tolerance = 5e-4)
      expect_equal(hier_vec$p2[k], ref_fit$par[2], tolerance = 5e-4)
      expect_equal(hier_vec$p1_u[k], ref_u[1], tolerance = 5e-4)
      expect_equal(hier_vec$p2_u[k], ref_u[2], tolerance = 5e-4)
    }
  }
})


test_that("Vectorized multi-label EstimateFractions matches fit_general_mixture (TILAC K=3 & dual K=4)", {
  set.seed(303)
  G <- 25
  make_dual_rows <- function(samp, g) {
    n_reads <- 120L
    pop <- sample(1:4, n_reads, replace = TRUE, prob = c(0.15, 0.25, 0.25, 0.35))
    nT <- sample(20:35, n_reads, replace = TRUE)
    nG <- sample(18:32, n_reads, replace = TRUE)
    tc <- rbinom(n_reads, nT, ifelse(pop %in% c(1, 2), 0.06, 0.002))
    ga <- rbinom(n_reads, nG, ifelse(pop %in% c(1, 3), 0.05, 0.002))
    dplyr::tibble(sample = samp, XF = paste0("G", g), TC = tc, nT = nT, GA = ga, nG = nG, n = 1L)
  }

  cB_dual <- dplyr::bind_rows(
    lapply(seq_len(G), function(g) make_dual_rows("s1", g)),
    dplyr::tibble(sample = "ctl", XF = c("G1", "G2"), TC = 0L, nT = 25L, GA = 0L, nG = 25L, n = 50L)
  )
  meta_dual <- dplyr::tibble(
    sample = c("s1", "ctl"),
    tl_TC = c(2, 0),
    tl_GA = c(2, 0),
    treatment = c("A", "A")
  )

  # Full 4-component dual label
  ez4 <- suppressMessages(
    EZbakRData(cB_dual, meta_dual) %>%
      EstimateFractions(Poisson = FALSE)
  )
  fns4 <- EZget(ez4, type = "fractions")

  # Compare against per-gene fit_general_mixture on labeled sample
  md4 <- EZbakR:::create_fraction_design(c("TC", "GA")) %>%
    dplyr::filter(present) %>%
    dplyr::select(-present)
  cB_agg <- data.table::as.data.table(cB_dual)[sample == "s1", .(n = sum(n)), by = .(sample, XF, TC, nT, GA, nG)]
  mr <- ez4$mutation_rates
  pnew_s1 <- c(TC = mr$TC$pnew[mr$TC$sample == "s1"], GA = mr$GA$pnew[mr$GA$sample == "s1"])
  pold_s1 <- c(TC = mr$TC$pold[mr$TC$sample == "s1"], GA = mr$GA$pold[mr$GA$sample == "s1"])

  for (g in c(1, 5, 12, 25)) {
    gname <- paste0("G", g)
    ref_out <- EZbakR:::fit_general_mixture(
      dataset = cB_agg[XF == gname],
      Poisson = FALSE,
      mutrate_design = md4,
      twocomp = FALSE,
      pnew = pnew_s1,
      pold = pold_s1,
      mutcols = c("TC", "GA"),
      basecols = c("nT", "nG"),
      sample = "s1",
      samples_with_no_label = "ctl"
    )
    row_vec <- fns4[fns4$sample == "s1" & fns4$XF == gname, names(ref_out)]
    expect_equal(as.numeric(row_vec[1, ]), unlist(ref_out, use.names = FALSE), tolerance = 5e-4)
  }

  # Also test 3-component TILAC fraction_design
  ez3 <- suppressMessages(
    EstimateFractions(ez4, fraction_design = EZbakR::tilac_fraction_design, Poisson = TRUE, overwrite = FALSE)
  )
  fns3 <- EZget(ez3, type = "fractions", fraction_design = EZbakR::tilac_fraction_design)
  expect_equal(nrow(fns3), G + 2)
})


test_that("EZbakRArrowData pushdown matches in-memory EZbakRData across MutRates, Fractions, and QC", {
  set.seed(404)
  sim <- EZSimulate(nfeatures = 25, ntreatments = 2, nreps = 1, nctlreps = 1, seqdepth = 25 * 400)
  cB_all <- sim$cB
  meta_all <- sim$metadf

  pq_dir <- tempfile("ezbakr_arrow_")
  dir.create(pq_dir)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  arrow::write_dataset(cB_all, pq_dir, partitioning = "sample")
  ds <- arrow::open_dataset(pq_dir)

  ez_mem <- suppressMessages(
    EZbakRData(cB_all, meta_all) %>%
      EstimateMutRates(pold_from_nolabel = TRUE) %>%
      EstimateFractions(Poisson = FALSE)
  )
  ez_arr <- suppressMessages(
    EZbakRArrowData(ds, meta_all) %>%
      EstimateMutRates(pold_from_nolabel = TRUE) %>%
      EstimateFractions(Poisson = FALSE)
  )

  f_mem <- EZget(ez_mem, type = "fractions") %>% dplyr::arrange(sample, feature)
  f_arr <- EZget(ez_arr, type = "fractions") %>% dplyr::arrange(sample, feature)
  expect_equal(f_arr$logit_fraction_highTC, f_mem$logit_fraction_highTC, tolerance = 1e-5)
  expect_equal(f_arr$se_logit_fraction_highTC, f_mem$se_logit_fraction_highTC, tolerance = 1e-5)
  expect_equal(f_arr$n, f_mem$n)

  raw_mem <- EZbakR:::check_raw_mutation_rates(ez_mem, "all")$table %>% dplyr::arrange(sample, muttype)
  raw_arr <- EZbakR:::check_raw_mutation_rates(ez_arr, "all")$table %>% dplyr::arrange(sample, muttype)
  expect_equal(raw_arr, raw_mem, tolerance = 1e-12)
})


test_that("Vectorized fit_ezbakR_linear_model_batch matches fit_ezbakR_linear_model to machine precision", {
  set.seed(505)
  G <- 60
  samps <- paste0("s", 1:6)
  meta <- dplyr::tibble(
    sample = samps,
    tl = rep(2, 6),
    treatment = factor(c("WT", "WT", "WT", "KO", "KO", "KO"), levels = c("WT", "KO")),
    batch = factor(c("B1", "B2", "B1", "B1", "B2", "B2"), levels = c("B1", "B2"))
  )

  kin <- tidyr::expand_grid(XF = paste0("Gene", seq_len(G)), sample = samps) %>%
    dplyr::inner_join(meta, by = "sample") %>%
    dplyr::mutate(
      log_kdeg = rnorm(dplyr::n(), -1.2, 0.35),
      se_log_kdeg = runif(dplyr::n(), 0.04, 0.18),
      normalized_reads = runif(dplyr::n(), 50, 2000)
    )

  # Test 1: sd_groups = NULL
  f1 <- stats::as.formula("log_kdeg ~ treatment + batch - 1")
  batch1 <- EZbakR:::fit_ezbakR_linear_model_batch(
    kinetics = kin, features_to_analyze = "XF", parameter = "log_kdeg",
    formula_mean = f1, sd_groups = NULL,
    coverage_col = "normalized_reads", uncertainties_col = "se_log_kdeg",
    error_if_singular = TRUE
  )
  ref1 <- kin %>%
    dplyr::group_by(XF) %>%
    dplyr::do(dplyr::tibble(parameters = list(
      EZbakR:::fit_ezbakR_linear_model(
        ., parameter = "log_kdeg", formula_mean = f1, sd_groups = NULL,
        coverage_col = "normalized_reads", uncertainties_col = "se_log_kdeg",
        error_if_singular = TRUE
      )
    ))) %>%
    tidyr::unnest_wider(parameters) %>%
    dplyr::ungroup()
  expect_equal(batch1, ref1, tolerance = 1e-12)

  # Test 2: sd_groups = "treatment"
  batch2 <- EZbakR:::fit_ezbakR_linear_model_batch(
    kinetics = kin, features_to_analyze = "XF", parameter = "log_kdeg",
    formula_mean = f1, sd_groups = "treatment",
    coverage_col = "normalized_reads", uncertainties_col = "se_log_kdeg",
    error_if_singular = TRUE
  )
  ref2 <- kin %>%
    dplyr::group_by(XF) %>%
    dplyr::do(dplyr::tibble(parameters = list(
      EZbakR:::fit_ezbakR_linear_model(
        ., parameter = "log_kdeg", formula_mean = f1, sd_groups = "treatment",
        coverage_col = "normalized_reads", uncertainties_col = "se_log_kdeg",
        error_if_singular = TRUE
      )
    ))) %>%
    tidyr::unnest_wider(parameters) %>%
    dplyr::ungroup()
  expect_equal(batch2, ref2, tolerance = 1e-12)

  # Test 3: Partial sample presence (subset of samples dropped for some features) + singular design
  kin_partial <- kin %>%
    dplyr::mutate(
      redundant = treatment,
      keep = !(XF %in% paste0("Gene", 1:10) & sample %in% c("s1", "s6"))
    ) %>%
    dplyr::filter(keep) %>%
    dplyr::select(-keep)
  f_sing <- stats::as.formula("log_kdeg ~ treatment + batch + redundant - 1")

  batch3 <- EZbakR:::fit_ezbakR_linear_model_batch(
    kinetics = kin_partial, features_to_analyze = "XF", parameter = "log_kdeg",
    formula_mean = f_sing, sd_groups = NULL,
    coverage_col = "normalized_reads", uncertainties_col = "se_log_kdeg",
    error_if_singular = FALSE
  )
  ref3 <- kin_partial %>%
    dplyr::group_by(XF) %>%
    dplyr::do(dplyr::tibble(parameters = list(
      EZbakR:::fit_ezbakR_linear_model(
        ., parameter = "log_kdeg", formula_mean = f_sing, sd_groups = NULL,
        coverage_col = "normalized_reads", uncertainties_col = "se_log_kdeg",
        error_if_singular = FALSE
      )
    ))) %>%
    tidyr::unnest_wider(parameters) %>%
    dplyr::ungroup()
  expect_equal(batch3, ref3, tolerance = 1e-12)
})


test_that("expm_ode handles distinct, defective Jordan block, stiff, complex, and 4x4 matrices accurately", {
  # 1. Exact 2x2 defective Jordan block: B = [[-2, 0], [3, -2]] -> exp(-2) * [[1, 0], [3, 1]]
  B_jordan <- matrix(c(-2, 3, 0, -2), nrow = 2, ncol = 2)
  expected_jordan <- exp(-2) * matrix(c(1, 3, 0, 1), nrow = 2, ncol = 2)
  expect_equal(EZbakR:::expm_ode(B_jordan), expected_jordan, tolerance = 1e-12)

  # 2. Extreme stiff 2x2 compartmental matrix (ratio > 10^4, previously overflowed cosh)
  B_stiff <- matrix(c(-2500, 2500, 0, -0.25), nrow = 2, ncol = 2)
  ev <- eigen(B_stiff)
  expected_stiff <- ev$vectors %*% diag(exp(ev$values)) %*% solve(ev$vectors)
  expect_equal(EZbakR:::expm_ode(B_stiff), expected_stiff, tolerance = 1e-12)

  # 3. 2x2 with complex conjugate eigenvalues: rotation-decay [[-1, -2], [2, -1]]
  B_cplx <- matrix(c(-1, 2, -2, -1), nrow = 2, ncol = 2)
  expected_cplx <- exp(-1) * matrix(c(cos(2), sin(2), -sin(2), cos(2)), nrow = 2, ncol = 2)
  expect_equal(EZbakR:::expm_ode(B_cplx), expected_cplx, tolerance = 1e-12)

  # 4. 4x4 subtlseq / preRNA+nuc2cyto transition matrix vs eigen decomposition
  B_4x4 <- matrix(c(
    -3.2,  1.8,  0.9,  0.0,
     0.0, -1.5,  0.0,  1.1,
     0.0,  0.0, -2.1,  1.6,
     0.0,  0.0,  0.0, -0.4
  ), nrow = 4, ncol = 4)
  ev4 <- eigen(B_4x4)
  expected_4x4 <- Re(ev4$vectors %*% diag(exp(ev4$values)) %*% solve(ev4$vectors))
  expect_equal(EZbakR:::expm_ode(B_4x4), expected_4x4, tolerance = 1e-12)
})


test_that("Fast fit_beta_regression, identify_groups, and closed-form NormalizeForDropout match references", {
  # 1. fit_beta_regression analytical gr + Hessian vs numerical optim(beta_r_likelihood)
  df_iso <- dplyr::tibble(
    group = c(1L, 1L, 2L, 3L, 3L),
    transcript_id = c("tx1", "tx2", "tx1", "tx1", "tx2"),
    p = c(0.7, 0.3, 1.0, 0.2, 0.8),
    fn = c(0.62, 0.62, 0.75, 0.32, 0.32),
    n = c(150, 150, 90, 210, 210)
  )
  fast_iso <- EZbakR:::fit_beta_regression(df_iso)

  X_iso <- matrix(c(0.7, 1.0, 0.2, 0.3, 0.0, 0.8), nrow = 3, ncol = 2)
  ref_opt <- stats::optim(
    par = c(0, 0),
    fn = EZbakR:::beta_r_likelihood,
    data = c(0.62, 0.75, 0.32),
    design_matrix = X_iso,
    v = c(150, 90, 210),
    method = "L-BFGS-B", lower = -9, upper = 9,
    prior_a = 0, prior_b = 1.25, hessian = TRUE
  )
  expect_equal(fast_iso$logit_fn, ref_opt$par, tolerance = 1e-5)
  expect_equal(fast_iso$se_logit_fn, sqrt(diag(solve(ref_opt$hessian))), tolerance = 1e-4)

  # 2. identify_groups Union-Find connected components
  df_conv <- dplyr::tibble(TEC = c("t1+t2", "t3", "t2+t4", "t5+t6"))
  ig <- EZbakR:::identify_groups(df_conv, "transcript_id")
  g_map <- stats::setNames(ig$feature_group_id, ig$transcript_id)
  expect_equal(unname(g_map["t1"]), unname(g_map["t2"]))
  expect_equal(unname(g_map["t1"]), unname(g_map["t4"]))
  expect_true(g_map["t3"] != g_map["t1"])
  expect_equal(unname(g_map["t5"]), unname(g_map["t6"]))
  expect_true(g_map["t5"] != g_map["t1"])

  # 3. Closed-form NormalizeForDropout vs do_norm_ll optim
  set.seed(606)
  reffn <- runif(80, 0.05, 0.92)
  true_pdo <- 0.35
  fndo <- plogis(qlogis(reffn) + log(1 - true_pdo) + rnorm(80, 0, 0.05))
  sig <- runif(80, 0.05, 0.2)

  opt_pdo <- plogis(stats::optim(
    par = -2, fn = EZbakR:::do_norm_ll,
    reffn = reffn, fndo = fndo, sig = sig,
    method = "L-BFGS-B", lower = -6, upper = 6
  )$par[1])

  w <- 1 / (sig^2)
  delta_hat <- sum(w * (qlogis(reffn) - qlogis(fndo))) / sum(w)
  delta_clamped <- min(max(delta_hat, log1p(exp(-6))), 6 + log1p(exp(-6)))
  cf_pdo <- -expm1(-delta_clamped)
  expect_equal(cf_pdo, opt_pdo, tolerance = 1e-6)
})
