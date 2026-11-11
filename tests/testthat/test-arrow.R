make_arrow_sim <- function(nfeatures = 30, seed = 42) {
  set.seed(seed)
  sim <- EZSimulate(nfeatures = nfeatures, nreps = 2)
  sim
}


test_that("WriteCBArrow writes partitioned integer Parquet from data.frame, EZbakRData, and CSV", {
  sim <- make_arrow_sim(20, seed = 101)
  cB_dbl <- as.data.frame(sim$cB)
  cB_dbl$TC <- as.double(cB_dbl$TC)
  cB_dbl$nT <- as.double(cB_dbl$nT)
  cB_dbl$n <- as.double(cB_dbl$n)

  pq_dir1 <- file.path(tempdir(), paste0("ezbakr_wcb_df_", Sys.getpid()))
  unlink(pq_dir1, recursive = TRUE)
  on.exit(unlink(pq_dir1, recursive = TRUE), add = TRUE)

  out_path <- WriteCBArrow(cB_dbl, pq_dir1)
  expect_true(dir.exists(out_path))
  # Hive subdirectories sample=... should exist
  subdirs <- list.dirs(pq_dir1, full.names = FALSE, recursive = FALSE)
  expect_true(all(grepl("^sample=", subdirs)))

  ds1 <- arrow::open_dataset(pq_dir1)
  sch1 <- arrow::schema(ds1)
  expect_equal(sch1[["TC"]]$type$ToString(), "int32")
  expect_equal(sch1[["nT"]]$type$ToString(), "int32")
  expect_equal(sch1[["n"]]$type$ToString(), "int32")

  # Overwrite protection
  expect_error(WriteCBArrow(cB_dbl, pq_dir1, overwrite = FALSE),
               class = "arrow_dir_exists")
  expect_silent(WriteCBArrow(cB_dbl, pq_dir1, overwrite = TRUE))

  # From EZbakRData
  ez_mem <- EZbakRData(sim$cB, sim$metadf)
  pq_dir2 <- file.path(tempdir(), paste0("ezbakr_wcb_ez_", Sys.getpid()))
  unlink(pq_dir2, recursive = TRUE)
  on.exit(unlink(pq_dir2, recursive = TRUE), add = TRUE)
  WriteCBArrow(ez_mem, pq_dir2)
  expect_equal(nrow( dplyr::collect(arrow::open_dataset(pq_dir2)) ), nrow(ez_mem$cB))

  # From CSV files
  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  utils::write.csv(sim$cB, csv_file, row.names = FALSE)

  pq_dir3 <- file.path(tempdir(), paste0("ezbakr_wcb_csv_", Sys.getpid()))
  unlink(pq_dir3, recursive = TRUE)
  on.exit(unlink(pq_dir3, recursive = TRUE), add = TRUE)
  WriteCBArrow(csv_file, pq_dir3)
  expect_equal(nrow(dplyr::collect(arrow::open_dataset(pq_dir3))), nrow(sim$cB))
})


test_that("EZbakRArrowData accepts directory paths, file vectors, Datasets, and arrow_dplyr_query", {
  sim <- make_arrow_sim(20, seed = 102)
  pq_dir <- file.path(tempdir(), paste0("ezbakr_inputs_", Sys.getpid()))
  unlink(pq_dir, recursive = TRUE)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  WriteCBArrow(sim$cB, pq_dir)

  # 1. Direct directory path
  ez_dir <- EZbakRArrowData(pq_dir, sim$metadf)
  expect_s3_class(ez_dir, "EZbakRArrowData")
  expect_true(ez_dir$metadata$arrow_spec$partitioned_by_sample)

  # 2. Vector of Hive-partitioned file paths (both full and subset preserve sample partition column)
  pq_files <- list.files(pq_dir, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  ez_files <- EZbakRArrowData(pq_files, sim$metadf)
  expect_s3_class(ez_files, "EZbakRArrowData")
  expect_true("sample" %in% names(arrow::schema(ez_files$cBds)))

  sub_samps <- sim$metadf$sample[1:2]
  pq_sub_files <- pq_files[grepl(paste0("sample=(", paste(sub_samps, collapse = "|"), ")/"),
                                 gsub("\\\\", "/", pq_files))]
  ez_sub_files <- EZbakRArrowData(pq_sub_files, sim$metadf[1:2, ])
  expect_s3_class(ez_sub_files, "EZbakRArrowData")
  expect_setequal(ez_sub_files$metadata$arrow_spec$samples, sub_samps)

  # 3. Open arrow::Dataset
  ds <- arrow::open_dataset(pq_dir)
  ez_ds <- EZbakRArrowData(ds, sim$metadf)
  expect_s3_class(ez_ds, "EZbakRArrowData")

  # 4. Filtered arrow_dplyr_query emits helpful note about filter persistence
  q_filt <- ds %>% dplyr::filter(n >= 1)
  expect_message(
    ez_q <- EZbakRArrowData(q_filt, sim$metadf),
    "in-memory `arrow_dplyr_query` filter"
  )
  expect_s3_class(ez_q, "EZbakRArrowData")
  expect_true(ez_q$metadata$arrow_spec$has_unserializable_filter)

  # Non-existent path throws clean error class
  expect_error(
    EZbakRArrowData(file.path(tempdir(), "definitely_missing_dir_xyz"), sim$metadf),
    class = "arrow_path_not_found"
  )
})


test_that("Validator handles special sample names, unpartitioned datasets, and float count columns", {
  # 1. Hyphenated, dotted, and space-containing sample names in Hive partitions
  cB_special <- data.frame(
    sample = c("WT-1", "WT.2", "KO-1", "KO 2"),
    TC = c(0L, 2L, 1L, 3L),
    nT = c(20L, 20L, 20L, 20L),
    feature = c("Gene1", "Gene1", "Gene1", "Gene1"),
    n = c(5L, 6L, 7L, 8L),
    stringsAsFactors = FALSE
  )
  meta_special <- data.frame(
    sample = c("WT-1", "WT.2", "KO-1", "KO 2"),
    tl = c(2, 2, 2, 2),
    condition = c("WT", "WT", "KO", "KO"),
    stringsAsFactors = FALSE
  )

  pq_special <- file.path(tempdir(), paste0("ezbakr_special_", Sys.getpid()))
  unlink(pq_special, recursive = TRUE)
  on.exit(unlink(pq_special, recursive = TRUE), add = TRUE)
  WriteCBArrow(cB_special, pq_special)

  expect_silent(ez_sp <- EZbakRArrowData(pq_special, meta_special))
  expect_setequal(ez_sp$metadata$arrow_spec$samples, meta_special$sample)

  # 2. Unpartitioned dataset emits informative performance tip & works
  pq_unpart <- file.path(tempdir(), paste0("ezbakr_unpart_", Sys.getpid()))
  unlink(pq_unpart, recursive = TRUE)
  on.exit(unlink(pq_unpart, recursive = TRUE), add = TRUE)
  WriteCBArrow(cB_special, pq_unpart, partitioning = NULL)

  expect_message(
    ez_unpart <- EZbakRArrowData(pq_unpart, meta_special),
    "not Hive-partitioned"
  )
  expect_false(ez_unpart$metadata$arrow_spec$partitioned_by_sample)

  # 3. Whole-number float64 columns pass; fractional/non-positive values fail
  pq_flt_ok <- file.path(tempdir(), paste0("ezbakr_flt_ok_", Sys.getpid()))
  unlink(pq_flt_ok, recursive = TRUE)
  on.exit(unlink(pq_flt_ok, recursive = TRUE), add = TRUE)
  dir.create(pq_flt_ok, recursive = TRUE)

  cB_flt_ok <- transform(cB_special,
                         TC = as.double(TC),
                         nT = as.double(nT),
                         n = as.double(n))
  arrow::write_dataset(cB_flt_ok, pq_flt_ok, format = "parquet", partitioning = "sample")
  expect_silent(EZbakRArrowData(pq_flt_ok, meta_special))

  # Fractional n -> cB_n_int
  pq_flt_bad_n <- file.path(tempdir(), paste0("ezbakr_flt_badn_", Sys.getpid()))
  unlink(pq_flt_bad_n, recursive = TRUE)
  on.exit(unlink(pq_flt_bad_n, recursive = TRUE), add = TRUE)
  cB_bad_n <- cB_flt_ok
  cB_bad_n$n[1] <- 2.5
  arrow::write_dataset(cB_bad_n, pq_flt_bad_n, format = "parquet", partitioning = "sample")
  expect_error(EZbakRArrowData(pq_flt_bad_n, meta_special), class = "cB_n_int")

  # Fractional TC -> mut_not_int
  pq_flt_bad_tc <- file.path(tempdir(), paste0("ezbakr_flt_badtc_", Sys.getpid()))
  unlink(pq_flt_bad_tc, recursive = TRUE)
  on.exit(unlink(pq_flt_bad_tc, recursive = TRUE), add = TRUE)
  cB_bad_tc <- cB_flt_ok
  cB_bad_tc$TC[1] <- 0.5
  arrow::write_dataset(cB_bad_tc, pq_flt_bad_tc, format = "parquet", partitioning = "sample")
  expect_error(EZbakRArrowData(pq_flt_bad_tc, meta_special), class = "mut_not_int")

  # Missing feature columns -> cBds_no_features
  pq_nofeat <- file.path(tempdir(), paste0("ezbakr_nofeat_", Sys.getpid()))
  unlink(pq_nofeat, recursive = TRUE)
  on.exit(unlink(pq_nofeat, recursive = TRUE), add = TRUE)
  cB_nofeat <- cB_special[, c("sample", "TC", "nT", "n")]
  arrow::write_dataset(cB_nofeat, pq_nofeat, format = "parquet", partitioning = "sample")
  expect_error(EZbakRArrowData(pq_nofeat, meta_special), class = "cBds_no_features")
})


test_that("saveRDS and readRDS transparently self-heal EZbakRArrowData across sessions", {
  sim <- make_arrow_sim(25, seed = 103)
  pq_dir <- file.path(tempdir(), paste0("ezbakr_rds_heal_", Sys.getpid()))
  unlink(pq_dir, recursive = TRUE)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  WriteCBArrow(sim$cB, pq_dir)

  ez_arrow <- EZbakRArrowData(pq_dir, sim$metadf)
  rds_path <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_path), add = TRUE)
  saveRDS(ez_arrow, rds_path)

  # Reload from RDS: raw C++ externalptr in cBds and cache$ds are dead (0x0)
  ez_reloaded <- readRDS(rds_path)
  expect_false(EZbakR:::arrow_is_live(ez_reloaded$cBds))
  expect_false(EZbakR:::arrow_is_live(ez_reloaded$metadata$arrow_spec$cache$ds))

  # Non-mutating print() auto-reopens and caches live handle inside spec$cache$ds!
  printed <- capture.output(print(ez_reloaded))
  expect_true(any(grepl("Arrow backend: connected", printed)))
  expect_true(EZbakR:::arrow_is_live(ez_reloaded$metadata$arrow_spec$cache$ds))

  # Fresh reload to test EstimateFractions directly on dead externalptr without prior print()
  ez_reloaded2 <- readRDS(rds_path)
  expect_false(EZbakR:::arrow_is_live(ez_reloaded2$cBds))

  ez_reloaded2 <- EstimateFractions(ez_reloaded2, overwrite = TRUE)
  expect_s3_class(ez_reloaded2, "EZbakRFractions")
  expect_true(EZbakR:::arrow_is_live(ez_reloaded2$cBds))

  # Compare against in-memory EZbakRData
  ez_mem <- EZbakRData(sim$cB, sim$metadf)
  ez_mem <- EstimateFractions(ez_mem)

  fn_arrow <- EZget(ez_reloaded2, type = "fractions") %>%
    dplyr::arrange(sample, feature)
  fn_mem <- EZget(ez_mem, type = "fractions") %>%
    dplyr::arrange(sample, feature)

  expect_equal(fn_arrow$fraction_highTC, fn_mem$fraction_highTC, tolerance = 1e-4)
  expect_equal(fn_arrow$n, fn_mem$n)
})


test_that("Offline mode runs post-fraction analyses cleanly and ReconnectArrow restores moved datasets", {
  sim <- make_arrow_sim(30, seed = 104)
  pq_orig <- file.path(tempdir(), paste0("ezbakr_orig_", Sys.getpid()))
  pq_moved <- file.path(tempdir(), paste0("ezbakr_moved_", Sys.getpid()))
  unlink(pq_orig, recursive = TRUE)
  unlink(pq_moved, recursive = TRUE)
  on.exit({
    unlink(pq_orig, recursive = TRUE)
    unlink(pq_moved, recursive = TRUE)
  }, add = TRUE)

  WriteCBArrow(sim$cB, pq_orig)
  ez_arrow <- EZbakRArrowData(pq_orig, sim$metadf)
  ez_arrow <- EstimateFractions(ez_arrow)

  rds_path <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_path), add = TRUE)
  saveRDS(ez_arrow, rds_path)

  # Move the Parquet folder to simulate moving to another path / machine
  expect_true(file.rename(pq_orig, pq_moved))
  expect_false(dir.exists(pq_orig))

  ez_offline <- readRDS(rds_path)

  # 1. print() works cleanly in offline mode
  out_off <- capture.output(print(ez_offline))
  expect_true(any(grepl("Arrow backend: offline", out_off)))

  # 2. Downstream post-cB workflow succeeds 100% offline!
  ez_offline <- EstimateKinetics(ez_offline)
  ez_offline <- AverageAndRegularize(ez_offline, formula_mean = ~ treatment - 1)
  ez_offline <- CompareParameters(ez_offline,
                                  parameter = "log_kdeg",
                                  design_factor = "treatment",
                                  reference = "treatment1",
                                  experimental = "treatment2")
  expect_true(length(ez_offline$comparisons) == 1L)

  # 3. EZQC on EZbakRFractions skips only raw mutation scan with a clear message
  expect_message(
    qc_off <- EZQC(ez_offline),
    "Arrow dataset is offline; skipping raw mutation rate scan"
  )
  expect_null(qc_off$Raw_mutrates)
  expect_false(is.null(qc_off$Fraction_labeled_dist))

  # 4. Raw-data methods throw ezbakr_arrow_disconnected with ReconnectArrow guidance
  expect_error(
    EstimateMutRates(ez_offline),
    class = "ezbakr_arrow_disconnected"
  )
  expect_error(
    CollectArrow(ez_offline),
    class = "ezbakr_arrow_disconnected"
  )

  # 5. ReconnectArrow() re-binds to moved folder while preserving all downstream analyses!
  ez_reconn <- ReconnectArrow(ez_offline, pq_moved)
  expect_s3_class(ez_reconn, "EZbakRKinetics")
  expect_s3_class(ez_reconn, "EZbakRFractions")
  expect_s3_class(ez_reconn, "EZbakRArrowData")
  expect_equal(ez_reconn$comparisons, ez_offline$comparisons)
  expect_equal(ez_reconn$kinetics, ez_offline$kinetics)

  out_on <- capture.output(print(ez_reconn))
  expect_true(any(grepl("Arrow backend: connected", out_on)))

  # Raw mutation rate QC now works again!
  qc_on <- suppressMessages(EZQC(ez_reconn))
  expect_false(is.null(qc_on$Raw_mutrates))

  # ReconnectArrow rejects dataset missing required samples
  pq_partial <- file.path(tempdir(), paste0("ezbakr_partial_", Sys.getpid()))
  unlink(pq_partial, recursive = TRUE)
  on.exit(unlink(pq_partial, recursive = TRUE), add = TRUE)
  WriteCBArrow(sim$cB[sim$cB$sample == sim$metadf$sample[1], ], pq_partial)

  expect_error(
    ReconnectArrow(ez_offline, pq_partial),
    class = "cB_metadf_samples"
  )
})


test_that("get_normalized_read_counts.EZbakRArrowData matches in-memory EZbakRData", {
  sim <- make_arrow_sim(25, seed = 105)
  pq_dir <- file.path(tempdir(), paste0("ezbakr_normrc_", Sys.getpid()))
  unlink(pq_dir, recursive = TRUE)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  WriteCBArrow(sim$cB, pq_dir)

  # Subset metadf to verify sample pushdown filtering matches EZbakRData on those samples
  sub_meta <- sim$metadf[1:3, ]
  ez_arrow <- EZbakRArrowData(pq_dir, sub_meta)
  ez_mem <- EZbakRData(sim$cB[sim$cB$sample %in% sub_meta$sample, ], sub_meta)

  rc_arrow <- get_normalized_read_counts(ez_arrow, features_to_analyze = "feature") %>%
    dplyr::arrange(sample, feature)
  rc_mem <- get_normalized_read_counts(ez_mem, features_to_analyze = "feature") %>%
    dplyr::arrange(sample, feature)

  expect_equal(rc_arrow$n, rc_mem$n)
  expect_equal(rc_arrow$normalized_reads, rc_mem$normalized_reads, tolerance = 1e-8)
  expect_equal(rc_arrow$scale_factor, rc_mem$scale_factor, tolerance = 1e-8)
})


test_that("CollectArrow materializes full and subsetted EZbakRData objects", {
  sim <- make_arrow_sim(25, seed = 106)
  pq_dir <- file.path(tempdir(), paste0("ezbakr_collect_", Sys.getpid()))
  unlink(pq_dir, recursive = TRUE)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  WriteCBArrow(sim$cB, pq_dir)

  ez_arrow <- EZbakRArrowData(pq_dir, sim$metadf)
  ez_arrow <- EstimateFractions(ez_arrow)
  ez_arrow <- EstimateKinetics(ez_arrow)

  # Full materialization
  ez_full <- CollectArrow(ez_arrow)
  expect_s3_class(ez_full, "EZbakRData")
  expect_false(inherits(ez_full, "EZbakRArrowData"))
  expect_s3_class(ez_full, "EZbakRKinetics")
  expect_equal(nrow(ez_full$cB), nrow(sim$cB))
  expect_null(ez_full$metadata$arrow_spec)

  # Filtered slice by samples and features
  keep_samps <- sim$metadf$sample[1:2]
  keep_feats <- unique(sim$cB$feature)[1:5]
  ez_sub <- CollectArrow(ez_arrow, samples = keep_samps, features = keep_feats)

  expect_s3_class(ez_sub, "EZbakRData")
  expect_true(all(ez_sub$cB$sample %in% keep_samps))
  expect_true(all(ez_sub$cB$feature %in% keep_feats))
  expect_equal(ez_sub$metadf$sample, keep_samps)
  expect_true(all(EZget(ez_sub, type = "fractions")$sample %in% keep_samps))
  expect_true(all(EZget(ez_sub, type = "fractions")$feature %in% keep_feats))
  expect_true(all(EZget(ez_sub, type = "kinetics")$sample %in% keep_samps))
  expect_true(all(EZget(ez_sub, type = "kinetics")$feature %in% keep_feats))

  # copy_analyses = FALSE returns bare EZbakRData
  ez_bare <- CollectArrow(ez_arrow, samples = keep_samps, copy_analyses = FALSE)
  expect_identical(class(ez_bare), "EZbakRData")
  expect_null(ez_bare$fractions)
})


test_that("EstimateKinetics error message and EZconcat handle EZbakRArrowData properly", {
  sim <- make_arrow_sim(15, seed = 107)
  pq_dir <- file.path(tempdir(), paste0("ezbakr_misc_", Sys.getpid()))
  unlink(pq_dir, recursive = TRUE)
  on.exit(unlink(pq_dir, recursive = TRUE), add = TRUE)
  WriteCBArrow(sim$cB, pq_dir)

  ez_arrow <- EZbakRArrowData(pq_dir, sim$metadf)
  expect_error(
    EstimateKinetics(ez_arrow),
    "EstimateFractions"
  )

  ez1 <- EstimateFractions(ez_arrow)
  ez2 <- EstimateFractions(ez_arrow)
  ez_cat <- EZbakR:::EZconcat(ez1, ez2)
  expect_s3_class(ez_cat, "EZbakRArrowData")
  expect_false(is.null(ez_cat$metadata$arrow_spec))
})
