### CODE THAT SHOULD WORK

test_that("TILAC EZbakRDataobject creation works", {
  cB_TIL <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                       TC = c(1, 1, 1),
                       nT = c(20, 20, 12),
                       GA = c(0, 2, 1),
                       nG = c(12, 5, 15),
                       XF = "a",
                       GF = "b",
                       n = c(1, 2, 3))

  metadf <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                       tl_TC = c(1, 1, 0),
                       tl_GA = c(1, 1, 0),
                       treatment = c("WT_1", "WT_1", "WT_1"))


  expect_silent(EZbakRData(cB_TIL, metadf))
})


test_that("STL EZbakRDataobject creation works", {
  cB_STL <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                       TC = c(1, 1, 1),
                       nT = c(20, 20, 12),
                       TSS = "TSS1",
                       n = c(1, 2, 3))

  metadf <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                       tl = c(1, 1, 0),
                       treatment = c("WT", "WT", "WT"))


  expect_silent(EZbakRData(cB_STL, metadf))
})


test_that("Pulse-chase EZbakRDataobject creation works", {
  cB_pc <- data.table(sample = c("WT_pulse", "WT_chase", "WT_ctl"),
                      TC = c(1, 1, 1),
                      nT = c(20, 20, 12),
                      TSS = "GeneA",
                      n = c(1, 2, 3))

  metadf <- data.table(sample = c("WT_pulse", "WT_chase", "WT_ctl"),
                       tpulse = c(1, 1, 0),
                       tchase = c(0, 1, 0),
                       treatment = c("WT", "WT", "WT"))


  expect_silent(EZbakRData(cB_pc, metadf))
})


test_that("Dual label pulse-chase EZbakRDataobject creation works", {
  cB_pc <- data.table(sample = c("WT_pulse", "WT_chase", "WT_ctl"),
                      TC = c(1, 1, 1),
                      nT = c(20, 20, 12),
                      GA = c(1, 2, 1),
                      nG = c(10, 1, 2),
                      XF = "GeneA",
                      n = c(1, 2, 3))

  metadf <- data.table(sample = c("WT_pulse", "WT_chase", "WT_ctl"),
                       tpulse_TC = c(2, 1, 0),
                       tchase_TC = c(0, 1, 0),
                       tpulse_GA = c(1, 1, 0),
                       tchase_GA = c(0, 1, 0),
                       treatment = c("WT", "WT", "WT"))


  expect_silent(EZbakRData(cB_pc, metadf))
})


### CODE THAT SHOULD THROW ERRORS

test_that("TILAC EZbakRDataobject throws error if nucleotide counts aren't
          whole numbers.", {
  cB_TIL <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                       TC = c(1, 1, 1),
                       nT = c(20, 20, 12.5),
                       GA = c(0, 2, 1),
                       nG = c(12, 5, 15),
                       XF = "a",
                       GF = "b",
                       n = c(1, 2, 3))

  metadf <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                       tl_TC = c(1, 1, 0),
                       tl_GA = c(1, 1, 0),
                       treatment = c("WT_1", "WT_1", "WT_1"))


  expect_error(EZbakRData(cB_TIL, metadf),
               class = "basecounts_not_pos_whole")
})

test_that("TILAC EZbakRDataobject throws error if n isn't always
          whole numbers.", {
            cB_TIL <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                                 TC = c(1, 1, 1),
                                 nT = c(20, 20, 12),
                                 GA = c(0, 2, 1),
                                 nG = c(12, 5, 15),
                                 XF = "a",
                                 GF = "b",
                                 n = c(1.5, 2, 3))

            metadf <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                                 tl_TC = c(1, 1, 0),
                                 tl_GA = c(1, 1, 0),
                                 treatment = c("WT_1", "WT_1", "WT_1"))


            expect_error(EZbakRData(cB_TIL, metadf),
                         class = "cB_n_pos_whole")
})


test_that("TILAC EZbakRDataobject throws error if missing a tl column in metadf", {
            cB_TIL <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                                 TC = c(1, 1, 1),
                                 nT = c(20, 20, 12),
                                 GA = c(0, 2, 1),
                                 nG = c(12, 5, 15),
                                 XF = "a",
                                 GF = "b",
                                 n = c(1.5, 2, 3))

            metadf <- data.table(sample = c("WT_1", "WT_2", "WT_ctl"),
                                 tl_GA = c(1, 1, 0),
                                 treatment = c("WT_1", "WT_1", "WT_1"))


            expect_error(EZbakRData(cB_TIL, metadf), class = "metadf_tl_multimut")
          })



