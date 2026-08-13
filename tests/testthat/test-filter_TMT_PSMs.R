psm <- tmt_qf[["psms_raw"]]

test_that("update_average_sn adds a row-mean signal:noise column ignoring NA", {
  out <- update_average_sn(psm)

  expect_equal(
    SummarizedExperiment::rowData(out)$Average.Reporter.SN,
    rowMeans(SummarizedExperiment::assay(psm), na.rm = TRUE)
  )
})

test_that("update_average_sn writes to a custom column name", {
  out <- update_average_sn(psm, sn_col = "my_sn")
  expect_true("my_sn" %in% colnames(SummarizedExperiment::rowData(out)))
})

test_that("filter_TMT_PSMs removes Quan.Info-flagged and ambiguous PSMs and no-quant rows", {
  se <- psm
  rd <- SummarizedExperiment::rowData(se)
  rd$Quan.Info[1:3] <- "No Quan Values"
  rd$PSM.Ambiguity[4:6] <- "Rejected"
  SummarizedExperiment::rowData(se) <- rd

  out <- suppressMessages(filter_TMT_PSMs(se, from_PD = TRUE, verbose = FALSE))

  expect_false(any(SummarizedExperiment::rowData(out)$Quan.Info != ""))
  expect_true(all(SummarizedExperiment::rowData(out)$PSM.Ambiguity %in% c("Selected", "Unambiguous")))
  expect_true(all(rowSums(is.finite(SummarizedExperiment::assay(out))) > 0))
  expect_lt(nrow(out), nrow(se))
})

test_that("filter_TMT_PSMs applies interference, S:N and SPS-MM thresholds when set", {
  se <- update_average_sn(psm)
  se <- se[is.finite(SummarizedExperiment::rowData(se)$Isolation.Interference.in.Percent) &
             is.finite(SummarizedExperiment::rowData(se)$SPS.Mass.Matches.in.Percent), ]

  inter_thresh <- stats::median(SummarizedExperiment::rowData(se)$Isolation.Interference.in.Percent, na.rm = TRUE)
  sn_thresh <- stats::median(SummarizedExperiment::rowData(se)$Average.Reporter.SN, na.rm = TRUE)
  spsmm_thresh <- stats::median(SummarizedExperiment::rowData(se)$SPS.Mass.Matches.in.Percent, na.rm = TRUE)

  out <- suppressMessages(filter_TMT_PSMs(
    se,
    inter_thresh = inter_thresh,
    sn_thresh = sn_thresh,
    spsmm_thresh = spsmm_thresh,
    from_PD = FALSE,
    verbose = FALSE
  ))

  expect_true(all(SummarizedExperiment::rowData(out)$Isolation.Interference.in.Percent <= inter_thresh))
  expect_true(all(SummarizedExperiment::rowData(out)$Average.Reporter.SN >= sn_thresh))
  expect_true(all(SummarizedExperiment::rowData(out)$SPS.Mass.Matches.in.Percent >= spsmm_thresh))
})

test_that("filter_TMT_PSMs skips PD-specific filtering when from_PD = FALSE", {
  se <- psm
  rd <- SummarizedExperiment::rowData(se)
  rd$Quan.Info[1:3] <- "No Quan Values"
  SummarizedExperiment::rowData(se) <- rd

  out <- suppressMessages(filter_TMT_PSMs(se, from_PD = FALSE, verbose = FALSE))

  expect_true(any(SummarizedExperiment::rowData(out)$Quan.Info != ""))
})
