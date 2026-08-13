psm <- tmt_qf[["psms_raw"]]

test_that("get_medians returns per-sample column medians ignoring NA", {
  medians <- get_medians(psm)
  expect_equal(medians, robustbase::colMedians(SummarizedExperiment::assay(psm), na.rm = TRUE))
  expect_length(medians, ncol(psm))
})

test_that("center_normalise_to_ref divides by medians scaled to mean 1 on linear scale", {
  se <- psm
  medians <- get_medians(se)

  out <- center_normalise_to_ref(se, medians, center_to_zero = FALSE, on_log_scale = FALSE)

  expected_medians <- medians / mean(medians)
  expected_assay <- t(t(SummarizedExperiment::assay(se)) / expected_medians)
  expect_equal(SummarizedExperiment::assay(out), expected_assay)
})

test_that("center_normalise_to_ref subtracts medians centred on zero on log scale", {
  se <- psm
  medians <- get_medians(se)

  out <- center_normalise_to_ref(se, medians, center_to_zero = TRUE, on_log_scale = TRUE)

  expected_assay <- t(t(SummarizedExperiment::assay(se)) - medians)
  expect_equal(SummarizedExperiment::assay(out), expected_assay)
})

test_that("center_normalise_to_ref centers medians to mean zero when center_to_zero=FALSE and on_log_scale=TRUE", {
  se <- psm
  medians <- get_medians(se)

  out <- center_normalise_to_ref(se, medians, center_to_zero = FALSE, on_log_scale = TRUE)

  expected_medians <- medians - mean(medians)
  expected_assay <- t(t(SummarizedExperiment::assay(se)) - expected_medians)
  expect_equal(SummarizedExperiment::assay(out), expected_assay)
})
