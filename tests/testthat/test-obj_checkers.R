
test_that("check_null errors on NULL and passes otherwise", {
  expect_error(check_null(NULL), "null object")
  expect_silent(check_null(tmt_qf))
})

test_that("check_q errors on non-QFeatures and passes on QFeatures", {
  expect_error(check_q(data.frame()), "must be a QFeatures object")
  expect_error(check_q(NULL), "null object")
  expect_silent(check_q(tmt_qf))
})

test_that("check_se errors on non-SummarizedExperiment and passes otherwise", {
  expect_error(check_se(tmt_qf), "must be a SummarizedExperiment object")
  expect_silent(check_se(tmt_qf[["psms_raw"]]))
})

test_that("check_se_exists errors when the named experiment is missing", {
  expect_error(check_se_exists(tmt_qf, "nonexistent"), "null object")
  expect_silent(check_se_exists(tmt_qf, "psms_raw"))
})

test_that("check_colData_col errors when the column is absent", {
  expect_error(check_colData_col(tmt_qf, "nonexistent"), "not in the colData")
})

test_that("check_se_psm/peptide/protein error on non-SummarizedExperiment and pass otherwise", {
  expect_error(check_se_psm(tmt_qf), "must be a SummarizedExperiment object")
  expect_silent(check_se_psm(tmt_qf[["psms_raw"]]))

  expect_error(check_se_peptide(tmt_qf), "must be a SummarizedExperiment object")
  expect_silent(check_se_peptide(tmt_qf[["psms_raw"]]))

  expect_error(check_se_protein(tmt_qf), "must be a SummarizedExperiment object")
  expect_silent(check_se_protein(tmt_qf[["psms_raw"]]))
})
