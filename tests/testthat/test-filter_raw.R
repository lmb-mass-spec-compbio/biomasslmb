load(test_path("..", "..", "data", "tmt_qf.rda"))
psm <- tmt_qf[["psms_raw"]]

test_that("remove_no_master drops features with an empty master protein", {
  n_empty <- sum(SummarizedExperiment::rowData(psm)[["Master.Protein.Accessions"]] == "")
  expect_gt(n_empty, 0)

  out <- suppressMessages(remove_no_master(psm, "Master.Protein.Accessions"))
  expect_equal(nrow(out), nrow(psm) - n_empty)
  expect_false(any(SummarizedExperiment::rowData(out)[["Master.Protein.Accessions"]] == ""))
})

test_that("remove_non_unique_master_protein keeps only Number.of.Protein.Groups == 1", {
  out <- suppressMessages(remove_non_unique_master_protein(psm, "Master.Protein.Accessions"))
  expect_true(all(SummarizedExperiment::rowData(out)[["Number.of.Protein.Groups"]] == 1))
  expect_lt(nrow(out), nrow(psm))
})

test_that("remove_no_quant_assay drops features with no finite quantification", {
  n_no_quant <- sum(rowSums(is.finite(SummarizedExperiment::assay(psm))) == 0)
  expect_gt(n_no_quant, 0)

  out <- suppressMessages(remove_no_quant_assay(psm, "Master.Protein.Accessions"))
  expect_equal(nrow(out), nrow(psm) - n_no_quant)
})

test_that("filter_features_pd_dda applies all filtering steps together", {
  out <- suppressMessages(filter_features_pd_dda(
    psm,
    filter_contaminant = FALSE,
    unique_master = TRUE,
    remove_no_quant = TRUE
  ))

  expect_s4_class(out, "SummarizedExperiment")
  expect_lt(nrow(out), nrow(psm))
  expect_false(any(SummarizedExperiment::rowData(out)[["Master.Protein.Accessions"]] == ""))
  expect_true(all(SummarizedExperiment::rowData(out)[["Number.of.Protein.Groups"]] == 1))
  expect_true(all(rowSums(is.finite(SummarizedExperiment::assay(out))) > 0))
})
