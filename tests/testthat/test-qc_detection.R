suppressPackageStartupMessages(library(QFeatures))
sub_qf <- suppressWarnings(tmt_qf[, , c("psms_raw", "psms_filtered")])

test_that("get_samples_present tallies how many samples each protein was detected in per assay", {
  out <- suppressMessages(get_samples_present(sub_qf, rowVars = "Master.Protein.Accessions"))

  expect_true(is.data.frame(out))
  expect_setequal(colnames(out), c("Master.Protein.Accessions", "psms_raw", "psms_filtered"))
  expect_true(all(out$psms_raw >= 1 & out$psms_raw <= 10))
  # psms_filtered is NA for proteins entirely absent after filtering
  expect_true(all(out$psms_filtered[!is.na(out$psms_filtered)] >= 1 & out$psms_filtered[!is.na(out$psms_filtered)] <= 10))
  expect_true(any(is.na(out$psms_filtered)))

  # cross check one protein by hand against the raw assay
  protein <- out$Master.Protein.Accessions[[1]]
  rd <- SummarizedExperiment::rowData(sub_qf[["psms_raw"]])
  mat <- SummarizedExperiment::assay(sub_qf[["psms_raw"]])
  feature_rows <- rd$Master.Protein.Accessions == protein
  expected_n <- sum(colSums(is.finite(mat[feature_rows, , drop = FALSE])) > 0)
  expect_equal(out$psms_raw[out$Master.Protein.Accessions == protein], expected_n)
})

test_that("get_samples_present applies rename_cols to select and rename experiments", {
  out <- suppressMessages(get_samples_present(
    sub_qf, rowVars = "Master.Protein.Accessions",
    rename_cols = c(raw = "psms_raw")
  ))

  expect_setequal(colnames(out), c("Master.Protein.Accessions", "raw"))
})

test_that("plot_samples_present returns a ggplot object", {
  samples_present <- suppressMessages(get_samples_present(sub_qf, rowVars = "Master.Protein.Accessions"))
  p <- plot_samples_present(samples_present, rowvars = "Master.Protein.Accessions")

  expect_s3_class(p, "ggplot")
})
