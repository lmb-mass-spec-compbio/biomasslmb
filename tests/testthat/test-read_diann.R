# diann_report_subset.tsv is the first 3000 rows (6 runs) of the bundled
# inst/extdata/diann_report.tsv, trimmed to the columns readDIANNFilterQJoin
# actually uses, to keep the fixture small.
fixture_path <- test_path("read_diann", "diann_report_subset.tsv")

test_that("readDIANNFilterQJoin filters on Q-value thresholds and joins runs into one assay", {
  out <- suppressWarnings(suppressMessages(readDIANNFilterQJoin(fixture_path)))

  expect_s4_class(out, "QFeatures")
  expect_true(nrow(out[["peptides_fdr_cntrl"]]) > 0)
  expect_equal(ncol(out[["peptides_fdr_cntrl"]]), 6)
})

test_that("readDIANNFilterQJoin returns fewer rows with a stricter Q-value threshold", {
  lenient <- suppressWarnings(suppressMessages(readDIANNFilterQJoin(fixture_path, run_precursor_q = 1)))
  strict <- suppressWarnings(suppressMessages(readDIANNFilterQJoin(fixture_path, run_precursor_q = 1e-6)))

  expect_lt(nrow(strict[["peptides_fdr_cntrl"]]), nrow(lenient[["peptides_fdr_cntrl"]]))
})

test_that("readDIANNFilterQJoin can return per-sample quantification alongside the joined assay", {
  out <- suppressMessages(readDIANNFilterQJoin(fixture_path, return_sep_quant = TRUE))

  expect_true("peptides_fdr_cntrl" %in% names(out))
  expect_gt(length(out), 1) # separate per-run assays plus the joined one
})
