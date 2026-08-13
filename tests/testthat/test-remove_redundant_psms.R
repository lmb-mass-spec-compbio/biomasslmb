write_pep_file <- function(df) {
  path <- tempfile(fileext = ".txt")
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

test_that("update_peptide_assignments replaces psm-level assignments with peptide-level ones", {
  row_data <- S4Vectors::DataFrame(
    Sequence = c("PEPTIDEA", "PEPTIDEB"),
    Protein.Accessions = c("P1_old", "P2_old"),
    Master.Protein.Accessions = c("P1_old", "P2_old"),
    row.names = c("f1", "f2")
  )
  assay_mat <- matrix(1, nrow = 2, ncol = 1, dimnames = list(c("f1", "f2"), "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  pep_inf <- write_pep_file(data.frame(
    Sequence = c("PEPTIDEA", "PEPTIDEB"),
    Protein.Accessions = c("P1_new", "P2_new"),
    Master.Protein.Accessions = c("P1_new", "P2_new"),
    stringsAsFactors = FALSE
  ))

  out <- update_peptide_assignments(se, pep_inf)

  expect_equal(SummarizedExperiment::rowData(out)[c("f1", "f2"), "Master.Protein.Accessions"],
               c("P1_new", "P2_new"))
})

test_that("update_peptide_assignments drops psms with a sequence not in the peptide file", {
  row_data <- S4Vectors::DataFrame(
    Sequence = c("PEPTIDEA", "PEPTIDEB"),
    Master.Protein.Accessions = c("P1_old", "P2_old"),
    row.names = c("f1", "f2")
  )
  assay_mat <- matrix(1, nrow = 2, ncol = 1, dimnames = list(c("f1", "f2"), "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  pep_inf <- write_pep_file(data.frame(
    Sequence = "PEPTIDEA",
    Master.Protein.Accessions = "P1_new",
    stringsAsFactors = FALSE
  ))

  out <- update_peptide_assignments(se, pep_inf)

  expect_equal(rownames(out), "f1")
})

test_that("update_peptide_assignments errors when the peptide file has non-unique assignments", {
  row_data <- S4Vectors::DataFrame(
    Sequence = "PEPTIDEA",
    Master.Protein.Accessions = "P1_old",
    row.names = "f1"
  )
  assay_mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("f1", "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  pep_inf <- write_pep_file(data.frame(
    Sequence = c("PEPTIDEA", "PEPTIDEA"),
    Master.Protein.Accessions = c("P1_new", "P1_other"),
    stringsAsFactors = FALSE
  ))

  expect_error(update_peptide_assignments(se, pep_inf), "not unique")
})

test_that("update_peptide_assignments works on a bare data.frame input", {
  row_data <- data.frame(
    Sequence = c("PEPTIDEA", "PEPTIDEB"),
    Master.Protein.Accessions = c("P1_old", "P2_old"),
    row.names = c("f1", "f2"),
    stringsAsFactors = FALSE
  )

  pep_inf <- write_pep_file(data.frame(
    Sequence = c("PEPTIDEA", "PEPTIDEB"),
    Master.Protein.Accessions = c("P1_new", "P2_new"),
    stringsAsFactors = FALSE
  ))

  out <- update_peptide_assignments(row_data, pep_inf)

  expect_true(is.data.frame(out))
  expect_equal(out[c("f1", "f2"), "Master.Protein.Accessions"], c("P1_new", "P2_new"))
})

make_psm_se <- function() {
  row_data <- S4Vectors::DataFrame(
    Sequence = c("PEP1", "PEP1", "PEP2"),
    Protein.Accessions = c("P1", "P1", "P2"),
    First.Scan = c(100, 100, 200),
    Last.Scan = c(101, 101, 201),
    Master.Scans = c("100", "100", "200"),
    File.ID = c("F1", "F1", "F1")
  )
  assay_mat <- matrix(c(1, 1, 2, 1, 1, 2), nrow = 3, ncol = 2, dimnames = list(NULL, c("s1", "s2")))
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)
}

test_that("remove_redundant_psm_quant removes exact-duplicate PSM quantification rows", {
  se <- make_psm_se()
  out <- suppressMessages(remove_redundant_psm_quant(se))

  expect_equal(nrow(out), 2)
})

test_that("remove_redundant_psm_quant validates required spectra ID columns are present", {
  se <- make_psm_se()
  SummarizedExperiment::rowData(se)$First.Scan <- NULL

  expect_error(suppressMessages(remove_redundant_psm_quant(se)), "First.Scan")
})

test_that("remove_redundant_psm_quant with validate=TRUE errors on inconsistent Sequence/Protein.Accessions", {
  se <- make_psm_se()
  rd <- SummarizedExperiment::rowData(se)
  rd$Protein.Accessions[2] <- "P_different"
  SummarizedExperiment::rowData(se) <- rd

  expect_error(suppressMessages(remove_redundant_psm_quant(se, validate = TRUE)), "update_peptide_assignments")
})

test_that("remove_redundant_psm_quant with validate=FALSE skips the consistency check", {
  se <- make_psm_se()
  rd <- SummarizedExperiment::rowData(se)
  rd$Protein.Accessions[2] <- "P_different"
  rd$First.Scan[2] <- 999 # no longer a true duplicate of row 1
  SummarizedExperiment::rowData(se) <- rd

  out <- suppressMessages(remove_redundant_psm_quant(se, validate = FALSE))
  expect_equal(nrow(out), 3)
})
