make_peptide_se <- function() {
  row_data <- S4Vectors::DataFrame(
    Master.Protein.Accessions = c("P1", "P2", "P3", "P4")
  )
  assay_mat <- matrix(1, nrow = 4, ncol = 1, dimnames = list(NULL, "s1"))
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)
}

write_protein_fdr_file <- function() {
  path <- tempfile(fileext = ".txt")
  df <- data.frame(
    Accession = c("P1", "P2", "P3"),
    Protein.FDR.Confidence.Combined = c("High", "Low", "High"),
    stringsAsFactors = FALSE
  )
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

test_that("filter_by_protein_fdr retains only peptides from High-confidence proteins", {
  se <- make_peptide_se()
  fdr_path <- write_protein_fdr_file()

  out <- suppressMessages(filter_by_protein_fdr(se, fdr_path))

  # P1: High (kept), P2: Low (dropped), P3: High (kept), P4: not in FDR file -> NA (dropped)
  expect_setequal(SummarizedExperiment::rowData(out)$Master.Protein.Accessions, c("P1", "P3"))
})

test_that("filter_by_protein_fdr always retains proteins in retain_proteins", {
  se <- make_peptide_se()
  fdr_path <- write_protein_fdr_file()

  out <- suppressMessages(filter_by_protein_fdr(se, fdr_path, retain_proteins = c("P4")))

  expect_setequal(SummarizedExperiment::rowData(out)$Master.Protein.Accessions, c("P1", "P3", "P4"))
})
