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

test_that("filter_by_protein_fdr honours protein_FDR_col", {
  se <- make_peptide_se()

  path <- tempfile(fileext = ".txt")
  df <- data.frame(
    Accession = c("P1", "P2", "P3"),
    # the confidence lives in a differently named column, as it does when the
    # protein-group rather than the protein FDR is the one of interest
    Protein.Group.FDR.Confidence.Combined = c("High", "Low", "High"),
    stringsAsFactors = FALSE
  )
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)

  out <- suppressMessages(filter_by_protein_fdr(
    se, path, protein_FDR_col = "Protein.Group.FDR.Confidence.Combined"))

  expect_setequal(SummarizedExperiment::rowData(out)$Master.Protein.Accessions,
                  c("P1", "P3"))
})

test_that("filter_by_protein_fdr always retains proteins in retain_proteins", {
  se <- make_peptide_se()
  fdr_path <- write_protein_fdr_file()

  out <- suppressMessages(filter_by_protein_fdr(se, fdr_path, retain_proteins = c("P4")))

  expect_setequal(SummarizedExperiment::rowData(out)$Master.Protein.Accessions, c("P1", "P3", "P4"))
})

make_group_se <- function() {
  row_data <- S4Vectors::DataFrame(
    Master.Protein.Accessions = c(
      "P1", "P2", "P4",  # single: High, Low, absent
      "P1; P3",          # every member High
      "P1;P3",           # as above, no space after the separator
      "P1; P2",          # mixed confidence
      "P4; P1",          # one High, one absent
      "P2; P5",          # every member Low
      "P4; P6",          # every member absent
      "", NA)
  )
  assay_mat <- matrix(1, nrow = nrow(row_data), ncol = 1, dimnames = list(NULL, "s1"))
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)
}

write_group_fdr_file <- function(accession = c("P1", "P2", "P3", "P5"),
                                 confidence = c("High", "Low", "High", "Low")) {
  path <- tempfile(fileext = ".txt")
  df <- data.frame(Accession = accession,
                   Protein.FDR.Confidence.Combined = confidence)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

filter_groups <- function(..., fdr_path = write_group_fdr_file()) {
  out <- suppressMessages(filter_by_protein_fdr(make_group_se(), fdr_path, ...))
  SummarizedExperiment::rowData(out)$Master.Protein.Accessions
}

test_that("filter_by_protein_fdr drops multi-accession groups by default", {
  expect_identical(filter_groups(), "P1")
  expect_identical(filter_groups(multi_accession = "remove"), "P1")
})

test_that("filter_by_protein_fdr multi_accession = 'all' needs every member to be High", {
  expect_identical(filter_groups(multi_accession = "all"), c("P1", "P1; P3", "P1;P3"))
})

test_that("filter_by_protein_fdr multi_accession = 'any' needs one member to be High", {
  expect_identical(filter_groups(multi_accession = "any"),
                   c("P1", "P1; P3", "P1;P3", "P1; P2", "P4; P1"))
})

test_that("filter_by_protein_fdr treats single accessions the same under every multi_accession", {
  for (setting in c("remove", "all", "any")) {
    out <- suppressMessages(filter_by_protein_fdr(
      make_peptide_se(), write_protein_fdr_file(), multi_accession = setting))
    expect_identical(SummarizedExperiment::rowData(out)$Master.Protein.Accessions,
                     c("P1", "P3"))
  }
})

test_that("filter_by_protein_fdr honours retain_proteins with multi_accession", {
  expect_identical(
    filter_groups(retain_proteins = c("P4", "P2; P5"), multi_accession = "all"),
    c("P1", "P4", "P1; P3", "P1;P3", "P2; P5"))
})

test_that("filter_by_protein_fdr rejects an invalid multi_accession", {
  expect_error(filter_groups(multi_accession = "some"), "should be one of")
})

test_that("filter_by_protein_fdr errors on an accession with conflicting confidence", {
  conflicting <- write_group_fdr_file(c("P1", "P1", "P2"), c("High", "Low", "High"))
  for (setting in c("remove", "all", "any")) {
    expect_error(filter_groups(fdr_path = conflicting, multi_accession = setting),
                 "more than one value.*P1")
  }

  # exact duplicate rows do not conflict
  repeated <- write_group_fdr_file(c("P1", "P1", "P2"), c("High", "High", "Low"))
  expect_identical(filter_groups(fdr_path = repeated), "P1")
})

test_that("filter_by_protein_fdr removes blank and NA accessions, even if the file has blank rows", {
  blank <- write_group_fdr_file(c("P1", "", "", NA), c("High", "High", "Low", "High"))
  for (setting in c("remove", "all", "any")) {
    out <- filter_groups(fdr_path = blank, multi_accession = setting)
    expect_false(any(out %in% c("", NA)))
  }
})
