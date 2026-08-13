write_proteome_fasta <- function() {
  path <- tempfile(fileext = ".fasta")
  # Trypsin cleaves after K/R: MSAK | PEPTIDEAK | GGVPEPTIDEBR | AA
  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c("sp|P1|PROT1_HUMAN" = "MSAKPEPTIDEAKGGVPEPTIDEBRAA")),
    filepath = path
  )
  path
}

test_that("add_peptide_positions locates a uniquely-occurring peptide within its protein", {
  fasta_path <- write_proteome_fasta()
  row_data <- S4Vectors::DataFrame(
    Master.Protein.Accessions = c("P1", "P1"),
    Sequence = c("PEPTIDEAK", "GGVPEPTIDEBR")
  )
  assay_mat <- matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  out <- suppressWarnings(add_peptide_positions(se, fasta_path))
  rd <- SummarizedExperiment::rowData(out)

  expect_equal(rd$peptide_start, c(5, 14))
  expect_equal(rd$peptide_end, c(13, 25))
  expect_equal(rd$protein_length, c(27, 27))
})

test_that("add_peptide_positions marks a peptide with a protein missing from the fasta as NA", {
  fasta_path <- write_proteome_fasta()
  row_data <- S4Vectors::DataFrame(
    Master.Protein.Accessions = "P2", # not in the fasta
    Sequence = "PEPTIDEAK"
  )
  assay_mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  out <- suppressWarnings(add_peptide_positions(se, fasta_path))
  rd <- SummarizedExperiment::rowData(out)

  expect_true(is.na(rd$peptide_start))
  expect_equal(rd$peptide_position_info, "Protein not in fasta")
})

test_that("add_peptide_positions marks peptides with multiple locations in the protein as NA", {
  fasta_path <- write_proteome_fasta()
  row_data <- S4Vectors::DataFrame(
    Master.Protein.Accessions = "P1",
    Sequence = "A" # occurs twice at the very end (positions 26 and 27)
  )
  assay_mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  out <- suppressWarnings(add_peptide_positions(se, fasta_path))
  rd <- SummarizedExperiment::rowData(out)

  expect_true(is.na(rd$peptide_start))
  expect_equal(rd$peptide_position_info, "Multiple peptide locations")
})

test_that("extract_peptide_positions parses single and multi-location Positions.in.Master.Proteins", {
  qf <- QFeatures::QFeatures(list(peptides = SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      Positions.in.Master.Proteins = c(
        "P1 [5-13]",
        "P1 [5-13]; P2 [20-28]",
        NA_character_
      )
    )
  )))

  out <- extract_peptide_positions(qf, "peptides")
  rd <- SummarizedExperiment::rowData(out[["peptides"]])

  expect_equal(rd$start, c("5", "5;20", NA))
  expect_equal(rd$end, c("13", "13;28", NA))
})

test_that("extract_peptide_positions validates its arguments", {
  qf <- QFeatures::QFeatures(list(peptides = SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Positions.in.Master.Proteins = "P1 [1-2]")
  )))

  expect_error(extract_peptide_positions(qf), "must be specified")
  expect_error(extract_peptide_positions(qf, i = c(1, 2)), "single integer or character")
  expect_error(extract_peptide_positions(qf, i = "peptides", start_col = "x", end_col = "x"),
               "must be different")
  expect_error(extract_peptide_positions(qf, i = "nonexistent"), "null object")

  qf_missing_col <- QFeatures::QFeatures(list(peptides = SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(other_col = "x")
  )))
  expect_error(extract_peptide_positions(qf_missing_col, "peptides"), "Positions.in.Master.Proteins")
})

test_that("add_tryptic_termini_flags marks peptide termini against theoretical trypsin cleavage sites", {
  fasta_path <- write_proteome_fasta()
  row_data <- S4Vectors::DataFrame(
    Leading.razor.protein = c("P1", "P1"),
    start = c(5, 6),   # PEPTIDEAK: valid tryptic start; +1 is not
    end = c(13, 13)    # PEPTIDEAK: valid tryptic end
  )
  assay_mat <- matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)

  out <- add_tryptic_termini_flags(se, fasta_path, missed_cleavages = c(0, 1, 2))
  rd <- SummarizedExperiment::rowData(out)

  expect_true(rd$N_tryp[[1]])
  expect_true(rd$C_tryp[[1]])
  expect_false(rd$N_tryp[[2]])
  expect_true(rd$C_tryp[[2]])
})
