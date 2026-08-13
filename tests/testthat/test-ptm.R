write_protein_fasta <- function() {
  path <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c("sp|P1|PROT1_HUMAN" = "MSAKPEPTIDEAKGGVPEPTIDEBRAA")),
    filepath = path
  )
  path
}

test_that("filter_maxquant_ptm keeps only rows with a non-empty PTM probability column", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Phospho..STY..Probabilities = c("", "AAS(0.5)AA", ""))
  )
  qf <- QFeatures::QFeatures(list(pep = se))

  out <- filter_maxquant_ptm(qf, "pep")
  expect_equal(nrow(out), 1)
})

test_that("add_ptm_pos_rowdata_mq locates a 'p'-prefixed phospho tag before its residue", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Sequence = "AASLE", Modified.sequence = "_AApSLE_")
  )

  out <- add_ptm_pos_rowdata_mq(se, verbose = FALSE)
  rd <- SummarizedExperiment::rowData(out)

  expect_equal(rd$ptms, "p")
  expect_equal(rd$ptm_positions, "3")
  expect_equal(rd$ptm_amino_acids, "S")
})

test_that("add_ptm_pos_rowdata_mq locates an '(ox)'-suffixed tag after its residue", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Sequence = "AMPLES", Modified.sequence = "_AM(ox)PLES_")
  )

  out <- add_ptm_pos_rowdata_mq(se, ptms_to_retain = "(ox)", verbose = FALSE)
  rd <- SummarizedExperiment::rowData(out)

  expect_equal(rd$ptms, "(ox)")
  expect_equal(rd$ptm_positions, "2")
  expect_equal(rd$ptm_amino_acids, "M")
})

test_that("add_ptm_pos_rowdata_mq errors on a PTM with no defined encoding position", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Sequence = "AASLE", Modified.sequence = "_AApSLE_")
  )

  expect_error(
    add_ptm_pos_rowdata_mq(se, ptms_to_retain = "unknown_ptm", verbose = FALSE),
    "Need to define position"
  )
})

test_that("add_filter_ptm_pos_rowdata_mq extracts site, position and probability from PD-style columns", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      Sequence = c("AASAAK", "AASAAK"),
      Modified.sequence = c("_AAS(Phospho (STY))AAK_", "_AAS(Phospho (STY))AAK_"),
      Phospho..STY..Probabilities = c("AAS(0.995)AAK", "AAS(0.2)AAK")
    )
  )

  out <- add_filter_ptm_pos_rowdata_mq(se, verbose = FALSE, filter_pep_by_prob = FALSE)
  rd <- SummarizedExperiment::rowData(out)

  expect_equal(rd$ptms, c("Phospho (STY)", ""))
  expect_equal(rd$ptm_positions, c("3", ""))
  expect_equal(rd$ptm_amino_acids, c("S", ""))
  expect_equal(rd$n_ptms, c("1", "0"))
  expect_equal(rd$n_ptms_detected, c("1", "1"))
})

test_that("add_filter_ptm_pos_rowdata_mq with filter_pep_by_prob=TRUE drops peptides below min_prob", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      Sequence = c("AASAAK", "AASAAK"),
      Modified.sequence = c("_AAS(Phospho (STY))AAK_", "_AAS(Phospho (STY))AAK_"),
      Phospho..STY..Probabilities = c("AAS(0.995)AAK", "AAS(0.2)AAK")
    )
  )

  out <- add_filter_ptm_pos_rowdata_mq(se, verbose = FALSE, filter_pep_by_prob = TRUE)
  expect_equal(nrow(out), 1)
  expect_equal(SummarizedExperiment::rowData(out)$ptm_amino_acids, "S")
})

test_that("add_peptide_positions_from_cleavage locates a tryptic peptide's start/end in the protein", {
  fasta_path <- write_protein_fasta()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Leading.razor.protein = "P1", Sequence = "GGVPEPTIDEBR")
  )

  out <- add_peptide_positions_from_cleavage(se, fasta_path)
  rd <- SummarizedExperiment::rowData(out)

  expect_equal(rd$start, "14")
  expect_equal(rd$end, "25")
})

test_that("add_ptm_positions converts a peptide-relative PTM position into a protein-absolute site name", {
  fasta_path <- write_protein_fasta()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      Leading.razor.protein = "P1",
      Sequence = "GGVPEPTIDEBR",
      ptm_positions = "5", # 5th residue of GGVPEPTIDEBR = E
      ptm_amino_acids = "E"
    )
  )

  out <- add_ptm_positions(se, fasta_path)
  rd <- SummarizedExperiment::rowData(out)

  expect_equal(rd$ptm_positions_prot, "18")
  expect_equal(rd$ptm_name, "E18")
})

test_that("parse_PTM_scores_pd extracts amino acid, position and score, filtering on threshold", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      ptmRS.Best.Site.Probabilities = c("Inconclusive data", "S6(Phospho): 99.5", "S6(Phospho): 20")
    )
  )

  out <- parse_PTM_scores_pd(se, threshold = 95, verbose = FALSE)
  rd <- SummarizedExperiment::rowData(out)

  # "Inconclusive data" row is dropped entirely
  expect_equal(nrow(out), 2)
  expect_equal(rd$ptm_amino_acids, c("S", ""))
  expect_equal(rd$ptm_positions, c("6", ""))
  expect_equal(rd$ptms, c("Phospho", ""))
  expect_equal(rd$ptm_scores, c("99.5", ""))
})

test_that("get_sequence extracts a padded, lower-cased site sequence around a PTM", {
  proteome <- Biostrings::AAStringSet(c(P1 = paste(rep("A", 20), collapse = "")))
  proteome[["P1"]][3] <- Biostrings::AAString("M")

  expect_equal(get_sequence(proteome, "P1", 3, pad = 3), "_AAmAAA")
  expect_equal(get_sequence(proteome, "P1", 1, pad = 3), "___aAMA")
  expect_equal(get_sequence(proteome, "P1", 20, pad = 3), "AAAa___")
})

test_that("get_sequence returns NA for a missing protein, NA position, or out-of-range position", {
  proteome <- Biostrings::AAStringSet(c(P1 = paste(rep("A", 20), collapse = "")))

  expect_true(is.na(get_sequence(proteome, "P2", 3, pad = 3)))
  expect_true(is.na(get_sequence(proteome, "P1", NA, pad = 3)))
  expect_warning(result <- get_sequence(proteome, "P1", 100, pad = 3), "outside protein sequence")
  expect_true(is.na(result))
})

test_that("add_site_sequence adds a rowData column with the sequence around ptm_positions_prot", {
  proteome_path <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c("sp|P1|X" = paste(rep("A", 20), collapse = ""))),
    filepath = proteome_path
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(Master.Protein.Accessions = "P1", ptm_positions_prot = "10")
  )

  out <- add_site_sequence(se, proteome_path, sequence_pad = 2)
  expect_equal(SummarizedExperiment::rowData(out)$site_seq, "AAaAA")
})
