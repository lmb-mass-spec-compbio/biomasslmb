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

# Ambiguous PTM site grouping -------------------------------------------------

# Peptides are placed so that peptide-relative positions 3, 4 and 8 land on
# protein positions 32, 33 and 37, which are the residues every case below is
# built from.
mq_phospho_se <- function(sequence, probabilities, n_ptms, start, protein = "P1") {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      1, nrow = length(sequence), ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      Sequence = sequence,
      Phospho..STY..Probabilities = probabilities,
      Phospho..STY. = n_ptms,
      Leading.razor.protein = protein,
      start = start)
  )
}

# a peptide whose phosphate is tied between positions 32 and 33
pep_32_33 <- list(sequence = "AASSAAAAAK",
                  probabilities = "AAS(0.5)S(0.5)AAAAAK",
                  n_ptms = 1, start = "30")

# a peptide starting one residue later, tied between 33 and 37
pep_33_37 <- list(sequence = "AASAAASAAK",
                  probabilities = "AAS(0.5)AAAS(0.5)AAK",
                  n_ptms = 1, start = "31")

test_that("parse_ptm_candidates_mq reads position, residue and probability, ignoring unmodified rows", {
  se <- mq_phospho_se(
    sequence = c("AASSAAAAAK", "AAAAAAAAAK"),
    probabilities = c("AAS(0.75)S(0.25)AAAAAK", ""),
    n_ptms = c(1, 0),
    start = c("30", "30"))

  candidates <- parse_ptm_candidates_mq(se)

  expect_equal(candidates$row, c(1, 1))
  expect_equal(candidates$pep_pos, c(3, 4))
  expect_equal(candidates$residue, c("S", "S"))
  expect_equal(candidates$prob, c(0.75, 0.25))
  expect_equal(candidates$n_ptms, c(1L, 1L))
})

test_that("parse_ptm_candidates_pd rescales percentages and recovers the PTM count", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      1, nrow = 2, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      ptmRS.Phospho.Site.Probabilities = c(
        "S(3): 50.0; S(4): 30.0; T(8): 20.0",
        "S(3): 100.0; S(4): 100.0"))
  )

  candidates <- parse_ptm_candidates_pd(se)

  expect_equal(candidates$row, c(1, 1, 1, 2, 2))
  expect_equal(candidates$pep_pos, c(3, 4, 8, 3, 4))
  expect_equal(candidates$residue, c("S", "S", "T", "S", "S"))
  expect_equal(candidates$prob, c(0.5, 0.3, 0.2, 1, 1))
  expect_equal(candidates$n_ptms, c(1L, 1L, 1L, 2L, 2L))
})

test_that("parse_ptm_candidates_pd yields no candidates for ptmRS non-results", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      1, nrow = 4, ncol = 1, dimnames = list(NULL, "s1"))),
    rowData = S4Vectors::DataFrame(
      ptmRS.Phospho.Site.Probabilities = c(
        "Inconclusive data", "Too many isoforms", "", "S(3): 100.0"))
  )

  candidates <- parse_ptm_candidates_pd(se)

  expect_equal(nrow(candidates), 1)
  expect_equal(candidates$row, 4)

  # and the same columns come back when nothing at all parses
  empty <- parse_ptm_candidates_pd(se[1:3, ])
  expect_equal(nrow(empty), 0)
  expect_equal(colnames(empty),
               c("row", "pep_pos", "residue", "prob", "n_ptms"))
})

test_that("add_ambiguous_ptm_group_rowdata chains overlapping candidate sets into one group", {
  se <- mq_phospho_se(
    sequence = c(pep_32_33$sequence, pep_33_37$sequence),
    probabilities = c(pep_32_33$probabilities, pep_33_37$probabilities),
    n_ptms = c(1, 1),
    start = c(pep_32_33$start, pep_33_37$start))

  rd <- SummarizedExperiment::rowData(
    add_ambiguous_ptm_group_rowdata(se, parse_ptm_candidates_mq(se), verbose = FALSE))

  # 32-or-33 and 33-or-37 share position 33, so they are one group of three
  expect_equal(rd$ptm_group_members, c("32;33;37", "32;33;37"))
  expect_equal(rd$ptm_group_n_candidates, c(3L, 3L))
  expect_equal(rd$ptm_group_resolved, c(FALSE, FALSE))
  expect_equal(length(unique(rd$ptm_group_id)), 1)
})

test_that("add_ambiguous_ptm_group_rowdata keeps peptides with different PTM counts apart", {
  se <- mq_phospho_se(
    sequence = c(pep_32_33$sequence, pep_33_37$sequence, "AASSAAASAAK"),
    probabilities = c(pep_32_33$probabilities, pep_33_37$probabilities,
                      "AAS(0.7)S(0.7)AAAS(0.6)AAK"),
    n_ptms = c(1, 1, 2),
    start = c(pep_32_33$start, pep_33_37$start, "30"))

  rd <- SummarizedExperiment::rowData(
    add_ambiguous_ptm_group_rowdata(se, parse_ptm_candidates_mq(se), verbose = FALSE))

  # the doubly modified peptide covers the same residues but is a different
  # species, so it gets its own group rather than joining the singly modified one
  expect_equal(rd$ptm_group_members, rep("32;33;37", 3))
  expect_equal(rd$ptm_group_id[1], rd$ptm_group_id[2])
  expect_false(rd$ptm_group_id[3] == rd$ptm_group_id[1])
})

test_that("add_ambiguous_ptm_group_rowdata leaves a localised peptide on its own residue", {
  se <- mq_phospho_se(
    sequence = c(pep_32_33$sequence, pep_33_37$sequence, "AASSAAAAAK"),
    probabilities = c(pep_32_33$probabilities, pep_33_37$probabilities,
                      "AAS(0.99)S(0.01)AAAAAK"),
    n_ptms = c(1, 1, 1),
    start = c(pep_32_33$start, pep_33_37$start, "30"))

  rd <- SummarizedExperiment::rowData(
    add_ambiguous_ptm_group_rowdata(se, parse_ptm_candidates_mq(se), verbose = FALSE))

  # position 32 is confidently localised on row 3 and a candidate in the pooled
  # group; the resolved peptide keeps its own single-residue identity
  expect_equal(rd$ptm_group_members[3], "32")
  expect_equal(rd$ptm_group_n_candidates[3], 1L)
  expect_true(rd$ptm_group_resolved[3])
  expect_equal(rd$ptm_group_members[1:2], c("32;33;37", "32;33;37"))
})

test_that("add_ambiguous_ptm_group_rowdata drops candidates below min_candidate_prob", {
  se <- mq_phospho_se(
    sequence = "AASSAAASAAK",
    probabilities = "AAS(0.495)S(0.495)AAAS(0.01)AAK",
    n_ptms = 1,
    start = "30")
  candidates <- parse_ptm_candidates_mq(se)

  at_default <- SummarizedExperiment::rowData(
    add_ambiguous_ptm_group_rowdata(se, candidates, verbose = FALSE))
  expect_equal(at_default$ptm_group_members, "32;33")

  below <- SummarizedExperiment::rowData(add_ambiguous_ptm_group_rowdata(
    se, candidates, min_candidate_prob = 0.005, verbose = FALSE))
  expect_equal(below$ptm_group_members, "32;33;37")
})

test_that("add_ambiguous_ptm_group_rowdata leaves peptides with no single position unassigned", {
  se <- mq_phospho_se(
    sequence = c(pep_32_33$sequence, "AAAAAAAAAK"),
    # a peptide occurring twice in its protein, and a peptide with no PTM
    probabilities = c(pep_32_33$probabilities, ""),
    n_ptms = c(1, 0),
    start = c("30;120", "30"))

  rd <- SummarizedExperiment::rowData(
    add_ambiguous_ptm_group_rowdata(se, parse_ptm_candidates_mq(se), verbose = FALSE))

  expect_true(all(is.na(rd$ptm_group_id)))
  expect_true(all(is.na(rd$ptm_group_resolved)))
})

test_that("add_ambiguous_ptm_group_rowdata errors when peptide positions are missing", {
  se <- mq_phospho_se(pep_32_33$sequence, pep_32_33$probabilities,
                      pep_32_33$n_ptms, pep_32_33$start)
  SummarizedExperiment::rowData(se)$start <- NULL

  expect_error(
    add_ambiguous_ptm_group_rowdata(se, parse_ptm_candidates_mq(se), verbose = FALSE),
    "add_peptide_positions_from_cleavage")
})

test_that("split_wide_components cuts a component at its largest internal gap", {
  nodes <- c(10, 20, 30, 100, 110)

  membership <- split_wide_components(nodes, rep(1L, 5), max_group_span = 50)

  expect_equal(unname(split(nodes, membership)),
               list(c(10, 20, 30), c(100, 110)))
})

test_that("split_wide_components cuts repeatedly when one cut is not enough", {
  nodes <- c(0, 10, 20, 100, 110, 200)

  membership <- split_wide_components(nodes, rep(1L, 6), max_group_span = 50)

  # cutting at the widest gap leaves 0-110 still too wide, so it is cut again
  expect_equal(sort(unname(sapply(split(nodes, membership), min))), c(0, 100, 200))
  expect_true(all(tapply(nodes, membership, function(x) max(x) - min(x)) <= 50))
})

test_that("split_wide_components leaves a component within the span untouched", {
  nodes <- c(10, 20, 30)

  expect_equal(split_wide_components(nodes, rep(1L, 3), max_group_span = 50),
               rep(1L, 3))
})

test_that("add_ambiguous_ptm_group_rowdata unassigns a peptide split across two groups", {
  se <- mq_phospho_se(
    sequence = c(pep_32_33$sequence, pep_33_37$sequence),
    probabilities = c(pep_32_33$probabilities, pep_33_37$probabilities),
    n_ptms = c(1, 1),
    start = c(pep_32_33$start, pep_33_37$start))

  rd <- SummarizedExperiment::rowData(add_ambiguous_ptm_group_rowdata(
    se, parse_ptm_candidates_mq(se), max_group_span = 3, verbose = FALSE))

  # the 32;33;37 component spans 5 residues, so it is cut between 33 and 37.
  # The 33-or-37 peptide then straddles both pieces and has no single group
  expect_equal(rd$ptm_group_members[1], "32;33")
  expect_true(is.na(rd$ptm_group_id[2]))
})

test_that("summarise_ptm_groups reports pooled group size and the mass left behind", {
  se <- mq_phospho_se(
    sequence = "AASSAAASAAK",
    probabilities = "AAS(0.49)S(0.49)AAAS(0.02)AAK",
    n_ptms = 1,
    start = "30")
  candidates <- parse_ptm_candidates_mq(se)

  kept <- summarise_ptm_groups(se, candidates, min_candidate_prob = 0.005)
  expect_equal(kept$pooled_groups, 1)
  expect_equal(kept$max_size, 3L)
  expect_equal(kept$pc_miss_over_1, 0)

  # dropping the 0.02 candidate leaves 2% of the probability mass outside the
  # group, so the group now has a 2% chance of excluding the true site
  dropped <- summarise_ptm_groups(se, candidates, min_candidate_prob = 0.03)
  expect_equal(dropped$max_size, 2L)
  expect_equal(dropped$pc_miss_over_1, 100)
  expect_equal(dropped$pc_miss_over_5, 0)
})
