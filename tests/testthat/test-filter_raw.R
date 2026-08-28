psm <- tmt_qf[["psms_raw"]]
peptides <- dia_qf[["precursors"]]

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

test_that("filter_features_pd_dda errors without contaminant_proteins when filter_contaminant=TRUE", {
  expect_error(
    suppressMessages(filter_features_pd_dda(psm, filter_contaminant = TRUE)),
    "must supply the contaminant_proteins"
  )
})

test_that("filter_features_pd_dda errors if both crap_proteins and contaminant_proteins are given", {
  expect_error(
    filter_features_pd_dda(psm, crap_proteins = "P1", contaminant_proteins = "P2"),
    "Please use only contaminant_proteins"
  )
})

test_that("remove_contaminant drops features matching a contaminant accession", {
  contaminant_id <- SummarizedExperiment::rowData(psm)[["Master.Protein.Accessions"]][[1]]
  n_matching <- sum(SummarizedExperiment::rowData(psm)[["Master.Protein.Accessions"]] == contaminant_id)
  expect_gt(n_matching, 0)

  out <- suppressMessages(remove_contaminant(
    psm,
    contaminant_proteins = contaminant_id,
    filter_associated_contaminant = FALSE,
    master_protein_col = "Master.Protein.Accessions"
  ))

  expect_equal(nrow(out), nrow(psm) - n_matching)
  expect_false(any(SummarizedExperiment::rowData(out)[["Master.Protein.Accessions"]] == contaminant_id))
})

test_that("remove_contaminant errors without contaminant_proteins", {
  expect_error(
    remove_contaminant(psm, NULL, FALSE, "Master.Protein.Accessions"),
    "must supply the contaminant_proteins"
  )
})

test_that("filter_features_diann filters on Protein.Group and adds Number.of.Protein.Groups", {
  n_multi <- sum(grepl(";", SummarizedExperiment::rowData(peptides)[["Protein.Group"]]))
  expect_gt(n_multi, 0)

  out <- suppressMessages(filter_features_diann(
    peptides,
    filter_contaminant = FALSE,
    unique_master = TRUE,
    remove_no_quant = TRUE
  ))

  expect_s4_class(out, "SummarizedExperiment")
  expect_true(all(SummarizedExperiment::rowData(out)[["Number.of.Protein.Groups"]] == 1))
  expect_false(any(grepl(";", SummarizedExperiment::rowData(out)[["Protein.Group"]])))
  expect_true(all(rowSums(is.finite(SummarizedExperiment::assay(out))) > 0))
})

test_that("filter_features_diann errors without contaminant_proteins when filter_contaminant=TRUE", {
  expect_error(
    suppressMessages(filter_features_diann(peptides, filter_contaminant = TRUE)),
    "must supply the contaminant_proteins"
  )
})

## Synthetic Spectronaut-format data: filter_features_sn is not exercised by
## any vignette or bundled dataset, so we build a minimal SE by hand using the
## documented column names.
make_sn_se <- function() {
  row_data <- S4Vectors::DataFrame(
    PG.ProteinAccessions = c("P1", "P1;P2", "", "P3", "P4")
  )
  assay_mat <- matrix(
    c(1, 2, 3, NA, 5,
      1, 2, 3, NA, 5),
    nrow = 5, ncol = 2,
    dimnames = list(NULL, c("s1", "s2"))
  )
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)
}

test_that("filter_features_sn removes non-unique, missing-master and no-quant features", {
  se <- make_sn_se()
  out <- suppressMessages(filter_features_sn(se, filter_contaminant = FALSE))

  # "P1;P2" (non-unique), "" (no master) and "P3" (no finite quant) are removed;
  # "P1" and "P4" are unique, have a master protein, and have finite quant.
  expect_setequal(SummarizedExperiment::rowData(out)$PG.ProteinAccessions, c("P1", "P4"))
})

test_that("filter_features_sn errors without contaminant_proteins when filter_contaminant=TRUE", {
  se <- make_sn_se()
  expect_error(
    suppressMessages(filter_features_sn(se, filter_contaminant = TRUE)),
    "must supply the contaminant_proteins"
  )
})

## Synthetic MaxQuant-format data for filter_features_mq_dda / remove_contaminant_mq,
## which are not exercised by any vignette or bundled dataset.
make_mq_se <- function() {
  row_data <- S4Vectors::DataFrame(
    Leading.razor.protein = c("P1", "P2", "P3", "CON__P4", "P5"),
    Proteins = c("P1", "P2", "P3", "CON__P4", "P5"),
    Reverse = c("", "", "+", "", ""),
    Potential.contaminant = c("", "", "", "", ""),
    Unique..Proteins. = c("yes", "yes", "yes", "yes", "no")
  )
  assay_mat <- matrix(
    c(1, 2, 3, 4, 5,
      1, 2, 3, 4, NA),
    nrow = 5, ncol = 2,
    dimnames = list(NULL, c("s1", "s2"))
  )
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)
}

test_that("filter_features_mq_dda removes decoys, contaminants and non-unique/no-quant features", {
  se <- make_mq_se()
  out <- suppressMessages(filter_features_mq_dda(se, contaminant_proteins = character(0),
                                                  filter_associated_contaminant = FALSE))

  # Reverse == "+" (P3) and CON__ contaminant (P4) and Unique..Proteins.=="no" (P5) all removed
  expect_setequal(SummarizedExperiment::rowData(out)$Leading.razor.protein, c("P1", "P2"))
})

test_that("remove_contaminant_mq identifies CON__ prefixed proteins as contaminants", {
  se <- make_mq_se()
  out <- suppressMessages(remove_contaminant_mq(se, contaminant_proteins = character(0),
                                                 filter_associated_contaminant = FALSE))

  expect_false("CON__P4" %in% SummarizedExperiment::rowData(out)$Proteins)
  expect_equal(nrow(out), 4)
})

test_that("remove_contaminant_mq errors without contaminant_proteins", {
  se <- make_mq_se()
  expect_error(
    remove_contaminant_mq(se, NULL, FALSE),
    "must supply the contaminant_proteins"
  )
})
