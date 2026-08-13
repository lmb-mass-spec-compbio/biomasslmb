psm <- tmt_qf[["psms_filtered"]]
prot <- tmt_qf[["protein"]]

test_that("get_n_feature_per_prot tallies finite features per protein per sample", {
  out <- biomasslmb:::get_n_feature_per_prot(psm)

  expect_s3_class(out, "data.frame")
  expect_setequal(colnames(out), c("Master.Protein.Accessions", "sample", "n"))
  expect_true(all(out$n >= 1))

  # cross check one protein/sample combination by hand
  first_protein <- out$Master.Protein.Accessions[[1]]
  first_sample <- out$sample[[1]]
  rd <- SummarizedExperiment::rowData(psm)
  feature_rows <- rd$Master.Protein.Accessions == first_protein
  expected_n <- sum(is.finite(SummarizedExperiment::assay(psm)[feature_rows, first_sample]))
  expect_equal(out$n[out$Master.Protein.Accessions == first_protein & out$sample == first_sample], expected_n)
})

test_that("filter_features_per_protein keeps proteins with >= min_features in at least one sample", {
  n_per_prot <- biomasslmb:::get_n_feature_per_prot(psm)

  min_features <- 10
  expected_proteins <- unique(n_per_prot$Master.Protein.Accessions[n_per_prot$n >= min_features])

  out <- filter_features_per_protein(psm, min_features = min_features)

  expect_setequal(unique(SummarizedExperiment::rowData(out)$Master.Protein.Accessions), expected_proteins)
})

test_that("get_protein_no_quant_mask flags proteins with too few supporting features per sample as NA", {
  min_features <- 3
  mask <- get_protein_no_quant_mask(psm, min_features = min_features)

  expect_true(is.matrix(mask) || is.logical(mask))
  # values are either TRUE or NA, never FALSE
  expect_true(all(mask[!is.na(mask)] == TRUE))
})

test_that("mask_protein_level_quant replaces protein-level values with NA where the mask says so", {
  mask <- get_protein_no_quant_mask(psm, min_features = 3)

  out <- mask_protein_level_quant(prot, mask)

  common_proteins <- intersect(rownames(prot), rownames(mask))
  common_samples <- intersect(colnames(prot), colnames(mask))
  expected <- SummarizedExperiment::assay(prot)[common_proteins, common_samples] *
    mask[common_proteins, common_samples]

  expect_equal(SummarizedExperiment::assay(out), as.matrix(expected))
})

test_that("mask_protein_level_quant errors if obj has proteins missing from retain_mask", {
  mask <- get_protein_no_quant_mask(psm, min_features = 3)
  bad_prot <- prot[1, ]
  rownames(bad_prot) <- "NOT_A_PROTEIN"

  expect_error(mask_protein_level_quant(bad_prot, mask), "not in")
})
