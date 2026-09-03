# Smoke tests for plotting functions: verify they run without error and
# return the expected plot object class, using the real bundled tmt_qf data
# and small synthetic data.frames for the functions that take one directly.

suppressPackageStartupMessages(library(QFeatures))
psm <- tmt_qf[["psms_raw"]]
qf <- QFeatures(list(psms = psm))
colData(qf)$condition <- rep(c("A", "B"), length.out = ncol(qf))

test_that("theme_biomasslmb returns a ggplot theme and toggles aspect ratio / border", {
  th <- theme_biomasslmb()
  expect_s3_class(th, "theme")
  expect_true(!is.null(th$aspect.ratio))

  th_no_square <- theme_biomasslmb(aspect_square = FALSE)
  expect_null(th_no_square$aspect.ratio)

  th_border <- theme_biomasslmb(border = TRUE)
  expect_true(!is.null(th_border$panel.border))
})

test_that("plot_pca returns a ggplot, coloured and shaped by a colData column", {
  p <- plot_pca(qf, "psms", colour_by = "condition", shape_by = "condition")
  expect_s3_class(p, "ggplot")
})

test_that("plot_pca errors for an unknown colour_by/shape_by column", {
  expect_error(plot_pca(qf, "psms", colour_by = "nonexistent"), "not in the colData")
})

test_that("plot_rt_vs_delta and plot_rt_dist return ggplot objects", {
  p1 <- plot_rt_vs_delta(psm, rt_col = "RT.in.min", delta_ppm_col = "Delta.M.in.ppm")
  expect_s3_class(p1, "ggplot")

  p2 <- plot_rt_dist(psm, rt_col = "RT.in.min")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_volcano returns a ggplot and adds a fill aesthetic / title when requested", {
  volc_df <- data.frame(
    logFC = rnorm(50), adj.P.Val = runif(50), sig = sample(c(TRUE, FALSE), 50, TRUE)
  )
  p <- plot_volcano(volc_df, sig_col = "sig", title = "test")
  expect_s3_class(p, "ggplot")
  expect_true("fill" %in% names(p$mapping))
  expect_equal(p$labels$title, "test")
})

test_that("plot_protein_assays returns a ggplot for a set of proteins of interest", {
  poi <- unique(SummarizedExperiment::rowData(psm)$Master.Protein.Accessions)[1:2]
  p <- plot_protein_assays(qf, poi = poi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_quant returns a ggplot for each supported method", {
  expect_s3_class(plot_quant(psm, method = "box"), "ggplot")
  expect_s3_class(plot_quant(psm, method = "density"), "ggplot")
  expect_s3_class(plot_quant(psm, method = "histogram", log2transform = TRUE), "ggplot")
})

test_that("plot_missing_SN and plot_missing_SN_per_sample return ggplot objects", {
  p1 <- plot_missing_SN(psm)
  expect_s3_class(p1, "ggplot")

  p2 <- plot_missing_SN_per_sample(psm)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_missing_upset returns an upset plot", {
  p <- plot_missing_upset(qf, "psms")
  expect_s3_class(p, "upset")
})

test_that("plot_missing_upset passes ... on, including over its own defaults", {
  # nintersects, sets and keep.order are supplied inside the function, so
  # passing them at all depends on ... taking precedence over those defaults
  expect_s3_class(plot_missing_upset(qf, "psms", nintersects = 5), "upset")

  na_cols <- colnames(data.frame(assay(psm)))[colSums(is.na(assay(psm))) > 0]
  wanted <- paste0(sort(na_cols)[1:2], "_NA")
  p <- plot_missing_upset(qf, "psms", sets = wanted)
  expect_equal(p$Set_names, wanted)
})

test_that("plot_cor_samples runs without error and returns the sample correlation matrix", {
  # returns the corrplot() result (a list including a $corr matrix), not a
  # ggplot, despite the docs -- verify it runs cleanly and gives a sensible result.
  result <- plot_cor_samples(qf, "psms")
  expect_true(is.matrix(result$corr))
  expect_equal(dim(result$corr), c(ncol(psm), ncol(psm)))
})

test_that("plot_go_terms_upset returns an upset plot for GO terms of interest", {
  gene2cat <- data.frame(
    UNIPROTKB = c("P1", "P1", "P2"),
    TERM = c("apoptotic process", "cell cycle", "apoptotic process"),
    stringsAsFactors = FALSE
  )
  p <- plot_go_terms_upset(c("apoptotic process", "cell cycle"), c("P1", "P2"), gene2cat)
  expect_s3_class(p, "upset")
})
