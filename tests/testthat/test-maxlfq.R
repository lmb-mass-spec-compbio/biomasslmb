test_that("maxlfq_wrapper recovers the relative protein abundance pattern from peptide ratios", {
  skip_if_not_installed("iq")

  set.seed(1)
  protein_true <- c(s1 = 10, s2 = 10.5, s3 = 9.8, s4 = 11)
  # 3 peptides, each offset from the true protein level but consistent in
  # their sample-to-sample ratios -- MaxLFQ should recover the relative
  # differences between samples (not the absolute peptide-level offsets).
  mat <- rbind(
    protein_true,
    protein_true + 2,
    protein_true - 1
  )
  colnames(mat) <- names(protein_true)
  rownames(mat) <- paste0("pep", 1:3)

  out <- maxlfq_wrapper(mat)

  expect_length(out, 4)
  expect_equal(diff(out), unname(diff(protein_true)), tolerance = 1e-6)
})

test_that("maxlfq_wrapper errors on a non-numeric matrix", {
  skip_if_not_installed("iq")

  mat <- matrix(c("a", "b", "c", "d"), nrow = 2, ncol = 2)
  expect_error(maxlfq_wrapper(mat), "must be a numeric matrix")
})
