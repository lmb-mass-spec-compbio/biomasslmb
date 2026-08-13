suppressPackageStartupMessages(library(QFeatures))

make_impute_qf <- function() {
  samples <- c("s1", "s2", "s3", "s4")
  cd <- S4Vectors::DataFrame(condition = c("ctrl", "ctrl", "trt", "trt"), row.names = samples)

  mat_unimp <- matrix(
    c(1, NA, 3, 4,
      NA, NA, 7, 8,
      1, 2, 3, 4),
    nrow = 3, byrow = TRUE, dimnames = list(c("f1", "f2", "f3"), samples)
  )
  mat_imp <- matrix(
    c(1, 2.5, 3, 4,
      5, 6, 7, 8,
      1, 2, 3, 4),
    nrow = 3, byrow = TRUE, dimnames = list(c("f1", "f2", "f3"), samples)
  )
  rd <- S4Vectors::DataFrame(Protein = c("P1", "P2", "P3"))
  se_unimp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat_unimp), rowData = rd)
  se_imp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat_imp), rowData = rd)

  QFeatures(list(unimputed = se_unimp, imputed = se_imp), colData = cd)
}

test_that("restrict_imputation only uses imputed values where use_imputed_df says to", {
  qf <- make_impute_qf()
  use_imputed_df <- data.frame(condition = "ctrl", n_finite = 0)

  out <- suppressMessages(restrict_imputation(qf, "unimputed", "imputed", "restricted", use_imputed_df))
  restricted <- SummarizedExperiment::assay(out[["restricted"]])

  # f1/ctrl has 1 finite value (n_finite=1), not 0 -> not imputed, s2 stays NA
  expect_true(is.na(restricted["f1", "s2"]))
  # f2/ctrl has 0 finite values -> imputed values used for s1, s2
  expect_equal(unname(restricted["f2", c("s1", "s2")]), c(5, 6))
  # trt group not covered by use_imputed_df -> never imputed, unchanged from non-imputed data
  expect_equal(unname(restricted["f2", c("s3", "s4")]),
               unname(SummarizedExperiment::assay(qf[["unimputed"]])["f2", c("s3", "s4")]))
})

test_that("restrict_imputation errors if a grouping column is missing from colData", {
  qf <- make_impute_qf()
  use_imputed_df <- data.frame(nonexistent = "ctrl", n_finite = 0)

  expect_error(restrict_imputation(qf, "unimputed", "imputed", "restricted", use_imputed_df),
               "Not all grouping columns")
})

test_that("create_long_form_imputed_data flags which values were filled in by imputation", {
  qf <- make_impute_qf()

  out <- create_long_form_imputed_data(
    qf, "unimputed", "imputed",
    column_variables = "condition",
    row_variables = "Protein"
  )

  expect_true(is.data.frame(out))
  expect_true("is.imputed" %in% colnames(out))

  f2_s1 <- out[out$Protein == "P2" & out$colname == "s1", ]
  expect_true(f2_s1$is.imputed)
  expect_equal(f2_s1$value, 5)

  f3_s1 <- out[out$Protein == "P3" & out$colname == "s1", ]
  expect_false(f3_s1$is.imputed)
})
