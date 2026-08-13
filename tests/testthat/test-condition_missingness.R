suppressPackageStartupMessages(library(QFeatures))

make_test_qf <- function() {
  samples <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "IS")
  cd <- S4Vectors::DataFrame(condition = c(rep("A", 4), rep("B", 4), NA), row.names = samples)

  mat <- matrix(
    c(
      NA, NA, NA, NA, 10, 11, 12, 13, 5,   # f1: perfectly condition-structured
      1,  2,  3,  4,  5,  6,  7,  8,  9,   # f2: fully observed -> uninformative
      NA, NA, NA, NA, NA, NA, NA, 8,  NA,  # f3: only 1 observed -> uninformative
      1,  NA, 3,  NA, 5,  NA, 7,  NA, 9    # f4: missingness unrelated to condition
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("f1", "f2", "f3", "f4"), samples)
  )
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  QFeatures(list(psms = se), colData = cd)
}

test_that("filter_complete_groups drops samples with NA in the group column", {
  qf <- make_test_qf()
  out <- filter_complete_groups(qf, "psms", "condition")

  expect_equal(out$removed, "IS")
  expect_false("IS" %in% colnames(out$mat))
  expect_false("IS" %in% rownames(out$cd_assay))
  expect_equal(ncol(out$mat), 8)
})

test_that("filter_complete_groups errors on an unknown group column", {
  qf <- make_test_qf()
  expect_error(filter_complete_groups(qf, "psms", "nonexistent"), "not found in colData")
})

test_that("condition_miss_score gives a perfect score to fully condition-structured missingness", {
  qf <- make_test_qf()
  res <- suppressMessages(condition_miss_score(qf, "psms", "condition"))

  expect_equal(unname(res$scores["f1"]), 1)
  expect_equal(res$summary["f1", "condition_miss_class"], "condition_structured")
})

test_that("condition_miss_score marks fully-observed and near-empty features as uninformative", {
  qf <- make_test_qf()
  res <- suppressMessages(condition_miss_score(qf, "psms", "condition"))

  expect_true(is.na(res$scores["f2"])) # no missingness at all
  expect_true(is.na(res$scores["f3"])) # < min_observed
  expect_equal(res$summary["f2", "condition_miss_class"], "uninformative")
  expect_equal(res$summary["f3", "condition_miss_class"], "uninformative")
})

test_that("condition_miss_score gives a low score to condition-independent missingness", {
  qf <- make_test_qf()
  res <- suppressMessages(condition_miss_score(qf, "psms", "condition"))

  expect_lt(res$scores["f4"], 0.5)
})

test_that("condition_miss_score writes results into rowData when store_results = TRUE", {
  qf <- make_test_qf()
  res <- suppressMessages(condition_miss_score(qf, "psms", "condition", store_results = TRUE))

  rd <- SummarizedExperiment::rowData(res$obj[["psms"]])
  expect_true("condition_miss_score" %in% colnames(rd))
  expect_true("condition_miss_class" %in% colnames(rd))
})

test_that("condition_miss_score validates min_observed/min_missing and group_cols", {
  qf <- make_test_qf()
  expect_error(condition_miss_score(qf, "psms", "nonexistent"), "not found in colData")
  expect_error(condition_miss_score(qf, "psms", "condition", min_observed = 0), "min_observed")
  expect_error(condition_miss_score(qf, "psms", "condition", min_missing = 0), "min_missing")
})

test_that("global_condition_miss_score decomposes Tjur R2 into intensity and incremental parts", {
  qf <- make_test_qf()
  res <- suppressMessages(global_condition_miss_score(qf, "psms", "condition", log_transform = FALSE))

  expect_equal(res$tjur_full, res$tjur_intensity_only + res$tjur_incremental, tolerance = 1e-8)
  expect_true(res$tjur_incremental > 0) # f1's condition-structured missingness should show up
  expect_type(res$lrt_pvalue, "double")
})

test_that("condition_miss_index aggregates per-feature scores with miss_frac weighting and coverage penalty", {
  qf <- make_test_qf()
  res <- suppressMessages(condition_miss_score(qf, "psms", "condition"))
  idx <- suppressMessages(condition_miss_index(res$summary))

  expect_equal(idx$n_informative, 2) # f1 and f4
  expect_equal(idx$n_total, 4)
  expect_equal(idx$coverage, 0.5)
  expect_equal(idx$index, idx$weighted_mean_score * idx$coverage, tolerance = 1e-8)
})

test_that("condition_miss_index without coverage_penalty ignores the informative-feature fraction", {
  qf <- make_test_qf()
  res <- suppressMessages(condition_miss_score(qf, "psms", "condition"))
  idx <- suppressMessages(condition_miss_index(res$summary, coverage_penalty = FALSE))

  expect_equal(idx$index, idx$weighted_mean_score, tolerance = 1e-8)
})

test_that("condition_miss_index errors if summary_df is missing required columns", {
  # The message must name the columns that are actually missing, not just say
  # that some are
  expect_error(
    condition_miss_index(data.frame(condition_miss_score = 1, miss_frac = 0)),
    "missing required column\\(s\\): condition_miss_class"
  )
  expect_error(condition_miss_index(data.frame(x = 1)), "missing required column")
})

test_that("condition_miss_index returns NA when no features are informative", {
  summary_df <- data.frame(
    condition_miss_score = NA_real_,
    miss_frac = 0,
    condition_miss_class = "uninformative"
  )
  expect_warning(idx <- condition_miss_index(summary_df), "No informative features")
  expect_true(is.na(idx$index))
})

test_that("mnar_score/mnar_global_score/mnar_index are deprecated wrappers around the new names", {
  qf <- make_test_qf()

  expect_warning(res_old <- suppressMessages(mnar_score(qf, "psms", "condition")), "deprecated")
  res_new <- suppressMessages(condition_miss_score(qf, "psms", "condition"))
  expect_equal(res_old$scores, res_new$scores)

  expect_warning(g_old <- suppressMessages(mnar_global_score(qf, "psms", "condition")), "deprecated")
  g_new <- suppressMessages(global_condition_miss_score(qf, "psms", "condition"))
  expect_equal(g_old$tjur_full, g_new$tjur_full)

  expect_warning(i_old <- suppressMessages(mnar_index(res_new$summary)), "deprecated")
  i_new <- suppressMessages(condition_miss_index(res_new$summary))
  expect_equal(i_old$index, i_new$index)
})
