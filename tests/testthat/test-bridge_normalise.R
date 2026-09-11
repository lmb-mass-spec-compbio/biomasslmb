# Build a SummarizedExperiment with a known per-feature, per-plex offset.
#
# `n_bridge` gives the number of bridge channels in each plex and `offsets` the
# per-feature offset applied to every channel of each plex, so that a correct
# bridge normalisation must recover `base_quant` up to a per-feature constant.
make_plex_se <- function(base_quant,
                         offsets,
                         n_bridge = NULL,
                         n_sample = NULL){

  plexes <- names(offsets)

  if(is.null(n_bridge)) n_bridge <- setNames(as.list(rep(1, length(plexes))), plexes)
  if(is.null(n_sample)) n_sample <- setNames(as.list(rep(3, length(plexes))), plexes)

  quant <- do.call(cbind, lapply(plexes, function(p){
    n <- n_bridge[[p]] + n_sample[[p]]
    matrix(base_quant + offsets[[p]], nrow = length(base_quant), ncol = n)
  }))

  col_data <- do.call(rbind, lapply(plexes, function(p){
    data.frame(
      Plex = p,
      Type = c(rep("Bridge", n_bridge[[p]]), rep("Sample", n_sample[[p]])))
  }))

  colnames(quant) <- paste(col_data$Plex, col_data$Type,
                           sequence(rle(paste0(col_data$Plex, col_data$Type))$lengths),
                           sep = "_")
  rownames(quant) <- names(base_quant)
  rownames(col_data) <- colnames(quant)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(quant = quant),
    colData = S4Vectors::DataFrame(col_data))
}

base_quant <- c(f1 = 20, f2 = 18, f3 = 22, f4 = 19)

test_that("bridge_normalise recovers known per-plex offsets", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = c(0.5, -1, 0.25, 2),
                                    "2" = c(-0.5, 1, -0.25, -2)))

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          verbose = FALSE)

  # every channel of every plex should now hold the underlying value
  expect_equal(SummarizedExperiment::assay(out),
               matrix(base_quant, nrow = 4, ncol = 8,
                      dimnames = dimnames(SummarizedExperiment::assay(se))))
})

test_that("bridge_normalise centres on the mean of the plexes by default", {
  # An offset applied to every plex is not a plex effect and must be retained
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(3, 4), "2" = rep(3, 4)))

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          verbose = FALSE)

  expect_equal(unname(SummarizedExperiment::assay(out)[, 1]),
               unname(base_quant + 3))
})

test_that("bridge_normalise aligns to a nominated reference plex", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          reference = "2", verbose = FALSE)

  # everything is put onto plex 2's scale, rather than the mean of the two
  expect_equal(unname(SummarizedExperiment::assay(out)[, 1]),
               unname(base_quant - 1))
})

test_that("bridge_normalise weights plexes, not bridge channels, when counts are unequal", {
  # Plex 1 carries two bridge channels and plex 2 one. Pooling the three
  # channels would put the common scale at (1 + 1 + -1)/3 rather than
  # (1 + -1)/2, leaving a residual offset.
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)),
                     n_bridge = list("1" = 2, "2" = 1),
                     n_sample = list("1" = 3, "2" = 3))

  expect_warning(
    bridge_normalise(se, plex_col = "Plex",
                     bridge_cols = se$Type == "Bridge",
                     verbose = FALSE),
    "unequal numbers of bridge channels")

  out <- suppressWarnings(
    bridge_normalise(se, plex_col = "Plex",
                     bridge_cols = se$Type == "Bridge",
                     verbose = FALSE))

  expect_equal(unname(SummarizedExperiment::assay(out)[, 1]), unname(base_quant))
})

test_that("bridge_normalise collapses several bridge channels with fun", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)),
                     n_bridge = list("1" = 2, "2" = 2),
                     n_sample = list("1" = 2, "2" = 2))

  # disagreeing bridge channels within plex 1: mean is 1, so the correction is
  # unchanged, but a min/max would not be
  quant <- SummarizedExperiment::assay(se)
  quant[, 1] <- quant[, 1] - 0.5
  quant[, 2] <- quant[, 2] + 0.5
  SummarizedExperiment::assay(se) <- quant

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          verbose = FALSE)

  expect_equal(unname(SummarizedExperiment::assay(out)[, 3]), unname(base_quant))
})

test_that("bridge_normalise on a single plex is a no-op", {
  se <- make_plex_se(base_quant, offsets = list("1" = c(0.5, -1, 0.25, 2)),
                     n_bridge = list("1" = 1), n_sample = list("1" = 3))

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          verbose = FALSE)

  expect_equal(SummarizedExperiment::assay(out),
               SummarizedExperiment::assay(se))
})

test_that("bridge_normalise divides rather than subtracts when not log transformed", {
  se <- make_plex_se(base_quant, offsets = list("1" = rep(0, 4), "2" = rep(0, 4)))

  quant <- SummarizedExperiment::assay(se)
  quant[, 1:4] <- quant[, 1:4] * 2
  quant[, 5:8] <- quant[, 5:8] / 2
  SummarizedExperiment::assay(se) <- quant

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          on_log_scale = FALSE, verbose = FALSE)

  # both plexes are put onto the mean of the two bridge values
  expected <- base_quant * (2 + 0.5) / 2
  expect_equal(unname(SummarizedExperiment::assay(out)[, 1]), unname(expected))
  expect_equal(unname(SummarizedExperiment::assay(out)[, 5]), unname(expected))
})


# ---- Missing bridge values ------------------------------------------------

# f2 has no bridge value in plex 2, but does have sample data there
se_missing_bridge <- function(){
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))
  quant <- SummarizedExperiment::assay(se)
  quant["f2", "2_Bridge_1"] <- NA
  SummarizedExperiment::assay(se) <- quant
  se
}

test_that("on_missing='na' blanks the plex that has no usable bridge value", {
  se <- se_missing_bridge()

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          on_missing = "na", verbose = FALSE)

  quant <- SummarizedExperiment::assay(out)

  # plex 2 is unusable for f2, plex 1 is corrected onto plex 1's own bridge
  expect_true(all(is.na(quant["f2", se$Plex == "2"])))
  expect_equal(unname(quant["f2", "1_Sample_1"]), unname(base_quant[["f2"]] + 1))

  # other features are unaffected
  expect_equal(unname(quant["f1", ]), rep(unname(base_quant[["f1"]]), 8))
})

test_that("on_missing='drop' removes the feature entirely", {
  se <- se_missing_bridge()

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          on_missing = "drop", verbose = FALSE)

  expect_equal(rownames(out), c("f1", "f3", "f4"))
  expect_false(anyNA(SummarizedExperiment::assay(out)))
})

test_that("on_missing='ignore' leaves the plex uncorrected", {
  se <- se_missing_bridge()

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          on_missing = "ignore", verbose = FALSE)

  quant <- SummarizedExperiment::assay(out)

  # plex 2 keeps its original values for f2, so the plex offset survives
  expect_equal(unname(quant["f2", "2_Sample_1"]), unname(base_quant[["f2"]] - 1))
  # while plex 1 is corrected onto its own bridge, since it is the only
  # reference available for this feature
  expect_equal(unname(quant["f2", "1_Sample_1"]), unname(base_quant[["f2"]] + 1))
})

test_that("a feature absent from a plex is not treated as uncorrectable", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))
  quant <- SummarizedExperiment::assay(se)
  quant["f2", se$Plex == "2"] <- NA
  SummarizedExperiment::assay(se) <- quant

  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          on_missing = "drop", verbose = FALSE)

  # f2 has no data in plex 2 at all, so nothing was lost by not correcting it
  expect_true("f2" %in% rownames(out))
  expect_equal(unname(SummarizedExperiment::assay(out)["f2", "1_Sample_1"]),
               unname(base_quant[["f2"]] + 1))
})

test_that("min_bridge treats a partially observed reference as missing", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)),
                     n_bridge = list("1" = 2, "2" = 2),
                     n_sample = list("1" = 2, "2" = 2))
  quant <- SummarizedExperiment::assay(se)
  quant["f2", "2_Bridge_1"] <- NA
  SummarizedExperiment::assay(se) <- quant

  # one bridge value is enough by default
  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge", verbose = FALSE)
  expect_equal(unname(SummarizedExperiment::assay(out)["f2", "2_Sample_1"]),
               unname(base_quant[["f2"]]))

  # but not when two are demanded
  out <- bridge_normalise(se, plex_col = "Plex",
                          bridge_cols = se$Type == "Bridge",
                          min_bridge = 2, verbose = FALSE)
  expect_true(all(is.na(SummarizedExperiment::assay(out)["f2", se$Plex == "2"])))
})


# ---- Input validation ------------------------------------------------------

test_that("bridge_normalise accepts bridge_cols as column names", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))

  by_name <- bridge_normalise(se, plex_col = "Plex",
                              bridge_cols = c("1_Bridge_1", "2_Bridge_1"),
                              verbose = FALSE)
  by_logical <- bridge_normalise(se, plex_col = "Plex",
                                 bridge_cols = se$Type == "Bridge",
                                 verbose = FALSE)

  expect_equal(SummarizedExperiment::assay(by_name),
               SummarizedExperiment::assay(by_logical))
})

test_that("bridge_normalise rejects an unusable bridge specification", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))

  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = c(TRUE, FALSE)),
               "one value per sample")
  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = c("2_Bridge_1", "no_such_sample")),
               "not found in the samples")
  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = rep(FALSE, ncol(se))),
               "identifies no bridge channels")
  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = c(NA, rep(FALSE, ncol(se) - 1))),
               "contains missing values")
})

test_that("bridge_normalise requires a bridge in every plex", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))

  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = colnames(se) == "1_Bridge_1"),
               "no bridge channel in plex")
})

test_that("bridge_normalise rejects samples with no plex", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))
  se$Plex[3] <- NA

  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = se$Type == "Bridge"),
               "cannot be assigned to a plex")
})

test_that("bridge_normalise rejects an unknown reference plex or plex column", {
  se <- make_plex_se(base_quant,
                     offsets = list("1" = rep(1, 4), "2" = rep(-1, 4)))

  expect_error(bridge_normalise(se, plex_col = "Plex",
                                bridge_cols = se$Type == "Bridge",
                                reference = "3"),
               "is not present in column")
  expect_error(bridge_normalise(se, plex_col = "Batch",
                                bridge_cols = se$Type == "Bridge"),
               "is not in the colData")
})

test_that("bridge_normalise reports the correction and what was dropped", {
  se <- se_missing_bridge()

  expect_message(
    bridge_normalise(se, plex_col = "Plex",
                     bridge_cols = se$Type == "Bridge", on_missing = "na"),
    "median absolute correction")
  expect_message(
    bridge_normalise(se, plex_col = "Plex",
                     bridge_cols = se$Type == "Bridge", on_missing = "na"),
    "no usable bridge value")
})
