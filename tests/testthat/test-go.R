# These tests use real, stable GO.db terms rather than synthetic ones, since
# determine_ancestor_function() etc. dispatch on real GO.db mapping objects.
# GO:0006915 = apoptotic process (BP); ancestors include GO:0012501
# (programmed cell death), GO:0008219 (cell death), GO:0009987 (cellular
# process), GO:0008150 (biological_process root).
# GO:0007049 = cell cycle (BP), verified unrelated (no ancestor/offspring
# relationship) to the apoptosis terms above.

test_that("determine_ancestor_function/determine_offspring_function dispatch on ontology", {
  expect_identical(determine_ancestor_function("x", "MF"), GO.db::GOMFANCESTOR)
  expect_identical(determine_ancestor_function("x", "BP"), GO.db::GOBPANCESTOR)
  expect_identical(determine_ancestor_function("x", "CC"), GO.db::GOCCANCESTOR)
  expect_true(is.na(determine_ancestor_function("x", NA)))
  expect_true(is.na(determine_ancestor_function("x", "weird")))

  expect_identical(determine_offspring_function("x", "BP"), GO.db::GOBPOFFSPRING)
  expect_true(is.na(determine_offspring_function("x", NA)))
})

test_that("get_all_mappings returns real ancestor terms for a known GO ID", {
  res <- suppressMessages(get_all_mappings(
    "GO:0006915", c("GO:0006915" = "BP"), verbose = FALSE, direction = "ancestor"
  ))

  expect_setequal(res[["GO:0006915"]],
                   c("GO:0008150", "GO:0008219", "GO:0009987", "GO:0012501", "all"))
})

test_that("expand_terms unions a protein's GO terms with their ancestors, dropping 'all'", {
  ancestor_map <- list("GO:0006915" = c("GO:0008150", "GO:0012501", "all"))
  go_df <- data.frame(GOID = "GO:0006915", stringsAsFactors = FALSE)

  out <- expand_terms(go_df, "GOID", ancestor_map)

  expect_setequal(out$GOID, c("GO:0006915", "GO:0008150", "GO:0012501"))
  expect_false("all" %in% out$GOID)
})

test_that("get_ancestor_go expands a protein's GO terms to include ancestors, terms and ontology", {
  go_df <- data.frame(UNIPROTKB = "P1", GO.ID = "GO:0006915", stringsAsFactors = FALSE)

  out <- suppressMessages(get_ancestor_go(go_df, verbose = FALSE))

  expect_setequal(out$GO.ID,
                   c("GO:0006915", "GO:0008150", "GO:0008219", "GO:0009987", "GO:0012501"))
  expect_true(all(out$ONTOLOGY == "BP"))
  expect_equal(out$TERM[out$GO.ID == "GO:0006915"], "apoptotic process")
})

test_that("remove_redundant_go drops an ancestor term in favour of its more significant descendant", {
  obj <- data.frame(
    category = c("GO:0006915", "GO:0012501", "GO:0007049"),
    over_represented_pvalue = c(0.001, 0.01, 0.5),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(suppressWarnings(remove_redundant_go(obj)))

  # GO:0012501 is an ancestor of GO:0006915 -> redundant, dropped.
  # GO:0007049 is unrelated -> retained.
  expect_setequal(out$category, c("GO:0006915", "GO:0007049"))
})

test_that("estimate_overrep computes overrep and weighted adj_overrep scores", {
  pwf <- data.frame(
    DEgenes = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    bias.data = rep(1, 6),
    pwf = c(0.9, 0.9, 0.1, 0.1, 0.1, 0.1),
    row.names = paste0("G", 1:6)
  )
  gene2cat <- data.frame(
    gene = paste0("G", 1:6),
    term = c("GO:X", "GO:X", "GO:X", "GO:Y", "GO:Y", "GO:Y"),
    stringsAsFactors = FALSE
  )
  obj <- data.frame(
    category = c("GO:X", "GO:Y"),
    over_represented_pvalue = c(0.01, 0.5),
    numDEInCat = c(2, 0),
    numInCat = c(3, 3),
    stringsAsFactors = FALSE
  )

  out <- estimate_overrep(obj, pwf, gene2cat)

  expect_equal(out$overrep[out$category == "GO:X"], 2)
  expect_equal(out$adj_overrep[out$category == "GO:X"], 0.3157895, tolerance = 1e-6)
  expect_equal(out$adj_overrep[out$category == "GO:Y"], 0)
})

test_that("add_independent_filtering_padj adds a padj_if column", {
  set.seed(42)
  n <- 200
  filter_vals <- sample(2:50, n, replace = TRUE)
  pvals <- runif(n)
  pvals[filter_vals > 30] <- pvals[filter_vals > 30] * 0.01
  obj <- data.frame(over_represented_pvalue = pvals, numInCat = filter_vals)

  out <- add_independent_filtering_padj(obj, plot = FALSE)

  expect_true("padj_if" %in% colnames(out))
  expect_equal(nrow(out), n)
  expect_true(all(out$padj_if >= 0 & out$padj_if <= 1))
})

test_that("get_enriched_go enriches a GO term using a supplied gene2cat and adjusts p-values", {
  suppressPackageStartupMessages(library(goseq))
  genes <- paste0("G", 1:30)
  pwf <- data.frame(
    DEgenes = c(rep(TRUE, 10), rep(FALSE, 20)),
    bias.data = rep(1, 30),
    pwf = rep(0.5, 30),
    row.names = genes
  )
  gene2cat <- data.frame(
    gene = genes,
    term = c(rep("GO:0006915", 10), rep("GO:0007049", 10), rep("GO:0043065", 10)),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(get_enriched_go(pwf, gene2cat, use_genes_without_cat = TRUE))

  # GO:0043065 and GO:0007049 have 0 DE genes -> dropped by filter_no_DE
  expect_equal(out$category, "GO:0006915")
  expect_equal(out$numDEInCat, 10)
  expect_true("over_represented_adj_pval" %in% colnames(out))
  expect_true("under_represented_adj_pval" %in% colnames(out))
  expect_equal(out$term_short, "apoptotic process")
})

test_that("get_enriched_go with filter_no_DE=FALSE keeps zero-DE categories", {
  suppressPackageStartupMessages(library(goseq))
  genes <- paste0("G", 1:30)
  pwf <- data.frame(
    DEgenes = c(rep(TRUE, 10), rep(FALSE, 20)),
    bias.data = rep(1, 30),
    pwf = rep(0.5, 30),
    row.names = genes
  )
  gene2cat <- data.frame(
    gene = genes,
    term = c(rep("GO:0006915", 10), rep("GO:0007049", 10), rep("GO:0043065", 10)),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(get_enriched_go(pwf, gene2cat, use_genes_without_cat = TRUE, filter_no_DE = FALSE))

  expect_equal(nrow(out), 3)
})

test_that("plot_go returns a ggplot object and errors if required columns are missing", {
  df <- data.frame(
    term_short = c("A GO term", "Another GO term", "One more GO term"),
    ontology = c("BP", "MF", "CC"),
    over_represented_adj_pval = c(0.0001, 1, 0.01),
    adj_overrep = c(15, 3, 1),
    numDEInCat = c(304, 22, 78)
  )

  p <- plot_go(df)
  expect_s3_class(p, "ggplot")

  expect_error(plot_go(df[, -1]), "missing the following required columns")
})
