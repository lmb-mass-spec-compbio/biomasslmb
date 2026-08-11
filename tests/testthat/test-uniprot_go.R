httptest2::.mockPaths(testthat::test_path("uniprot_go"))

# Mock files under uniprot_go/ were recorded from the real UniProt ID mapping
# API via httptest2::capture_requests(simplify = FALSE) against
# c("O76024", "Q03135"), expand_terms = FALSE. simplify = FALSE is required
# for the same reason as tests/testthat/test-query_protein_tm_topology.R.
#
# Uses httptest2::with_mock_api() explicitly -- see test-uniparc_fallback.R
# for why (httptest, loaded by other test files, would otherwise shadow it).

test_that("get_go_terms returns one row per GO annotation for real accessions", {
  httptest2::with_mock_api({
    result <- get_go_terms(c("O76024", "Q03135"), expand_terms = FALSE)

    expect_equal(colnames(result), c("UNIPROTKB", "GO_desc", "GO.ID"))
    expect_setequal(unique(result$UNIPROTKB), c("O76024", "Q03135"))
    expect_true(all(grepl("^GO:[0-9]{7}$", result$GO.ID)))
    expect_true("GO:0030425" %in% result$GO.ID[result$UNIPROTKB == "O76024"])
  })
})
