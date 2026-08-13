httptest2::.mockPaths(testthat::test_path("..", "mocks", "tm"))

# Mock files under tests/mocks/tm/ were recorded from the real UniProt ID
# mapping API via httptest2::capture_requests(simplify = FALSE) against
# c("O76024", "Q03135"). simplify = FALSE is required here (rather than the
# plain-text default) because uniprot_map()'s internal polling logic reads the
# job status response's redirected `url` and its `x-total-results` header,
# not just its body, and those are only preserved in the full serialized
# response object.
#
# Uses httptest2::with_mock_api() explicitly (rather than library(httptest2)
# + unqualified with_mock_api()) because other test files load httptest,
# whose with_mock_api() would otherwise shadow this one and silently let real
# network calls through instead of mocking them.

test_that("query_protein_tm_topology returns TM/topology columns for real accessions", {
  httptest2::with_mock_api({
    result <- query_protein_tm_topology(c("O76024", "Q03135"))

    expect_equal(result$UniprotID, c("O76024", "Q03135"))
    expect_equal(result$Length, c(890, 178))
    expect_true(grepl("^TRANSMEM 314\\.\\.334", result$Transmembrane[1]))
    expect_equal(result$Transmembrane[2], "")
    expect_true(grepl("^INTRAMEM 105\\.\\.125", result$Intramembrane[2]))
  })
})

test_that("get_protein_tm_topology chains query_protein_tm_topology with the TM/topology parsers", {
  httptest2::with_mock_api({
    result <- get_protein_tm_topology(c("O76024", "Q03135"))

    expect_equal(result$UniprotID, c("O76024", "Q03135"))
    expect_equal(result$n_tms, c(11, 0))
    expect_equal(result$tm_start[[1]], "314;340;402;427;465;496;529;563;589;632;870")
  })
})
