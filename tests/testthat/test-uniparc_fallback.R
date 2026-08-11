httptest::.mockPaths(testthat::test_path("uniparc_fallback"))

# Mock files under uniparc_fallback/ were recorded from the real UniProt API
# via httptest::capture_requests() against:
#   - A0A5G2QPJ4: a retired UniProtKB accession (deleted, ENSEMBL source
#     withdrawn) whose UniParc entry has exactly one cross-reference
#     (UniProtKB/TrEMBL) carrying both a gene name and a protein name.
#   - P04406 (GAPDH, human): a stable reviewed UniProtKB accession, used to
#     check the early-return path for entries that are still active. The
#     recorded response was trimmed to its top-level fields only, since
#     get_uniparc_fallback_one only reads `entryType` for this case.
#
# Uses httptest::with_mock_api() explicitly (rather than library(httptest) +
# unqualified with_mock_api()) because other test files load httptest2, whose
# with_mock_api() would otherwise shadow this one and silently let real
# network calls through instead of mocking them.

test_that("get_uniparc_fallback_one recovers annotation for a retired accession", {
  httptest::with_mock_api({
    result <- get_uniparc_fallback_one("A0A5G2QPJ4")

    expect_equal(result$UniprotID, "A0A5G2QPJ4")
    expect_equal(result$Gene.Names, "CADM1")
    expect_equal(result$Gene.Names.First, "CADM1")
    expect_equal(result$Protein.names, "Cell adhesion molecule 1")
    expect_equal(result$Annotation.Source, "UniParc (retired UniProtKB entry)")
  })
})

test_that("get_uniparc_fallback_one returns NULL for an active accession", {
  httptest::with_mock_api({
    result <- get_uniparc_fallback_one("P04406")
    expect_null(result)
  })
})

test_that("get_uniparc_fallback combines results and drops unresolved accessions", {
  httptest::with_mock_api({
    result <- get_uniparc_fallback(c("A0A5G2QPJ4", "P04406"))

    expect_equal(nrow(result), 1)
    expect_equal(result$UniprotID, "A0A5G2QPJ4")
  })
})
