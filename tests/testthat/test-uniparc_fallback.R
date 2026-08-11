library(httptest)

.mockPaths(test_path("uniparc_fallback"))

# Mock files under uniparc_fallback/ were recorded from the real UniProt API
# via httptest::capture_requests() against:
#   - A0A5G2QPJ4: a retired UniProtKB accession (deleted, ENSEMBL source
#     withdrawn) whose UniParc entry has exactly one cross-reference
#     (UniProtKB/TrEMBL) carrying both a gene name and a protein name.
#   - P04406 (GAPDH, human): a stable reviewed UniProtKB accession, used to
#     check the early-return path for entries that are still active. The
#     recorded response was trimmed to its top-level fields only, since
#     get_uniparc_fallback_one only reads `entryType` for this case.

test_that("get_uniparc_fallback_one recovers annotation for a retired accession", {
  with_mock_api({
    result <- get_uniparc_fallback_one("A0A5G2QPJ4")

    expect_equal(result$UniprotID, "A0A5G2QPJ4")
    expect_equal(result$Gene.Names, "CADM1")
    expect_equal(result$Gene.Names.First, "CADM1")
    expect_equal(result$Protein.names, "Cell adhesion molecule 1")
    expect_equal(result$Annotation.Source, "UniParc (retired UniProtKB entry)")
  })
})

test_that("get_uniparc_fallback_one returns NULL for an active accession", {
  with_mock_api({
    result <- get_uniparc_fallback_one("P04406")
    expect_null(result)
  })
})

test_that("get_uniparc_fallback combines results and drops unresolved accessions", {
  with_mock_api({
    result <- get_uniparc_fallback(c("A0A5G2QPJ4", "P04406"))

    expect_equal(nrow(result), 1)
    expect_equal(result$UniprotID, "A0A5G2QPJ4")
  })
})
