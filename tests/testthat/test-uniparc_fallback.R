httptest::.mockPaths(testthat::test_path("..", "mocks", "up"))
httptest2::.mockPaths(testthat::test_path("..", "mocks", "up"))

# Mock files under tests/mocks/up/rest.uniprot.org/{uniprotkb,uniparc}/ were
# recorded from the real UniProt API (record_uniparc_fixtures.R via httr v1,
# refreshed by record_uniparc_fallback_parallel_fixture.R via httr2 -- both
# resolve to the exact same "<path>.json" filenames for these plain GETs, so
# one fixture set serves get_uniparc_fallback_one() (httr v1) and
# get_uniparc_fallback() (httr2, see below) alike) against:
#   - A0A5G2QPJ4: a retired UniProtKB accession (deleted, ENSEMBL source
#     withdrawn) whose UniParc entry has exactly one cross-reference
#     (UniProtKB/TrEMBL) carrying both a gene name and a protein name.
#   - P04406 (GAPDH, human): a stable reviewed UniProtKB accession, used to
#     check the early-return path for entries that are still active. The
#     recorded response was trimmed to its top-level fields only, since
#     get_uniparc_fallback_one()/get_uniparc_fallback() only read `entryType`
#     for this case.
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

# get_uniparc_fallback() (unlike get_uniparc_fallback_one()) fetches both
# accessions' UniProtKB entries in parallel via httr2::req_perform_parallel(),
# then any UniParc entries needed, so it's mocked with httptest2 rather than
# httptest -- see the fixture note above for why the same files serve both.

test_that("get_uniparc_fallback combines results and drops unresolved accessions", {
  httptest2::with_mock_api({
    result <- get_uniparc_fallback(c("A0A5G2QPJ4", "P04406"))

    expect_equal(nrow(result), 1)
    expect_equal(result$UniprotID, "A0A5G2QPJ4")
  })
})

# get_uniprot_details() calls uniprotREST::uniprot_map() and, for any
# accession needing fallback, get_uniparc_fallback() -- both httr2-based, so
# a single httptest2::with_mock_api() covers the whole call. The uniprot_map()
# leg was recorded by data-raw/record_get_uniprot_details_fixture.R into this
# same tests/mocks/up/ directory as a single call against c("O76024",
# "A0A5G2QPJ4") -- both accessions must be queried together in one test to
# match that recording. The get_uniparc_fallback() leg reuses the fixtures
# already recorded for A0A5G2QPJ4.

test_that("get_uniprot_details annotates a live accession and falls back to UniParc for a retired one", {
  httptest2::with_mock_api({
    result <- get_uniprot_details(c("O76024", "A0A5G2QPJ4"))

    live_row <- result[result$UniprotID == "O76024", ]
    expect_equal(live_row$Annotation.Source, "UniProtKB (live)")
    expect_false(live_row$Gene.Names == "")

    fallback_row <- result[result$UniprotID == "A0A5G2QPJ4", ]
    expect_equal(fallback_row$Annotation.Source, "UniParc (retired UniProtKB entry)")
    expect_equal(fallback_row$Gene.Names, "CADM1")
  })
})
