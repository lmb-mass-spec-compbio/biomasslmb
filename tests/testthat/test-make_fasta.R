httptest::.mockPaths(testthat::test_path("make_fasta"))

test_that("sub_crap inserts sequential cRAP numbers at the default start/width", {
  result <- sub_crap(c("|sometext|", "|moretext|"))
  expect_equal(result, c("|cRAP001|sometext|", "|cRAP002|moretext|"))
})

test_that("sub_crap respects a custom start", {
  result <- sub_crap(c("|sometext|", "|moretext|"), start = 88)
  expect_equal(result, c("|cRAP088|sometext|", "|cRAP089|moretext|"))
})

test_that("sub_crap respects a custom width", {
  result <- sub_crap(c("|sometext|", "|moretext|"), start = 1111, width = 4)
  expect_equal(result, c("|cRAP1111|sometext|", "|cRAP1112|moretext|"))
})

test_that("sub_crap only substitutes the first '|'", {
  result <- sub_crap("sp|XXXXXX|YYYY_YYYY Text goes here")
  expect_equal(result, "sp|cRAP001|XXXXXX|YYYY_YYYY Text goes here")
})

test_that("sub_crap errors when an element has no '|...|' pair", {
  expect_error(sub_crap("no pipes here"), "only works when each element")
})

# Mock file under make_fasta/ was built by hand (not via capture_requests(),
# which doesn't trace HEAD -- see data-raw/record_uniprot_release_fixture.R)
# from a real response to https://rest.uniprot.org/uniprot/P60709.fasta,
# which redirects to /uniprotkb/P60709.fasta.
#
# Uses httptest::with_mock_api() explicitly -- see test-uniparc_fallback.R
# for why (httptest2, loaded by other test files, would otherwise shadow it).

test_that("check_uniprot_release returns the release from the x-uniprot-release header", {
  httptest::with_mock_api({
    expect_equal(check_uniprot_release(), "2026_02")
  })
})
