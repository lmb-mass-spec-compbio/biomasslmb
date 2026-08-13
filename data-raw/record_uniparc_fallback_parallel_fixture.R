# Records the httptest2 fixtures used by tests/testthat/test-uniparc_fallback.R
# for get_uniparc_fallback() (the batch/parallel path -- httr2-based). These
# are plain 200-JSON GETs with no query string, so httptest2's default
# simplify = TRUE writes them to the exact same "<path>.json" filenames that
# record_uniparc_fixtures.R already wrote for get_uniparc_fallback_one()'s
# httr v1 fixtures (both packages resolve mocks by URL the same way for this
# case) -- so this script reuses/refreshes those files rather than creating a
# separate httr2-only copy. Requires real network access to rest.uniprot.org.
#
# Requests are issued one at a time via httr2::req_perform() rather than by
# calling get_uniparc_fallback() itself, because httptest2::capture_requests()
# works by tracing req_perform() -- it never fires for the internal per-request
# calls httr2::req_perform_parallel() makes, so nothing gets recorded that
# way. Replay is unaffected: get_uniparc_fallback()'s test only exercises the
# mocked *response*, which httptest2 looks up purely by URL, regardless of
# whether the original recording used req_perform() or req_perform_parallel().
#
# Same two accessions as record_uniparc_fixtures.R, for the same reasons:
#   A0A5G2QPJ4: retired UniProtKB accession -> exercises the UniParc lookup.
#   P04406 (GAPDH, human): stable reviewed accession -> exercises the
#     early-drop path (no UniParc lookup needed). The recorded response is
#     ~300kB; get_uniparc_fallback() only reads `entryType` for this case, so
#     after recording, trim it down to just the top-level fields (see
#     record_uniparc_fixtures.R for the trimmed content).

library(httptest2)

.mockPaths("tests/mocks/up")

capture_requests({
  urls <- c(
    "https://rest.uniprot.org/uniprotkb/A0A5G2QPJ4.json",
    "https://rest.uniprot.org/uniprotkb/P04406.json",
    "https://rest.uniprot.org/uniparc/UPI0020B0616E.json"
  )
  for (url in urls) {
    httr2::req_perform(uniprotREST::uniprot_request(url, rate = 10))
  }
})
