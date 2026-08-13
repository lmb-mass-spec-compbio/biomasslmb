# Re-records the httptest fixture used by tests/testthat/test-make_fasta.R
# for check_uniprot_release(). Requires real network access to
# rest.uniprot.org.
#
# check_uniprot_release() uses httr::HEAD(), which redirects
# (/uniprot/... -> /uniprotkb/...) and reads the x-uniprot-release header
# from the final response. httptest's capture_requests() only traces
# GET/POST/PUT/PATCH/DELETE/VERB/RETRY -- not HEAD, since httr::HEAD() calls
# request_perform() directly rather than going through VERB() -- so the
# fixture has to be built by hand from a real response instead of via
# capture_requests(). Replay via with_mock_api() works fine for HEAD, since
# it patches request_perform() directly.

library(httr)

r <- HEAD("https://rest.uniprot.org/uniprot/P60709.fasta")
r$handle <- NULL
r$request <- NULL

dst <- "tests/mocks/fa/rest.uniprot.org/uniprot/P60709.fasta-HEAD.R"
dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
f <- file(dst, "wb", encoding = "UTF-8")
dput(r, file = f)
close(f)
