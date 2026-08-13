# Records the httptest fixture used by tests/testthat/test-make_fasta.R for
# make_fasta() and its internal check_uniprot_job() polling helper. Requires
# real network access to rest.uniprot.org.
#
# make_fasta() drives a 3-step UniProt ID mapping job via httr (v1), not
# httr2:
#   1. POST  /idmapping/run                                  (submit job)
#   2. GET   /idmapping/status/{jobId}                        (poll until done)
#   3. GET   /idmapping/uniprotkb/results/stream/{jobId}?format=fasta
#
# This mirrors the flow already recorded for tests/mocks/go and tests/mocks/tm,
# but those use uniprotREST's httr2 client -- make_fasta()'s own POST/GET
# calls are plain httr, so they need their own httr-format fixtures (same
# httptest package/format as record_uniparc_fixtures.R).
#
# Uses 2 small, stable accessions so the job completes quickly.
# check_uniprot_job() polls every 5s until the job finishes, so this script
# will take at least one 5s sleep (usually just one poll) to run.

library(httptest)
library(httr)

capture_requests(path = "tests/mocks/fa", {
  payload <- list(ids = "P60709 P04406", from = "UniProtKB_AC-ID", to = "UniProtKB")

  post_payload <- httr::POST(
    url = "https://rest.uniprot.org/idmapping/run",
    body = payload, encode = "multipart", httr::accept_json()
  )
  post_resp <- httr::content(post_payload, as = "parsed")
  job_id <- post_resp[["jobId"]]

  # Mirrors check_uniprot_job()'s own polling loop.
  repeat {
    status <- httr::content(
      httr::GET(sprintf("https://rest.uniprot.org/idmapping/status/%s", job_id), httr::accept_json()),
      as = "parsed"
    )
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) break
    Sys.sleep(5)
  }

  httr::GET(
    url = sprintf("https://rest.uniprot.org/idmapping/uniprotkb/results/stream/%s?format=fasta", job_id),
    config = httr::write_disk(tempfile(fileext = ".fasta"), overwrite = TRUE)
  )
})
