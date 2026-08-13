# Re-records the httptest fixtures used by tests/testthat/test-uniparc_fallback.R.
# Only needs to be re-run if the UniProt/UniParc JSON schema changes, or if the
# test accessions themselves change (e.g. A0A5G2QPJ4 is re-merged into an
# active entry). Requires real network access to rest.uniprot.org.
#
# A0A5G2QPJ4: retired UniProtKB accession (deleted, ENSEMBL source withdrawn),
#   whose UniParc entry currently has exactly one cross-reference with both a
#   gene name and a protein name.
# P04406 (GAPDH, human): stable reviewed UniProtKB accession, used to test the
#   early-return path for entries that are still active.

library(httptest)
library(httr)

capture_requests(path = "tests/mocks/up", {
  GET("https://rest.uniprot.org/uniprotkb/A0A5G2QPJ4.json")
  GET("https://rest.uniprot.org/uniparc/UPI0020B0616E.json")
  GET("https://rest.uniprot.org/uniprotkb/P04406.json")
})

# The P04406 uniprotkb response is ~300kB but get_uniparc_fallback_one() only
# reads `entryType` from it (it returns early once the entry is active), so
# after recording, trim it down to just the top-level fields:
#
# tests/mocks/up/rest.uniprot.org/uniprotkb/P04406.json.json
# {
#     "entryType": "UniProtKB reviewed (Swiss-Prot)",
#     "primaryAccession": "P04406",
#     "uniProtkbId": "G3P_HUMAN"
# }
