# Records the httptest2 fixture used by tests/testthat/test-uniparc_fallback.R
# for get_uniprot_details(). Requires real network access to rest.uniprot.org.
#
# get_uniprot_details() calls uniprotREST::uniprot_map() (httr2-based, no
# `fields=` argument, so it uses uniprotREST's default field set) and, for any
# accession that comes back with a blank Gene.Names or Protein.names=="deleted",
# separately calls get_uniparc_fallback() (httr v1-based). This script only
# records the uniprot_map() leg -- the fallback leg reuses the existing
# httr fixtures under tests/mocks/up/, recorded by
# record_uniparc_fixtures.R, since it queries the same accession
# (A0A5G2QPJ4) the same way.
#
# Accessions used:
#   O76024: stable reviewed accession -> exercises the plain
#     'UniProtKB (live)' annotation path with no fallback needed.
#   A0A5G2QPJ4: retired accession (see record_uniparc_fixtures.R) -> comes
#     back from uniprot_map() with a blank Gene.Names/Protein.names, which is
#     what triggers get_uniprot_details()'s needs_fallback branch.

library(httptest2)

.mockPaths("tests/mocks/up")

capture_requests(simplify = FALSE, {
  uniprotREST::uniprot_map(
    ids = c("O76024", "A0A5G2QPJ4"),
    from = "UniProtKB_AC-ID",
    to = "UniProtKB",
    method = "stream",
    verbosity = 0
  )
})
