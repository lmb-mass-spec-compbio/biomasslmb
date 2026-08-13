# Recorded HTTP fixtures

Mock responses replayed by `httptest::with_mock_api()` (httr v1) and
`httptest2::with_mock_api()` (httr2). Test files point at these with
`.mockPaths(testthat::test_path("..", "mocks", "<dir>"))`; the recording
scripts that produced them live in `data-raw/record_*.R`.

| dir  | used by                                                | recorded by                          |
|------|--------------------------------------------------------|--------------------------------------|
| `fa` | `test-make_fasta.R`                                    | `record_make_fasta_fixture.R`, `record_uniprot_release_fixture.R` |
| `go` | `test-uniprot_go.R`                                    | recorded alongside the `tm` fixtures |
| `tm` | `test-query_protein_tm_topology.R`                     | recorded alongside the `go` fixtures |
| `up` | `test-uniparc_fallback.R`                              | `record_uniparc_fixtures.R`, `record_get_uniprot_details_fixture.R` |

## Why the directory names are two characters

Both mocking libraries derive the fixture's path from the request URL, and
`R CMD build` warns about any path over 100 bytes once `biomasslmb/` is
prepended. The deepest UniProt endpoint used here already accounts for 83 of
those bytes:

```
biomasslmb/tests/mocks/fa/rest.uniprot.org/idmapping/uniprotkb/results/stream/81nmQqWHyV-cb5419.txt
^-- 11 --^ ^-- 12 --^ ^2^ ^------------------- 52 -------------------------^ ^-------- 21 -------^
```

That leaves almost no slack — the longest path above is 99 bytes. These
fixtures therefore sit in `tests/mocks/` rather than `tests/testthat/`
(`testthat/` alone would push them over), and the directory names are kept to
two characters.

If a new fixture comes from a deeper endpoint or a longer job ID, check
`R CMD build` output rather than assuming it fits.
