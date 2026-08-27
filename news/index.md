# Changelog

## biomasslmb (development version)

- New
  [`collapse_uniprot_details_multi_accession()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/collapse_uniprot_details_multi_accession.md)
  maps a
  [`get_uniprot_details()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
  result onto protein IDs that bundle several accessions (e.g. Proteome
  Discoverer’s `Master.Protein.Accessions`), collapsing back to one row
  per original ID.

- `Remotes:` temporarily points `uniprotREST` at
  `TomSmithCGAT/uniprotREST@check-complete-response` instead of
  `csdaw/uniprotREST`, to pick up the `check_complete` argument to
  `uniprot_map()` (detects UniProt ID mapping jobs that silently return
  a truncated result). This is an unmerged fork branch, not a released
  feature. **Switch `Remotes:` back to `csdaw/uniprotREST` once that PR
  is accepted upstream.**
