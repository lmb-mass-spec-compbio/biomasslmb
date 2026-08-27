# biomasslmb (development version)

* `get_crap_fasta_accessions()` is renamed to `get_contaminant_fasta_accessions()`,
  for consistency with `remove_contaminant()`/`remove_contaminant_mq()` and because
  the function isn't specific to the cRAP database. The old name is kept as a
  wrapper that forwards to the new function and emits a deprecation message.

* `Remotes:` temporarily points `uniprotREST` at
  `TomSmithCGAT/uniprotREST@check-complete-response` instead of
  `csdaw/uniprotREST`, to pick up the `check_complete` argument to
  `uniprot_map()` (detects UniProt ID mapping jobs that silently return a
  truncated result). This is an unmerged fork branch, not a released
  feature. **Switch `Remotes:` back to `csdaw/uniprotREST` once that PR is
  accepted upstream.**
