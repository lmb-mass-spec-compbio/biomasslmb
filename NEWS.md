# biomasslmb (development version)

* New `lfq_qf_turboid` dataset: an LFQ-DDA TurboID pulldown (three biotin
  against three control), processed to protein level in Part B of the
  `Enrichment designs` vignette, which now carries an enrichment experiment
  through to a tested list of candidate interactors.

* **`tmt_qf_mq` has changed.** Its `protein` assay is now masked with
  `get_protein_no_quant_mask()`/`mask_protein_level_quant()`, so that every
  value rests on at least two PSMs, as the LFQ QC vignettes already did. This
  removes 130 proteins that were quantified from a single PSM per channel,
  leaving 442. Analyses depending on the previous 572 rows will need re-running.

* `get_crap_fasta_accessions()` is renamed to `get_contaminant_fasta_accessions()`,
  for consistency with `remove_contaminant()`/`remove_contaminant_mq()` and because
  the function isn't specific to the cRAP database. The old name is kept as a
  wrapper that forwards to the new function and emits a deprecation message.

* New `collapse_uniprot_details_multi_accession()` maps a
  `get_uniprot_details()` result onto protein IDs that bundle several
  accessions (e.g. Proteome Discoverer's `Master.Protein.Accessions`),
  collapsing back to one row per original ID.

* `Remotes:` temporarily points `uniprotREST` at
  `TomSmithCGAT/uniprotREST@check-complete-response` instead of
  `csdaw/uniprotREST`, to pick up the `check_complete` argument to
  `uniprot_map()` (detects UniProt ID mapping jobs that silently return a
  truncated result). This is an unmerged fork branch, not a released
  feature. **Switch `Remotes:` back to `csdaw/uniprotREST` once that PR is
  accepted upstream.**
