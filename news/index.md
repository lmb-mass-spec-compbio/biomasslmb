# Changelog

## biomasslmb 0.1.0

### Breaking changes

- **The `filter_features_*()` functions now take two separate arguments
  where they previously took one.** `unique_master` asks whether the
  search engine resolved the feature to a single protein accession;
  `proteotypic` asks whether the peptide sequence occurs in only one
  protein. These are different questions, and `unique_master` had come
  to mean whichever of them the format happened to expose:
  `Number.of.Protein.Groups` for Proteome Discoverer, the accession
  count for DIA-NN and Spectronaut, but `Unique..Proteins.` — a
  proteotypic test — for MaxQuant.

  [`filter_features_pd_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_pd_dda.md),
  [`filter_features_diann()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_diann.md)
  and
  [`filter_features_sn()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_sn.md)
  are unchanged for existing calls, and gain `proteotypic` (default
  `FALSE`), which tests `Number.of.Proteins`, `Proteotypic` and
  `PEP.IsProteotypic` respectively.

  **[`filter_features_mq_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_mq_dda.md)
  changes behaviour, and its `unique_master` default changes with it.**
  MaxQuant assigns every feature a single razor protein and so never
  reports a tie between protein groups; there is nothing for
  `unique_master` to filter. It now defaults to `FALSE` and raises an
  error if passed `TRUE`, where it previously defaulted to `TRUE` and
  silently applied a proteotypic filter. A call that named the argument
  fails loudly; a call that relied on the default keeps features it used
  to remove, and nothing announces it. Pass `proteotypic = TRUE` to
  reproduce the previous behaviour, which it does exactly.

- **`dia_qf`, `lfq_qf` and `tmt_qf` hold different data under the same
  names.** Each is rebuilt by the rewritten workflow vignette it comes
  from, and analyses depending on any of them need re-running. Where an
  assay name is unchanged, the old code runs and quietly returns
  different numbers; where one is renamed, it errors.

  - `dia_qf` is a different experiment: a 31-sample plasma DIA dataset
    (`monkeypox_plasma_proteomes.parquet`) summarised to 239 proteins,
    where it previously held 407 proteins across 6 samples. Its assays
    are renamed with it — `peptides_fdr_cntrl` to `precursors`,
    `peptides_filtered_forRobust` to `peptides_for_summarisation`, and a
    `peptides_filtered_norm` assay is added — so code indexing an assay
    by name errors rather than quietly returning different numbers.
  - `lfq_qf` is rebuilt from a replaced, larger
    `lfq_dda_pd_PeptideGroups.txt`: 230 proteins across the same 6
    samples, where it previously held 339. Its `peptides_norm` assay is
    renamed `peptides_filtered_norm`.
  - `tmt_qf` is built from `psm_tmt_clock`, a TMT12plex experiment,
    rather than from the TMT10-plex `psm_tmt_total`: 402 proteins across
    12 channels, where it previously held 891 across 10. Its assay names
    are unchanged.

  `psm_tmt_total` itself is unchanged and still backs the function
  examples. The vignettes use `psm_tmt_clock`, which carries an
  experimental design table with it.

- [`get_crap_fasta_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_crap_fasta_accessions.md)
  is renamed to
  [`get_contaminant_fasta_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_contaminant_fasta_accessions.md),
  for consistency with `remove_contaminant()`/`remove_contaminant_mq()`
  and because the function isn’t specific to the cRAP database. The old
  name is kept as a wrapper that forwards to the new function and emits
  a deprecation message.

### New features

- New
  [`bridge_normalise()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/bridge_normalise.md)
  puts TMT plexes onto a common scale using their pooled bridge
  channels. The correction is per feature and per plex, because that is
  what a plex effect is — a feature is quantified from a different set
  of PSMs in each plex — which is what distinguishes it from
  [`center_normalise_to_ref()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/center_normalise_to_ref.md)
  and its single offset per sample. Features with no usable bridge value
  in a plex cannot be placed on the scale at all, and `on_missing`
  decides whether they are set to `NA`, dropped, or left uncorrected.

- New
  [`sync_coldata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/sync_coldata.md)
  copies the object-level `colData` onto named assays, matching on
  sample name. It replaces the `colData(obj[[i]]) <- colData(obj)`
  idiom, which appeared 32 times across the vignettes.

- New
  [`collapse_uniprot_details_multi_accession()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/collapse_uniprot_details_multi_accession.md)
  maps a
  [`get_uniprot_details()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
  result onto protein IDs that bundle several accessions (e.g. Proteome
  Discoverer’s `Master.Protein.Accessions`), collapsing back to one row
  per original ID.

- [`get_n_feature_per_prot()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_n_feature_per_prot.md),
  [`message_parse()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/message_parse.md)
  and
  [`qfeatures_long()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/qfeatures_long.md)
  are now exported. Each was already used internally and is called
  directly by the vignettes:
  [`get_n_feature_per_prot()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_n_feature_per_prot.md)
  tallies the finite values per protein per sample that summarisation
  needs,
  [`message_parse()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/message_parse.md)
  reports a feature and master protein count in the same format the
  `filter_features_*()` functions use, so a count can be reported at a
  step where no filter function ran, and
  [`qfeatures_long()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/qfeatures_long.md)
  dispatches to the right QFeatures long-format function for the
  installed Bioconductor version.
  [`message_parse()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/message_parse.md)’s
  third argument is named `note`.

- [`plot_missing_upset()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_missing_upset.md)
  gains a `...` argument, passed onto
  [`naniar::gg_miss_upset()`](https://naniar.njtierney.com/reference/gg_miss_upset.html)
  and from there onto
  [`UpSetR::upset()`](https://rdrr.io/pkg/UpSetR/man/upset.html).
  Arguments given this way take precedence over the ones the function
  sets itself, so `nintersects`, `sets` and `keep.order` are now
  defaults rather than fixed. Samples are also ordered by name before
  plotting.

### Bug fixes

- [`filter_by_protein_fdr()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_by_protein_fdr.md)
  honours `protein_FDR_col`. It took the argument but read
  `Protein.FDR.Confidence.Combined` regardless, so a column named
  anything else was silently ignored.

### New example data

Each PSM-level dataset is a real facility experiment, subsetted to a few
hundred proteins and anonymised, and comes with a design table whose
rownames match its sample columns.

- `psm_tmt_clock` with `tmt_clock_design`: Proteome Discoverer output
  for a TMT12plex control-against-mutant comparison, six replicates
  each.

- `psm_tmt_2plex` with `tmt_2plex_design`: one experiment spread over
  two TMTpro 18plex plexes, each carrying a pooled bridge channel, in a
  3 x 3 factorial of cell line against treatment timepoint. The
  abundance columns are named by TMT tag, since the same tag means a
  different sample in each plex.

- `psm_tmt_factorial` with `tmt_factorial_design`: MaxQuant
  `evidence.txt` output for a TMT18plex whole-proteome experiment
  crossing three cell lines with vehicle or compound treatment.

- `psm_tmt_per2_mq` with `tmt_per2_mq_design`: MaxQuant `evidence.txt`
  output for a TMT18plex immunoprecipitation, bait pulldown against
  control, six replicates each.

- `psm_tmt_phospho` and `psm_tmt_phospho_total` with
  `tmt_phospho_design`: the phospho-enriched and matched total fractions
  of one labelled TMTpro pool, a drug treatment against vehicle at two
  timepoints. The enriched fraction carries the `ptmRS` columns
  [`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
  needs.

- `tmt_qf_mq` holds `psm_tmt_per2_mq` processed to protein level. Its
  `protein` assay is masked with
  [`get_protein_no_quant_mask()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_protein_no_quant_mask.md)/[`mask_protein_level_quant()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mask_protein_level_quant.md),
  so that every value rests on at least two PSMs, as the LFQ QC
  vignettes already did; 130 of the 572 proteins were quantified from a
  single PSM per channel and are masked out, leaving 442.

- `tmt_qf_factorial` holds the protein-level abundances summarised from
  `psm_tmt_factorial`, read back in by the statistical testing article.
  Only the `protein` assay is kept — retaining the PSM-level assays
  would make the object an order of magnitude larger.

- `lfq_qf_turboid` is an LFQ-DDA TurboID pulldown (three biotin against
  three control), processed to protein level in Part B of the
  `Enrichment designs` vignette, which now carries an enrichment
  experiment through to a tested list of candidate interactors.

### Documentation

- The vignettes go from four articles to nineteen, covering the four
  acquisition and search-engine routes through QC and summarisation,
  then exploration and statistical testing, enrichment designs, PTM site
  quantification, multiplexed TMT, missing values, functional
  enrichment, iBAQ, and a set of cautionary `gotcha_*` articles on
  contaminants and FDR, peptide-to-protein inference, and technical
  replicates.

- New `Getting started` and `Working with QFeatures objects` articles
  cover the object model and the design-table contract, with a complete
  forty-line worked analysis at the top of the former.

- **Two vignettes are renamed, so their URLs change.**
  `LFQ_DIA_Peptide_QC` becomes `LFQ_DIA_Precursor_QC_Summarisation`, and
  `TMT_summarisation_methods` becomes `summarisation_methods`, which now
  covers LFQ as well as TMT.

- `_pkgdown.yml` groups the function reference by task rather than
  listing it alphabetically, and orders the articles for the sidebar.
  Nineteen of the commonly used functions gain runnable examples.

### Packaging

- The maintainer address is now `tsmith@mrclmb.ac.uk`.

- New `Suggests`: `arrow`, `fgsea`, `ggrepel`, `hexbin`, `imputeLCMD`,
  `limma`, `limpa`, `openxlsx`, `readxl` and `xml2`, all used by the new
  vignettes. None is needed to use the package itself.

- `Remotes:` temporarily points `uniprotREST` at
  `TomSmithCGAT/uniprotREST@check-complete-response` instead of
  `csdaw/uniprotREST`, to pick up the `check_complete` argument to
  `uniprot_map()` (detects UniProt ID mapping jobs that silently return
  a truncated result). This is an unmerged fork branch, not a released
  feature; it is proposed upstream as
  [csdaw/uniprotREST#7](https://github.com/csdaw/uniprotREST/pull/7).
  **Switch `Remotes:` back to `csdaw/uniprotREST` once that PR is
  accepted upstream.**
