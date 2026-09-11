# PSM-level PD output for the phospho-enriched fraction of a TMTpro experiment

Proteome Discoverer PSM-level output for the phospho-enriched fraction
of a real TMTpro experiment in mouse fibroblasts, comparing a drug
treatment against a vehicle control at two timepoints (four replicates
of each combination), subsetted to a random selection of proteins plus a
sample of contaminant and non-unique-master-protein PSMs for use in the
`PTM site quantification` vignette. Includes the `ptmRS` columns needed
by
[`parse_PTM_scores_pd`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md).
Sample columns are named to match the rownames of
[`tmt_phospho_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_phospho_design.md).

[`psm_tmt_phospho_total`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_phospho_total.md)
is the matched total (non-enriched) fraction of the same labelled pool.

## Usage

``` r
psm_tmt_phospho
```

## Format

An object of class `data.frame` with 13354 rows and 80 columns.
