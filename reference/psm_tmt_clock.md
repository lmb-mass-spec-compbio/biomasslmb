# PSM-level PD output for a TMT12plex Control vs Mutant experiment

Proteome Discoverer PSM-level output for a real TMT12plex experiment
comparing Control and Mutant samples (6 replicates each), subsetted to
~600 randomly-selected proteins plus a sample of contaminant and
non-unique-master-protein PSMs for use in the
`TMT workflow: PSM QC and protein summarisation` vignette. Sample
columns are named to match the rownames of
[`tmt_clock_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_clock_design.md).

## Usage

``` r
psm_tmt_clock
```

## Format

An object of class `data.frame` with 11281 rows and 67 columns.
