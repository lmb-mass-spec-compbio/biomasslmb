# PSM-level MaxQuant output for a TMT18plex IP vs control experiment

MaxQuant `evidence.txt` PSM-level output for a real TMT18plex
immunoprecipitation experiment comparing a bait pulldown (`IP`) with a
control pulldown (`Control`), 6 replicates each, subsetted to ~600
randomly-selected proteins plus a sample of contaminant and
non-unique-master-protein PSMs for use as the MaxQuant worked example in
the `Enrichment designs: IP, BioID and TurboID` vignette. Sample columns
are named to match the rownames of
[`tmt_per2_mq_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_per2_mq_design.md).
The original condition labels and the `Raw.file` column (which
identified the specific source project) have been anonymised to generic
values, and PSMs identifying the bait construct itself have been
excluded.

## Usage

``` r
psm_tmt_per2_mq
```

## Format

An object of class `data.frame` with 4960 rows and 76 columns.
