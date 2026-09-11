# PSM-level PD output for a two-plex TMTpro 18plex experiment

Proteome Discoverer PSM-level output for a real experiment spread over
two TMTpro 18plex plexes, each carrying a pooled bridge channel. The
design is a 3 x 3 factorial: three cell lines (a wild type and two
knockouts) crossed with three drug treatment timepoints, three
replicates each. Subsetted to ~300 randomly-selected proteins quantified
in both plexes, plus a sample of contaminant and
non-unique-master-protein PSMs, for use in the
`TMT workflow: multiple plexes and bridge correction` vignette.

The `Plex` column identifies which plex each PSM came from, and the
abundance columns are named by TMT tag, since the same tag denotes a
different sample in each plex.
[`tmt_2plex_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_2plex_design.md)
maps tag and plex to sample.

## Usage

``` r
psm_tmt_2plex
```

## Format

An object of class `data.frame` with 11619 rows and 73 columns.
