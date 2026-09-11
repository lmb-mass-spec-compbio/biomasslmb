# Experimental design for `psm_tmt_2plex`

A `data.frame` with one row per labelled channel, giving the Genotype,
Time and Replicate for each sample in
[`psm_tmt_2plex`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_2plex.md).
Channels that were not labelled in a given plex are absent. `quantCols`
gives the TMT tag and `runCol` the plex, as required by
`QFeatures::readQFeatures(runCol = 'Plex')`. The pooled bridge channel
present in both plexes has a `Genotype` of `Bridge`.

## Usage

``` r
tmt_2plex_design
```

## Format

An object of class `data.frame` with 29 rows and 7 columns.
