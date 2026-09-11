# PSM-level PD output for total proteome TMT10-plex data

Proteome Discoverer PSM-level output for a total proteome TMT10-plex
experiment, truncated to 5000 PSMs. Used in the examples for the QC and
plotting functions, where a small PSM-level table with reporter ion
quantification and the PD-specific QC columns (signal:noise,
interference) is all that is needed.

The vignettes use `psm_tmt_clock` instead, which carries an accompanying
experimental design in `tmt_clock_design`.

## Usage

``` r
psm_tmt_total
```

## Format

An object of class `data.frame` with 5000 rows and 50 columns.
