# LFQ-DDA TurboID pulldown data

`QFeatures` object holding the TurboID proximity labelling dataset
(`lfq_dda_pd_turboid_PeptideGroups.txt`), three biotin-treated samples
against three untreated controls, processed through peptide QC and
summarisation to protein level with `robustSummary`, as produced in Part
B of the `Enrichment designs: IP, BioID and TurboID` vignette.

The `protein` assay is masked so that every value rests on at least two
peptides, and so retains the presence/absence structure a pulldown
produces. `protein_imputed` and `protein_restricted` are the blanket and
restricted imputations of it, the latter filling only control samples
with at most one quantified replicate.

## Usage

``` r
lfq_qf_turboid
```

## Format

An object of class `QFeatures` of length 7.
