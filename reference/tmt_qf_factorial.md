# Protein-level abundances for the TMT factorial experiment

`QFeatures` object containing the protein-level abundances produced from
[`psm_tmt_factorial`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_factorial.md)
by the PSM filtering and summarisation in the
`TMT workflow: PSM QC and protein summarisation` vignette, read back in
by `Data exploration and statistical testing`. Only the `protein` assay
is retained: the PSM-level assays it was built from are available as
[`psm_tmt_factorial`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_factorial.md),
and keeping them here would make the object an order of magnitude
larger.

## Usage

``` r
tmt_qf_factorial
```

## Format

An object of class `QFeatures` of length 1.
