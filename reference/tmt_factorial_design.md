# Experimental design for `psm_tmt_factorial`

A `data.frame` giving the `Genotype` (`Line_A`, `Line_B` or `Line_AB`),
`Treatment` (`Control` or `Treated`) and `Replicate` for each sample in
[`psm_tmt_factorial`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_factorial.md),
with rownames matching the quantification column names. The two factors
are fully crossed, with 3 replicates in each of the 6 combinations.

## Usage

``` r
tmt_factorial_design
```

## Format

An object of class `data.frame` with 18 rows and 4 columns.
