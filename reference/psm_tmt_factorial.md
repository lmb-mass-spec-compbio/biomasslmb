# PSM-level MaxQuant output for a TMT18plex whole-proteome factorial experiment

MaxQuant `evidence.txt` PSM-level output for a real TMT18plex
whole-proteome experiment with a crossed design: three cell lines
(`Line_A`, `Line_B` and the double `Line_AB`) each treated with a
vehicle (`Control`) or a compound (`Treated`), in 3 replicates.
Subsetted to ~1500 randomly-selected proteins plus samples of
contaminant, decoy and non-unique-master-protein PSMs, and used as the
MaxQuant worked example in the
`TMT workflow: PSM QC and protein summarisation` vignette and as the
factorial worked example in `Data exploration and statistical testing`.
Sample columns are named to match the rownames of
[`tmt_factorial_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_factorial_design.md).

The proteins were sampled at random rather than chosen for effect size,
so the number of differentially abundant proteins reported in the
vignettes is what a random slice of this experiment gives, not what the
full experiment gives.

The design labels, the `Raw.file` column and PSMs identifying the tagged
constructs themselves have been anonymised or excluded, since those
identify the source project rather than describing the experiment.
Columns cross-referencing MaxQuant output tables that are not packaged
here, and columns constant across every retained PSM, have been dropped.

## Usage

``` r
psm_tmt_factorial
```

## Format

An object of class `data.frame` with 23167 rows and 65 columns.
