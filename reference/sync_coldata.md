# Copy the experimental design onto one or more assays

[`QFeatures::readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html)
attaches `colData` at the `QFeatures` level only, so the individual
`SummarizedExperiment` assays inside it start with an empty `colData`.
Functions which plot or model a single assay —
[`plot_quant()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_quant.md),
[`plot_pca()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_pca.md),
[`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md)
and anything using `limma` — need the design on the assay itself.

`sync_coldata()` copies the object-level `colData` onto the named
assays, matching on column name so that an assay holding a subset of the
samples gets the rows belonging to it.

Assays created by
[`QFeatures::aggregateFeatures()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
and
[`QFeatures::joinAssays()`](https://rformassspectrometry.github.io/QFeatures/reference/joinAssays.html)
also start without `colData`, so this is worth calling after those as
well as after reading the data in.

## Usage

``` r
sync_coldata(obj, i = names(obj))
```

## Arguments

- obj:

  `QFeatures`. Proteomics dataset

- i:

  `character` or `numeric`. Assays to copy the `colData` onto. Defaults
  to every assay in `obj`.

## Value

`QFeatures` with `colData` attached to the named assays

## Examples

``` r
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_clock,
  colData = tmt_clock_design,
  quantCols = rownames(tmt_clock_design),
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

# the assay starts without the design attached
dim(SummarizedExperiment::colData(tmt_qf[["psms_raw"]]))
#> [1] 12  0

tmt_qf <- sync_coldata(tmt_qf)

dim(SummarizedExperiment::colData(tmt_qf[["psms_raw"]]))
#> [1] 12  4
```
