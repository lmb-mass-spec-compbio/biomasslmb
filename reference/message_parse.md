# Report how many features and master proteins remain

Prints the number of rows and the number of distinct master proteins in
a feature-level annotation table, followed by a short note saying which
step the count describes. The `filter_features_*` functions call this
after each filter they apply, so the same message format can be used to
report a count at a point where no filter function ran, such as after
[`filterNA()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
or a manual subset.

## Usage

``` r
message_parse(x, column, note)
```

## Arguments

- x:

  `data.frame` or `DataFrame`. Feature-level annotations, normally the
  output of `rowData()` on an assay.

- column:

  `string`. Name of the column in `x` holding the master protein
  accession, e.g. `"Master.Protein.Accessions"`.

- note:

  `string`. Short description of the step being reported.

## Value

Invisibly `NULL`; called for the message it prints.

## Examples

``` r
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
  quantCols = 36:45,
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

message_parse(SummarizedExperiment::rowData(tmt_qf[["psms_raw"]]),
              "Master.Protein.Accessions",
              "Input")
#> 5000 features found from 2364 master proteins => Input
```
