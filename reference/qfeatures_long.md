# Compatibility wrapper for QFeatures long-format extraction

Dispatches to the appropriate QFeatures long-format function depending
on the Bioconductor version and object class.

## Usage

``` r
qfeatures_long(object, ...)
```

## Arguments

- object:

  A QFeatures object (or object class supported by longForm/longFormat)

- ...:

  Additional arguments passed to the underlying long-format function,
  typically including `colvars` and `rowvars`.

## Value

A data structure identical to the output of the underlying long-format
function (usually a `LongTable` or `DataFrame`).

## Details

Specifically, this function:

- Uses
  [`longForm()`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html)
  (BiocGenerics generic) if a method exists for the input object,

- Falls back to `QFeatures::longFormat()` otherwise.

This allows package or script code to remain compatible across QFeatures
releases without explicitly checking Bioconductor versions.

The function first checks whether
[`longForm()`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html) is a
generic and whether there is a registered method for the input object
class. If so, it calls
[`longForm()`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html).
Otherwise, it falls back to the legacy `longFormat()` function. An error
is thrown if neither is available.

## See also

`longFormat`,
[`longForm`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html),
[`isGeneric`](https://rdrr.io/r/methods/GenericFunctions.html),
[`findMethods`](https://rdrr.io/r/methods/findMethods.html)

## Examples

``` r
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
  quantCols = 36:45,
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

# long format, one row per feature per sample, for plotting with ggplot2
head(qfeatures_long(tmt_qf[["psms_raw"]]))
#> DataFrame with 6 rows and 4 columns
#>       rowname       colname     value assayName
#>   <character>   <character> <numeric> <integer>
#> 1           1 Abundance....     465.7         1
#> 2           2 Abundance....      96.8         1
#> 3           3 Abundance....     196.8         1
#> 4           4 Abundance....     609.2         1
#> 5           5 Abundance....     361.7         1
#> 6           6 Abundance....      52.6         1
```
