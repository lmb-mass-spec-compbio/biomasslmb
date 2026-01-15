# Plot the most common missing value patterns

The patterns in missing values can be informative with respect to
whether the experiment has worked, or if particular samples are
outliers. This function uses an 'upset' plot to show the top 50 most
common missing value patterns across the samples

## Usage

``` r
plot_missing_upset(obj, i)
```

## Arguments

- obj:

  `QFeatures`. Proteomics dataset

- i:

  `string`. Index for the SummarisedExperiment you wish to plot

## Value

Returns a *ggplot* object.

## Examples

``` r
set.seed(11)
library(ggplot2)

df <- diamonds[sample(nrow(diamonds), 1000), ]

tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
  quantCols = 36:45,
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.





```
