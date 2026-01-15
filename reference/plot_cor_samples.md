# Plot the correlation between sample

The correlation between feature quantifications in each sample can be a
useful QC step to assess whether the experiment has worked. This
function plots the correlation values in a heatmap. Samples are not
clustered and are ordered on the axes in the order they are in the
SummarizedExperiment

## Usage

``` r
plot_cor_samples(obj, i, order = "original", ...)
```

## Arguments

- obj:

  `QFeatures`. Proteomics dataset

- i:

  `string`. Index for the SummarizedExperiment you wish to plot

- ...:

  addiional arguments passed onto corrplot::corrplot

## Value

Returns a *ggplot* object.

## Examples

``` r
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
  quantCols = 36:45,
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

cor_sample(tmt_qf, 'psms_raw')
#> Error in cor_sample(tmt_qf, "psms_raw"): could not find function "cor_sample"
```
