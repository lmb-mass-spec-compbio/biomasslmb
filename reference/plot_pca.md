# Create a Principal Component plot from the feature quantification

A PCA visualisation of feature quantifications in each sample can allow
one to see how the experimental conditions relate to the sources of
variance (principal components). This function plots a PCA, with the
option to colour and/or shape the points by experimental conditions. The
percentage values indicated on the axes are the variance explained by
the PCs.

## Usage

``` r
plot_pca(
  obj,
  i,
  allowing_missing = FALSE,
  colour_by = NULL,
  shape_by = NULL,
  x = 1,
  y = 2,
  ...
)
```

## Arguments

- obj:

  `QFeatures`. Proteomics dataset

- i:

  `string`. Index for the SummarizedExperiment you wish to plot

- allowing_missing:

  `logical`. If TRUE, will use pcaMethods::pca to allow for missing
  values. If FALSE (default), will use stats::prcomp and remove any
  features with missing values

- colour_by:

  `string`. ColData column to colour points by

- shape_by:

  `string`. ColData column to shape points by

- x:

  `numeric`. Principal component to plot on x-axis

- y:

  `numeric`. Principal component to plot on y-axis

- ...:

  additional arguments passed onto
  [`pcaMethods::pca`](https://rdrr.io/pkg/pcaMethods/man/pca.html) (when
  `allowing_missing=TRUE`) or
  [`stats::prcomp`](https://rdrr.io/r/stats/prcomp.html) (when
  `allowing_missing=FALSE`)

## Value

Returns a *ggplot* object.

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

plot_pca(tmt_qf, "psms_raw", colour_by = "Condition")
```
