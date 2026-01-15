# Create a Principle Component plot from the feature quantification

A PCA visualisation of feature quantifications in each sample can allow
one to see how the experimental conditions relate to the sources of
variance (principle components). This function plots a PCA, with the
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

  `numeric`. Principle component to plot on x-axis

- y:

  `numeric`. Principle component to plot on x-axis

## Value

Returns a *ggplot* object.
