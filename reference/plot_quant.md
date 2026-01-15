# Plot distributions for feature intensities per sample.

Given a `SummarizedExperiment`, return a plot of the feature
quantifications per sample.

## Usage

``` r
plot_quant(
  obj,
  method = c("box", "density", "histogram"),
  log2transform = FALSE,
  facet_by_sample = FALSE
)
```

## Arguments

- obj:

  `SummarizedExperiment`.

- method:

  `string`. Plot style. Choice of box, density or histogram plot.

- log2transform:

  `logical`. Should feature quantifications be log-transformed?

- facet_by_sample:

  `logical`. Facet the plot by sample.

## Value

`ggplot` object.
