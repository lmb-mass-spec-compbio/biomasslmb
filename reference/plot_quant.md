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

plot_quant(tmt_qf[["psms_raw"]], log2transform = TRUE, method = "density")
#> Warning: Removed 49066 rows containing non-finite outside the scale range
#> (`stat_density()`).


# boxplot instead, coloured by a colData column
plot_quant(tmt_qf[["psms_raw"]], log2transform = TRUE, method = "box")
#> Warning: Removed 49066 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
```
