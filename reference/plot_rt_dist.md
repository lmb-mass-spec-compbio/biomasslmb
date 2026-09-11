# Plot the retention time distribution

It's useful to assess the retention time distribution as a quality
control step. The peptides should be spread across the retention times,
without any clear 'gaps' or very clear trend towards one end of the
gradient.

## Usage

``` r
plot_rt_dist(obj, rt_col = "RT.in.min.by.Search.Engine.Sequest.HT")
```

## Arguments

- obj:

  `SummarisedExperiment` containing peptide-level output from Proteome
  Discoverer.

- rt_col:

  `string`. Name of column with retention time

## Value

Returns a `ggplot` with the RT vs delta PPM

## Examples

``` r
lfq_inf <- system.file("extdata", "lfq_dda_pd_PeptideGroups.txt",
                       package = "biomasslmb")

lfq_qf_raw <- QFeatures::readQFeatures(
  assayData = read.delim(lfq_inf),
  quantCols = grep("^Abundance", colnames(read.delim(lfq_inf))),
  name = "peptides")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

plot_rt_dist(lfq_qf_raw[["peptides"]],
             rt_col = "PSM.RT.in.min.by.Search.Engine.CHIMERYS")
```
