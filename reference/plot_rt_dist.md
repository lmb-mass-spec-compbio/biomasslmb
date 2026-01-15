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
