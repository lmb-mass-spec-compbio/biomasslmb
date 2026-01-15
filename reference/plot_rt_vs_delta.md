# Plot the retention time vs delta precursor mass

It's useful to assess the relationship between retention time and delta
precursor mass as a quality control step. There should be little or no
relationship

## Usage

``` r
plot_rt_vs_delta(
  obj,
  rt_col = "RT.in.min.by.Search.Engine.Sequest.HT",
  delta_ppm_col = "Delta.M.in.ppm.by.Search.Engine.Sequest.HT"
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing peptide-level output from Proteome
  Discoverer.

- rt_col:

  `string`. Name of column with retention time

- delta_ppm_col:

  `string`. Name of column with the Delta precursor mass

## Value

Returns a `ggplot` with the RT vs delta PPM
