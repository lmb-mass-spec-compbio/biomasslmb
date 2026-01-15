# Plot the missing values vs signal:noise for each sample

Missing values are more frequent with low signal:noise (S:N). This
function visualises this relationship for each sample to aid selection
of thresholds for minimal S:N filtering.

## Usage

``` r
plot_missing_SN_per_sample(obj, sn_column = "Average.Reporter.SN", bins = 20)
```

## Arguments

- obj:

  `SummarizedExperiment` containing PSM level TMT intensities

- sn_column:

  `character` column name for Signal:noise values

- bins:

  `numeric` Number of bins to plot

## Value

`ggplot` tile plot to show S:N vs \# missing values for each sample
