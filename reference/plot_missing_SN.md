# Plot the missing values vs signal:noise

Missing values are more frequent with low signal:noise (S:N). This
function visualises this relationship to aid selection of thresholds for
minimal S:N filtering

## Usage

``` r
plot_missing_SN(obj, sn_column = "Average.Reporter.SN", bins = 20)
```

## Arguments

- obj:

  `SummarizedExperiment` containing PSM level TMT intensities

- sn_column:

  `character` column name for signal:noise values

- bins:

  `numeric` Number of bins to plot

## Value

`ggplot` stacked bar plot to show S:N vs \# missing values
