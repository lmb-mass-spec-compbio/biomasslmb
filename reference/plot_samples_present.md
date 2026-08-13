# Plots the number of samples each feature was detected in for each experiment in a Qfeatures object

When exploring the proteome coverage, it is instructive to consider how
many samples a protein was present in for each level of filtering. This
function plots the output of `get_samples_present`

## Usage

``` r
plot_samples_present(samples_present, rowvars, breaks = NULL)
```

## Arguments

- samples_present:

  `data.frame` output from `get_samples_present`

- rowvars:

  `character vector` row variables on which the tallying was grouped by

- breaks:

  `numeric vector` breaks for the 'Samples' fill scale. If NULL
  (default), breaks are determined automatically

## Value

`ggplot` object.
