# Adds new feature describing the average reporter Signal/Noise ratio.

PD column `Average.Reporter.SN` is NA when all tags have missing values
and imputes a value of zero for NA when there not all tags are missing.
So, with a 10-plex TMT experiment, if a single tag has a SN of 10 and
all others are NA, PD will report an average SN of 10/10 = 1.

This function adds an average reporter SN column which ignores missing
values. In the example above, the average SN value reported would be 10.
Where all values are NA, the average SN will remain NA.

The assumption is that the intensity values in the `exprs` matrix are
Signal/Noise ratios, as reported by PD by default

## Usage

``` r
update_average_sn(obj, sn_col = "Average.Reporter.SN")
```

## Arguments

- obj:

  `SummarizedExperiment`. Should contain PSMs-level TMT quantification

- sn_col:

  `string`. Name of output column containing the average signal:noise

## Value

Returns an `MSnSet` with the average SN included as a new feature
column.
