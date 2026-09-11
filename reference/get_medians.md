# Extract the assay column medians from an MSnSet

Extract the assay column medians from an MSnSet

## Usage

``` r
get_medians(obj)
```

## Arguments

- obj:

  `SummarisedExperiment`

## Value

`vector` of assay column medians in colnames order

## Examples

``` r
# per-sample median abundance, used as a normalisation reference
get_medians(tmt_qf[["protein"]])
#>        M1        C6        C5        M4        C3        C1        M6        C4 
#>  9.875348  9.933488 10.007576  9.901078  9.950782  9.915917  9.964775  9.934982 
#>        M3        C2        M2        M5 
#> 10.002450  9.943597  9.848576  9.852961 
```
