# Extract the number of samples each feature was detected in for each experiment in a Qfeatures object

When exploring the proteome coverage, it is instructive to consider how
many samples a protein was present in for each level of filtering. This
function takes a Qfeatures object and tallies how many samples a
variable in the `rowData` was present in. When the row variables are the
protein names, the tally is thus the number of samples the protein was
present in.

## Usage

``` r
get_samples_present(obj, rowVars, rename_cols = NULL)
```

## Arguments

- obj:

  `QFeatures` object

- rowVars:

  `character vector` row variables to group by, e.g
  'Master.Protein.Accessions'

- rename_cols:

  `named character vector` optional named list to rename the
  experiments. List values should be current experiment names and list
  names should be updated experiment names

## Value

`data.frame` object.

## Examples

``` r
rename_cols <- c("All PSMs" = "psms_raw", "Protein" = "protein")

samples_present <- get_samples_present(
  tmt_qf, rowVars = "Master.Protein.Accessions", rename_cols = rename_cols)

head(samples_present)
#> # A tibble: 6 × 3
#>   Master.Protein.Accessions `All PSMs` Protein
#>   <chr>                          <int>   <int>
#> 1 A0A286YCX6                        12      12
#> 2 A0A5F8MQ13                        12      12
#> 3 A2A8Z1                            12      12
#> 4 A2AHC3                            12      12
#> 5 A2AHL1                            12      12
#> 6 A2AT37                            12      12
```
