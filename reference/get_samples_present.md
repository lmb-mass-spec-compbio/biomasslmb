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
