# Get all ancestor GO terms

Given a data.frame with with a column containing GO terms, this function
will output a data.frame with with the ontologies of those GO terms,
their description, and all their ancestor terms.

## Usage

``` r
get_ancestor_go(
  go_df,
  feature_col = "UNIPROTKB",
  go_col = "GO.ID",
  return_early_debug = FALSE,
  verbose = TRUE
)
```

## Arguments

- go_df:

  `data.frame`. Contains all initial GO terms for proteins of interest
  with `columns == [feature_col, go_col]`.

- feature_col:

  `string`. Name of the column with the features, e.g. `"UNIPROTKB"`.

- go_col:

  `string`. Name of the column with the GO ids, e.g. `"GO.ID"`.

- return_early_debug:

  `logical`. Stop function early and return a list of GO terms and
  ancestor GO terms for debugging purposes.

- verbose:

  `logical`.

## Value

Returns a `data.frame` containing all ancestor GO terms for all proteins
plus GO term descriptions and ontologies.
