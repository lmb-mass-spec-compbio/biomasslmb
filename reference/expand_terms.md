# Expand data.frame GO terms

Input is a data.frame with all the GO terms annotated to a single
protein in a single column. Output is a new data.frame with all of the
GO terms for that protein (annotated and ancestor).

## Usage

``` r
expand_terms(go_df, go_col, go2Ancestor)
```

## Arguments

- go_df:

  `data.frame` for a single protein with a single column where each row
  is a single GO term.

- go_col:

  `variable`. Name of column from the data.frame that contains the GO
  terms.

- go2Ancestor:

  `named list`. Returned by
  [`get_all_mappings`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_all_mappings.md).
  Names == GO terms and values == vectors containing all ancestor GO
  terms for the particular input GO term.

## Value

Returns a `data.frame` containing all ancestor GO terms for a single
protein.
