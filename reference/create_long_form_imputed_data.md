# Create long format data with column defining if quantification value is imputed

This function takes a `Qfeatures` object with unimputed data and imputed
data in separate assays and creates a long format data.frame with a
column called `is.imputed` which describes whether the quantification
value was imputed

## Usage

``` r
create_long_form_imputed_data(
  obj,
  i_unimputed,
  i_imputed,
  column_variables,
  row_variables
)
```

## Arguments

- obj:

  `QFeatures`. Proteomics dataset

- i_unimputed:

  `string`. Index for the SummarizedExperiment with data without
  imputation

- i_imputed:

  `string`. Index for the SummarizedExperiment with data with imputation

- column_variables:

  `character vector`. Variables in `colData(obj)` to include

- row_variables:

  `character vector`. Variables in `rowData(obj)` to include

## Value

Returns a `data.frame` with quantification data in long form and
including a `is.imputed` column
