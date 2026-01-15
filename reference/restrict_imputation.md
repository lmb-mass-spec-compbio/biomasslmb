# Restrict imputed values to specific conditions

This function takes a `Qfeatures` object with unimputed data and imputed
data in separate assays and creates a new assay where the imputed values
are only used in specified circumstances. Note that this restriction
occurs post imputation, which may not be suitable for imputation methods
where the imputed values are not independently derived

## Usage

``` r
restrict_imputation(
  obj,
  i_unimputed,
  i_imputed,
  i_restricted_imputed,
  use_imputed_df,
  verbose = TRUE
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

- i_restricted_imputed:

  `string`. Index for the output assay with restricted imputation

- use_imputed_df:

  `data.frame` with conditions where imputation should be used. see
  @details for more information

- verbose:

  `logical` Describe how many missing values in input and output assays

## Value

Returns a `QFeatures` with restricted imputation with specified assay
name

## Details

The `data.frame` provided to `use_imputed_df` should specify the
conditions under which imputation should be performed, with respect to
the experimental conditions and the number of missing values.
`data.frame` must contain a column called `n_finite` and all other
columns must be columns in `colData(obj)`

For example, in an enrichment vs control experiment, to only input in
control samples where zero or one replicate are quantified, this should
be specified thusly

condition n_finite control 0 control 1
