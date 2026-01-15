# Filter Spectronaut output

This function filters the precursor .txt file from Spectronaut, based on
various criteria:

1.  Remove features without a master protein (PG.ProteinAccessions
    column)

2.  Remove features without a unique master protein

3.  Remove features matching a contaminant protein

4.  Remove features without quantification values

## Usage

``` r
filter_features_sn(
  obj,
  master_protein_col = "PG.ProteinAccessions",
  unique_master = TRUE,
  filter_contaminant = TRUE,
  contaminant_proteins = NULL,
  remove_no_quant = TRUE,
  cont_string = "Cont_"
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing output from Proteome Discoverer. Use
  [`readQFeatures`](https://rdrr.io/pkg/QFeatures/man/readQFeatures.html)
  to read in .txt file

- master_protein_col:

  `string`. Name of column containing master proteins.

- unique_master:

  `logical`. Filter out features without a unique master protein.

- filter_contaminant:

  `logical`. Filter out features which match a contaminant protein.

- contaminant_proteins:

  `character vector`. The protein IDs form the contaminant proteins

- remove_no_quant:

  `logical`. Remove features with no quantification

- cont_string:

  `string`. string to search for contaminants

## Value

Returns a `SummarisedExperiment` with the filtered Proteome Discoverer
output.
