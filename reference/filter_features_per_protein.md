# Remove features which are assigned to a protein with too few supporting features in total

For summarisation of PSM or peptide to protein, we need a minimum number
of finite values per protein per sample.

## Usage

``` r
filter_features_per_protein(
  obj,
  min_features,
  master_protein_col = "Master.Protein.Accessions",
  plot = FALSE
)
```

## Arguments

- obj:

  `SummarizedExperiment` with PSM or peptide-level quantification

- min_features:

  `numeric` Threshold for minimum features per protein

- master_protein_col:

  `character` Column name for master protein

- plot:

  Set TRUE to plot histogram of features per protein per sample

## Value

`SummarizedExperiment`
