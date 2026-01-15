# Filter PTM data from MaxQuant to retain those with a given PTM

This function takes a `Qfeatures` object with the evidence.txt output
from MaxQuant and filters it to retain peptides with a given PTM.

## Usage

``` r
filter_maxquant_ptm(obj, i, ptm_prob_col = "Phospho..STY..Probabilities")
```

## Arguments

- obj:

  `QFeatures`. Proteomics dataset

- i:

  `string`. Index for the SummarizedExperiment with data without
  imputation

- ptm_prob_col:

  `character` column name for PTMs.

## Value

Returns a `SummarizedExperiment` where the value for the PTM probability
column is not empty
