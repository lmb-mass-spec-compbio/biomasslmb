# Add rowData columns with positions of PTMs with respect to peptide sequence

This function takes a `SummarizedExperiment` object with the
evidence.txt output from MaxQuant and adds rowData columns to describe
the PTMs present in the peptide and their positions within the peptide'
'

## Usage

``` r
add_filter_ptm_pos_rowdata_mq(
  obj,
  ptms_to_retain = c("Phospho (STY)"),
  prob_col = "Phospho..STY..Probabilities",
  min_prob = 0.501,
  filter_pep_by_prob = TRUE,
  verbose = TRUE
)
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset

- ptms_to_retain:

  `character` vector with PTMs to retain. Inspect values in
  rowData(obj)\$Modified.sequence to get correct PTM name

- prob_col:

  `character` name of column containing PTM probabilities

- min_prob:

  `numeric` Minimum acceptable probability for PTM localisation

- filter_pep_by_prob:

  `logical` Filter the output to only return cases where the number of
  sites passing the probability threshold equals the number of PTMs in
  the peptide

- verbose:

  `logical` Describe the number of PTMs detected

- ptm_encoding_pos:

  `character` name vector describing whether modifification comes before
  (-1) or after (1) amino acid in Modified.sequence column values

## Value

Returns a `SummarizedExperiment` with an additional column in the
RowData describing the position of the PTMs with respect to the peptide
sequence
