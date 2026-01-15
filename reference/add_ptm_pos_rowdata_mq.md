# Add rowData columns with positions of PTMs with respect to peptide sequence

This function takes a `SummarizedExperiment` object with the
evidence.txt output from MaxQuant and adds rowData columns to describe
the PTMs present in the peptide and their positions within the peptide'
'

## Usage

``` r
add_ptm_pos_rowdata_mq(
  obj,
  ptms_to_retain = c("p", "(ph)", "(ac)p"),
  ptm_encoding_pos = c(p = -1, `(ph)` = 1, `(ac)p` = -1, `(ox)` = 1),
  verbose = TRUE
)
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset

- ptms_to_retain:

  `character` vector with PTMs to retain. Inspect values in
  rowData(obj)\$Modified.sequence to get correct PTM name

- ptm_encoding_pos:

  `character` name vector describing whether modifification comes before
  (-1) or after (1) amino acid in Modified.sequence column values.

## Value

Returns a `SummarizedExperiment` with an additional column in the
RowData describing the position of the PTMs with respect to the peptide
sequence
