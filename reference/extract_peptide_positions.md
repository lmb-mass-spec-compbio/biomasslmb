# Extract Peptide Position Start and End Columns from a QFeatures Object

Parses the `Positions.in.Master.Proteins` column from the `rowData` of a
specified `SummarizedExperiment` within a `QFeatures` object, and
extracts the start and end positions of each peptide into two new
`rowData` columns. When multiple entries are present
(semicolon-delimited), the function returns semicolon-delimited lists of
start and end values.

## Usage

``` r
extract_peptide_positions(qf, i, start_col = "start", end_col = "end")
```

## Arguments

- qf:

  A `QFeatures` object.

- i:

  A single integer or character string specifying the index or name of
  the `SummarizedExperiment` within `qf` to operate on.

- start_col:

  A single character string giving the name of the new column to store
  peptide start positions. Defaults to `"start"`.

- end_col:

  A single character string giving the name of the new column to store
  peptide end positions. Defaults to `"end"`.

## Value

The input `QFeatures` object with two additional columns in the
`rowData` of the specified experiment:

- start_col:

  Character. Semicolon-delimited start position(s) for the peptide.

- end_col:

  Character. Semicolon-delimited end position(s) for the peptide.

## Examples

``` r
if (FALSE) { # \dontrun{
library(QFeatures)
# Using default column names
qf <- extract_peptide_positions(qf, i = "peptides")

# Using custom column names
qf <- extract_peptide_positions(qf, i = 1L,
                                start_col = "pep_start",
                                end_col   = "pep_end")
} # }
```
