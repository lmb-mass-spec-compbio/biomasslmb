# Remove redundant Peptide Spectrum Matches

When multiple search engines are used in Proteome Discoverer for TMT
proteomics, the PSM-level quantification contains duplicate rows for
each search engine match to the same spectrum. We can use the scan
number columns and quantitative data to identify the duplicate rows and
exclude them.

The removal of duplicate quantification should be performed following
initial filtering to 1) remove rejected PSMs (see filter_TMT_PSMs) and
2) retain only rank 1 matches (Filter on Rank and Search.Engine.Rank
columns)

## Usage

``` r
remove_redundant_psm_quant(
  obj,
  verbose = FALSE,
  master_protein_col = "Master.Protein.Accessions",
  validate = TRUE
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing PSM-level output from Proteome
  Discoverer.

- verbose:

  `boolean`. Don't output tallies of PSMs/proteins; Default=FALSE

- master_protein_col:

  `string`. Name of column containing master

- validate:

  `logical`. Check that sequence to protein assignments are consistent.
  Default=TRUE proteins. (Only used when verbose=TRUE)

## Value

Returns a `SummarisedExperiment` with the redundant quantification
removed
