# Parse the PTM probabilities from Proteome Discoverer and add new columns with PTM information

Extract PTM information from the PTM probabilities column in the output
from PD and add new columns with this information. Also optionally
prints a summary of how many features (PSMs/peptides) pass a given
probability threshold

## Usage

``` r
parse_PTM_scores_pd(
  obj,
  threshold = 95,
  ptm_col = "ptmRS.Best.Site.Probabilities",
  prob_split = "; |: ",
  collapse_delimiter = "; ",
  verbose = TRUE
)
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset

- threshold:

  `numeric` If any score is below a set threshold, disregard all
  putative PTM sites

- ptm_col:

  `character` Columm name for PTM probabilities

- prob_split:

  `character` regex to split PTM probabilities

- collapse_delimiter:

  `character` delimiter for multiple values in output columns

- verbose:

  Set TRUE to print log of events to console

## Value

`SummarizedExperiment`
