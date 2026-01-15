# Read in data from DIA-NN

This function reads in the report.tsv file from DIA-NN, performs
filtering using the Q-Value columns and then joins together the
quantification from the individual samples

## Usage

``` r
readDIANNFilterQJoin(
  diann_report_infile,
  global_protein_q = 0.01,
  run_protein_q = 0.01,
  run_precursor_q = 0.01,
  return_sep_quant = FALSE
)
```

## Arguments

- diann_report_infile:

  `character`. Filepath to DIA-NN report.tsv file

- global_protein_q:

  `numeric`. Threshold for global protein FDR

- run_protein_q:

  `numeric`. Threshold for run-specific protein FDR

- run_precursor_q:

  `numeric`. Threshold for run-specific precursor FDR

- return_sep_quant:

  `logical`. Return the individual sample-level quantification as well
