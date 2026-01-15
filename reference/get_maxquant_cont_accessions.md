# Extract the contaminants protein accessions from a MaxQuant contaminants fasta file

Extract the contaminants protein accessions from a MaxQuant contaminants
fasta file

## Usage

``` r
get_maxquant_cont_accessions(cont_fasta_inf = NULL)
```

## Arguments

- cont_fasta_inf:

  `character` Filepath to the contaminants fasta. Defaults to NULL in
  which case, contaminants.fasta.gz from this package will be used

## Value

Returns a `list` of contaminant IDs
