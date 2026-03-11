# Add peptide positions with respect to the protein

This function takes a SummarizedExperiment object and adds rowData
columns to describe the positions of the peptides within the proteins',
taking into account the possible peptides according to the digestion
enzyme used.

## Usage

``` r
add_peptide_positions_from_cleavage(
  obj,
  proteome_fasta,
  digest_enzyme = "trypsin-simple",
  missed_cleavages = c(0, 1, 2),
  master_protein_col = "Leading.razor.protein",
  sequence_col = "Sequence",
  start_col = "start",
  end_col = "end"
)
```

## Arguments

- obj:

  SummarizedExperiment. Proteomics dataset

- proteome_fasta:

  character filepath to fasta with protein sequences

- digest_enzyme:

  character. Enzyme used. See ?cleaver::cleave

- missed_cleavages:

  numeric. Vector of allowed number of missed cleavages

- master_protein_col:

  character. Name of column containing master proteins

- sequence_col:

  character. Name of column containing peptide sequences

- start_col:

  character. Name of output column containing peptide start positions

- end_col:

  character. Name of output column containing peptide end positions

## Value

Returns a SummarizedExperiment with additional rowData columns
describing peptide start and end positions with respect to the protein
