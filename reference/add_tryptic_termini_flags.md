# Evaluate peptide termini against theoretical cleavage positions

Evaluate peptide termini against theoretical cleavage positions

## Usage

``` r
add_tryptic_termini_flags(
  obj,
  proteome_fasta,
  digest_enzyme = "trypsin-simple",
  missed_cleavages = c(0, 1, 2),
  master_protein_col = "Leading.razor.protein",
  start_col = "start",
  end_col = "end"
)
```

## Arguments

- obj:

  SummarizedExperiment

- proteome_fasta:

  character. Path to FASTA file

- digest_enzyme:

  character. Enzyme for cleaver::cleavageRanges

- missed_cleavages:

  numeric vector. Allowed missed cleavages

- master_protein_col:

  character. Column in rowData containing protein IDs

- start_col:

  character. Column containing peptide start positions

- end_col:

  character. Column containing peptide end positions

## Value

SummarizedExperiment with additional rowData columns: N_tryp and C_tryp
