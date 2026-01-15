# Add a column with amino acid sequence around a PTM

Add a rowData column with amino acid sequence around a PTM. Value will
be NA if peptide maps to multiple proteins or has multiple PTMs. If
padding extends outside the protein AA sequence, padding will be
extended with '\_'. The PTM- centered AA sequence is useful to integrate
external databases. Input should have been been passed through
`add_ptm_positions` to add the 'ptm_positions_prot' column to the
rowData.

## Usage

``` r
add_site_sequence(
  obj,
  proteome_fasta,
  master_protein_col = "Master.Protein.Accessions",
  sequence_pad = 7
)
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset

- proteome_fasta:

  `character` Filepath for proteome fasta

- master_protein_col:

  `character` Column name for master protein

- sequence_pad:

  `numeric` Number of amino acids to include on either side of PTM

## Value

`SummarizedExperiment`
