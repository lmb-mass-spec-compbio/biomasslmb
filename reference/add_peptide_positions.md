# Add a column describing the position of the peptide sequence with respect to the protein

Identify the position of the peptide sequence in the protein.

The peptide position is undefined (NA) if:

1.  The peptide sequence has multiple master proteins

2.  The protein does not exist in the proteome fasta

3.  The peptide is repeated in the protein sequence.

The peptide_position_info column details why the peptide positions
columns are NA

## Usage

``` r
add_peptide_positions(
  obj,
  proteome_fasta,
  master_protein_col = "Master.Protein.Accessions",
  protein_col_split = "; ",
  sequence_col = "Sequence"
)
```

## Arguments

- obj:

  `SummarizedExperiment` with PD output at PSM/peptide level

- proteome_fasta:

  `character` Filepath for proteome fasta

- master_protein_col:

  `character` Column name for master protein

- protein_col_split:

  `character` Delimiter for multiple proteins in master_protein_col

- sequence_col:

  `character` Column name for peptide sequence

## Value

`SummarizedExperiment` with extra row columns detailing the peptide
position
