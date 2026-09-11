# Filter to remove peptides from proteins failing the protein FDR threshold

This function filters the content of a summarised experiment based on a
separate file detailing the proteins passing the FDR confidence
threshold

This function is used to filter peptides to retain only those from
proteins passing the FDR threshold

Note that peptides whose master protein has no row in
`protein_fdr_filename` are also removed, since no FDR confidence can be
established for them. This includes every peptide whose master protein
is a group of several accessions, which cannot match a single-accession
row in the protein-level output. Applying this function after filtering
to unique master proteins (see
`filter_features_pd_dda(unique_master = TRUE)`) keeps the two effects
separate. Accessions listed in `retain_proteins` are exempt.

## Usage

``` r
filter_by_protein_fdr(
  obj,
  protein_fdr_filename,
  protein_col_peptide = "Master.Protein.Accessions",
  protein_col_protein = "Accession",
  protein_FDR_col = "Protein.FDR.Confidence.Combined",
  retain_proteins = NULL
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing peptide-level output from Proteome
  Discoverer.

- protein_fdr_filename:

  `string`. Filepath for protein-level output from Proteome Discoverer

- protein_col_peptide:

  `string`. Name of column containing master proteins in the
  peptide-level `obj`

- protein_col_protein:

  `string`. Name of column containing master proteins in the
  protein-level `protein_fdr_filename`

- protein_FDR_col:

  `string`. Name of column containing FDR information in the
  protein-level `protein_fdr_filename`

- retain_proteins:

  `character vector`. Vector of protein accessions to always retain,
  even if they do not pass the FDR threshold. Default is NULL.

## Value

Returns a `SummarisedExperiment` with the filtered Proteome Discoverer
output.
