# Filter to remove peptides from proteins failing the protein FDR threshold

This function filters the content of a summarised experiment based on a
separate file detailing the proteins passing the FDR confidence
threshold

This function is used to filter peptides to retain only those from
proteins passing the FDR threshold

Note that peptides whose master protein has no row in
`protein_fdr_filename` are also removed, since no FDR confidence can be
established for them. With the default `multi_accession = "remove"`,
this includes every peptide whose master protein is a group of several
accessions, which cannot match a single-accession row in the
protein-level output. Applying this function after filtering to unique
master proteins (see `filter_features_pd_dda(unique_master = TRUE)`)
keeps the two effects separate. Where groups have been kept
(`unique_master = FALSE`), `multi_accession = "all"` or `"any"` assesses
them on the confidence of their member accessions. Accessions listed in
`retain_proteins` are exempt.

Groups retained by `multi_accession = "all"` or `"any"` keep their
compound accession string; they are not split or reassigned. On
aggregation to protein level they become separate compound rows, which
can overlap the rows for their member accessions.

## Usage

``` r
filter_by_protein_fdr(
  obj,
  protein_fdr_filename,
  protein_col_peptide = "Master.Protein.Accessions",
  protein_col_protein = "Accession",
  protein_FDR_col = "Protein.FDR.Confidence.Combined",
  retain_proteins = NULL,
  multi_accession = c("remove", "all", "any")
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

- multi_accession:

  `string`. How to treat features whose master protein is a group of
  several accessions, separated by ";". `"remove"` (default) matches the
  whole string, so a group is removed unless the protein-level output
  contains that same string. `"all"` retains the feature if every
  accession has 'High' confidence, with an accession absent from the
  protein-level output counting as not 'High'. `"any"` retains the
  feature if at least one accession has 'High' confidence. Features with
  a single accession are unaffected.

## Value

Returns a `SummarisedExperiment` with the filtered Proteome Discoverer
output.
