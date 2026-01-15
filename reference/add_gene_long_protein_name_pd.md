# Extract the gene name and long-form protein name from the master protein descriptions column

This function extracts the gene name and long-form protein name from the
master protein descriptions column in output from Proteome Discoverer
for DDA.

It assumes the master protein description column is in the format (.*)
.* GN=(.*) PE.*, where the first group is the long format of the protein
name and the second group is the gene name. Where the gene name is not
included, an empty string is returned An example of the expected input
in the column is: Serum albumin OS=Bos taurus GN=ALB PE=1 SV=4

## Usage

``` r
add_gene_long_protein_name_pd(
  obj,
  master_protein_desc_col = "Master.Protein.Descriptions",
  gene_name_out_col = "Master.Protein.gene.name",
  long_protein_name_out_col = "Master.Protein.long.name"
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing output from Proteome Discoverer

- master_protein_desc_col:

  `string`. Name of column containing master proteins descriptions.

- gene_name_out_col:

  `string`. Name of output column containing the gene names

- long_protein_name_out_col:

  `string`. Name of output column containing the long format protein
  names

## Value

Returns a `QFeatures` with the filtered Proteome Discoverer output.
