# Identify how many features (PSMs/Peptides) are quantified for each protein

For summarisation of PSM or peptide to protein, we need a minimum number
of finite values per protein per sample. This function simply tallies
how many we have.

## Usage

``` r
get_n_feature_per_prot(obj, master_protein_col = "Master.Protein.Accessions")
```

## Arguments

- obj:

  `SummarizedExperiment` with PSM or peptide-level quantification

- master_protein_col:

  `character` Column name for master protein

## Value

`data.frame` detailing how many features are present for each protein in
each sample
