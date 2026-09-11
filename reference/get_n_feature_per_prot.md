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

## Examples

``` r
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
  quantCols = 36:45,
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

head(get_n_feature_per_prot(tmt_qf[["psms_raw"]]))
#>   Master.Protein.Accessions         sample n
#> 1                            Abundance.126 2
#> 2                           Abundance.127C 2
#> 3                           Abundance.127N 2
#> 4                           Abundance.128C 2
#> 5                           Abundance.128N 2
#> 6                           Abundance.129C 2
```
