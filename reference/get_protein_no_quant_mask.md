# Identify proteins which have too few features to quantify protein abundance in each sample

Protein level abundances are more accurately quantified where there are
too more features (PSMs/peptides) to summarise from.

Usually, we are performing the summarisation from a matrix
(columns=samples, rows=features) with an associated feature to protein
ID mapping. Within the matrix, some values may be missing (NA). In order
to correctly identify which proteins can be quantified, we need to start
from the feature level object and create a mask which we can use to
replace protein-level quantification values with NA. This is what this
function does. This can then be combined with the `mask_protein_quant`
function to replace protein level quantification values with NA where
they were derived from too few quantification values

## Usage

``` r
get_protein_no_quant_mask(
  obj,
  min_features,
  master_protein_col = "Master.Protein.Accessions",
  plot = FALSE
)
```

## Arguments

- obj:

  `SummarizedExperiment` with PSM or peptide-level quantification

- min_features:

  `numeric` Threshold for minimum features per protein

- master_protein_col:

  `character` Column name for master protein

- plot:

  Set TRUE to plot how many proteins are quantified in each sample.
  Horizontal line represents total number of proteins quantified across
  all samples

## Value

`Matrix` defining whether the protein was quantified from sufficient
features
