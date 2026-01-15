# Replace protein-level quantifications with NA if they derive from too few lower feature level quantifications

Protein level abundances are more accurately quantified where there are
more features (PSMs/peptides) to summarise from.

Usually, we are performing the summarisation from a matrix
(columns=samples, rows=features) with an associated feature to protein
ID mapping. Within the matrix, some values may be missing (NA). In order
to correctly identify which proteins can be quantified, we need to start
from the feature level object and create a mask which we can use to
replace protein-level quantification values with NA. This is can be
obtained with `get_protein_no_quant_mask `. This can then be combined
with this function to replace protein level quantification values with
NA where they were derived from too few quantification values

## Usage

``` r
mask_protein_level_quant(obj, retain_mask)
```

## Arguments

- obj:

  `SummarizedExperiment` with PSM or peptide-level quantification

- retain_mask:

  `matrix` detailing whether a protein-level quantification had
  sufficient lower level quantification values for each sample. Can be
  obtained with `get_protein_no_quant_mask `

## Value

`SummarizedExperiment` with quantification values replaced by NA where
they derive from too few lower feature level quantifications
