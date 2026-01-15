# Plot quantification at multiple levels

For the exploration of specific proteins, it's often useful to visualise
the quantification at multiple levels, e.g PSM, peptide, filtered
peptides and protein. This function plots the quantitative data at
multiple levels (assays) for a given list of proteins of interest.

## Usage

``` r
plot_protein_assays(
  obj,
  poi,
  experiments_to_plot = NULL,
  protein_id_col = "Master.Protein.Accessions",
  label_col = "Master.Protein.Accessions",
  row_vars = NULL,
  rename_labels = NULL,
  log2transform_cols = "",
  norm_quant = FALSE,
  add_mean_summary = FALSE,
  colour_assays = FALSE,
  alpha_assays = TRUE
)
```

## Arguments

- obj:

  `QFeatures` object

- poi:

  `list`. Proteins of interest

- experiments_to_plot:

  `list`. experiments (assays) to plot. Defaults to all assays

- protein_id_col:

  `string`. Column with the protein ids to search for `poi` values

- label_col:

  `string`. Column with labels to use for proteins in plot

- row_vars:

  `character vector`. Additional row variables to include in the plot
  data (e.g peptide sequence)

- rename_labels:

  `named list`. Mapping from labels to renamed labels

- log2transform_cols:

  `list`. Assays which need to be log2-transforms (all values should
  ultimately be transformed)

- norm_quant:

  `logical`. Should the quantifications be normalised to the fold-change
  vs mean abundance

- add_mean_summary:

  `logical`. Add a line summarising the mean value over all features in
  the assay (Not recommended if norm_quant=FALSE)

- colour_assays:

  `logical`. Column each assay

- alpha_assays:

  `logical`. Set sensible alpha value for each assay

## Value

`ggplot` object.
