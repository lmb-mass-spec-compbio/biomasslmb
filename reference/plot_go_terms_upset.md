# Plot selected GO terms

This function plots a set of GO terms of interest. be run after
[`get_enriched_go`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_enriched_go.md)
and `estimate_go_overrep`. To avoid plotting too many terms, you may
wish to use
[`remove_redundant_go`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/remove_redundant_go.md)
first too. A warning is shown if you try and plot more than 100 GO
terms.

## Usage

``` r
plot_go_terms_upset(
  goi,
  foi,
  gene2cat,
  term_col = "TERM",
  gene_col = "UNIPROTKB"
)
```

## Arguments

- goi:

  `character` GO terms of interest

- foi:

  `character` features (genes/proteins etc) of interest

- gene2cat:

  `data.frame` with 2 columns containing the mapping between genes
  (usually UniProt accessions, Ensembl gene IDs or similar) and GO
  terms.

- term_col:

  `character` column name for GO term description.

- gene_col:

  `character` column name for genes (or protein etc).

## Value

Returns a `ggplot` object.
