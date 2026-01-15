# GO term enrichment using goseq

A wrapper function around `goseq` to perform GO term enrichment
analysis. See the `goseq` documentation for details. `pwf` can be made
using `nullp`.

Over/underrepresented p-values are automatically adjusted using
`method = "BH"`. If `gene2cat` is not provided then this function will
default to using the Homo sapiens genome `hg19` and will expect Ensembl
gene IDs to have been used to construct the `pwf` input.

## Usage

``` r
get_enriched_go(
  pwf,
  gene2cat = NULL,
  ...,
  shorten_term = TRUE,
  shorten_lims = c(1L, 30L),
  filter_no_DE = TRUE
)
```

## Arguments

- pwf:

  `data.frame` with 3 columns (`DEgenes` = logical, `bias.data` =
  numeric/integer, `pwf` = numeric) and row names (usually UniProt
  accessions, Ensembl gene IDs or similar). Typically constructed using
  `nullp`.

- gene2cat:

  `data.frame` with 2 columns containing the mapping between genes
  (usually UniProt accessions, Ensembl gene IDs or similar) and GO
  terms. Alternatively, a `named list` where the names are genes and
  each entry is a `character vector` of GO terms.

- ...:

  Other arguments to be passed to `goseq`.

- shorten_term:

  `logical`. Should an extra column with a substring of the output GO
  terms be added to the output data.frame? Default is `TRUE`.

- shorten_lims:

  `integer vector` of length 2. The start and stop coordinates of the
  substring.

- filter_no_DE:

  `logical`. Should terms without any features in the foreground be
  removed?

## Value

Returns a `data.frame` of over/underrepresented GO terms.
