# GO term enrichment using goseq

A wrapper function around
[`goseq`](https://rdrr.io/pkg/goseq/man/goseq.html) to perform GO term
enrichment analysis. See the
[`goseq`](https://rdrr.io/pkg/goseq/man/goseq.html) documentation for
details. `pwf` can be made using
[`nullp`](https://rdrr.io/pkg/goseq/man/nullp.html).

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
  [`nullp`](https://rdrr.io/pkg/goseq/man/nullp.html).

- gene2cat:

  `data.frame` with 2 columns containing the mapping between genes
  (usually UniProt accessions, Ensembl gene IDs or similar) and GO
  terms. Alternatively, a `named list` where the names are genes and
  each entry is a `character vector` of GO terms.

- ...:

  Other arguments to be passed to
  [`goseq`](https://rdrr.io/pkg/goseq/man/goseq.html).

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
