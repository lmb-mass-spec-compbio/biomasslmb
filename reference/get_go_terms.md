# Obtain GO term annotations for proteins

Given a set of UniProt IDs, this function queries UniProt to obtain the
annotated GO terms. Optionally, these GO terms can be expanded to
include all ancestors too, which can be helpful when using GO
over-representation/enrichment tools that do not consider the GO term
heirarchy. Note that expansion will significantly increase the run-time

## Usage

``` r
get_go_terms(UniprotID, expand_terms = FALSE, verbosity = 0)
```

## Arguments

- UniprotID:

  `character vector` Uniprot IDs

- expand_terms:

  `logical` Should GO terms be expanded to include all ancestors

- verbosity:

  `integer` Verbosity level for uniprotREST::uniprot_map

## Value

`data.frame` object.

## Examples

``` r
if (FALSE) { # \dontrun{
uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
get_go_terms(uniprotIDs, expand_terms=TRUE)
} # }
```
