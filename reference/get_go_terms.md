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

- expand_terms:

  `logical` Should GO terms be expanded to include all ancestors

- verbosity:

  `integer` Verbosity level for uniprotREST::uniprot_map

- uniprotIDs:

  `character vector` Uniprot IDs

## Value

`data.frame` object.

## Examples

``` r
uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
get_go_terms(uniprotIDs, expand_terms=TRUE)
#> Error in (function (cond) .Internal(C_tryCatchHelper(addr, 1L, cond)))(structure(list(message = "there is no package called ‘GO.db’",     call = loadNamespace(x), package = "GO.db", lib.loc = NULL), class = c("packageNotFoundError", "error", "condition"))): error in evaluating the argument 'x' in selecting a method for function 'select': there is no package called ‘GO.db’
```
