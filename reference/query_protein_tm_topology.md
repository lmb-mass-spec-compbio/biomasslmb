# Query UniProt to determine protein transmembrane-domains and topology

Given a set of UniProt IDs, this function returns the information
regarding transmembrane domains and topology.

## Usage

``` r
query_protein_tm_topology(uniprotIDs, verbosity = 0)
```

## Arguments

- uniprotIDs:

  `character vector` Uniprot IDs

- verbosity:

  `integer` Verbosity level for uniprotREST::uniprot_map

## Value

`data.frame` object.

## Examples

``` r
if (FALSE) { # \dontrun{
uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
query_protein_tm_topology(uniprotIDs)
} # }
```
