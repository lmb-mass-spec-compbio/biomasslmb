# Obtain protein transmembrane-domains and topology information

Given a set of UniProt IDs, this function queries UniProt to obtain the
information regarding transmembrane domains (TMDs) and topology and then
parses this information to define additional columns describining the
TMDs and topology

## Usage

``` r
get_protein_tm_topology(uniprotIDs, verbosity = 0)
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
get_protein_tm_topology(uniprotIDs)
} # }
```
