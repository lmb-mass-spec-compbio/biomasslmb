# Determine GO ancestor object

For a single GO term, get its ontology and then return the correct
ancestor mapping object for that ontology.

## Usage

``` r
determine_ancestor_function(term, ontology)
```

## Arguments

- term:

  `string`. A single GO.ID.

- ontology:

  `string`. A single ontology, one of: BP, MF, CC.

## Value

Returns an `AnnDbBimap` object which maps GO terms (BP, MF, or CC) to
all ancestor terms. Each GO term is mapped to a vector of ancestor GO
terms.
