# Get all mappings for GO terms

For a vector of GO terms, obtain all of the ancestor or offspring terms.

## Usage

``` r
get_all_mappings(go_ids, ontologies, verbose = TRUE, direction = "ancestor")
```

## Arguments

- go_ids:

  `character vector`. GO terms to use in the format `GO:1234567`.

- ontologies:

  `named character vector`. Names = `go_ids` and values = ontologies
  e.g. `BP`, `CC`, or `MF`.

- verbose:

  `logical`.

- direction:

  `string` Either `"ancestor"` or `"offspring"`.

## Value

Returns a `named list` of `character vectors`. Names == GO terms and
values == vectors containing all ancestor GO terms for the particular
input GO term.
