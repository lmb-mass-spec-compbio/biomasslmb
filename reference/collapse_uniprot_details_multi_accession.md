# Map per-accession UniProt details onto protein IDs with multiple accessions

Takes the `data.frame` returned by
[`get_uniprot_details`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
(one row per single accession) and a vector of protein IDs where an
element may itself be several accessions joined by `sep` – the
`Master.Protein.Accessions` column from Proteome Discoverer
`PeptideGroups.txt` output is the common source of this. Each
`protein_ids` element is split into its constituent accessions, the
matching rows of `uniprot2details` are looked up, then re-collapsed back
to one row per original (possibly multi-accession) `protein_ids`
element. Every constituent accession contributes one value to each
column, joined by `sep`, so a group of `n` accessions yields `n`-part
`sep`-separated strings that stay positionally aligned across columns.
Does not call
[`get_uniprot_details`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
itself, so the same `uniprot2details` result can also be used on its
own, without fetching annotations twice.

## Usage

``` r
collapse_uniprot_details_multi_accession(
  uniprot2details,
  protein_ids,
  sep = "; "
)
```

## Arguments

- uniprot2details:

  `data.frame`, the output of
  [`get_uniprot_details`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md),
  covering (at least) every accession that appears in `protein_ids` once
  split by `sep`.

- protein_ids:

  `character vector` protein IDs, each optionally containing several
  UniProt accessions joined by `sep`.

- sep:

  `string` delimiter separating accessions within a `protein_ids`
  element, and also used to join differing values in the collapsed
  output. Default `"; "`, matching Proteome Discoverer's
  `Master.Protein.Accessions` column.

## Value

`data.frame` with one row per unique `protein_ids` element and the same
columns as `uniprot2details`. For a `protein_ids` element made of `n`
accessions, each column holds the `n` per-accession values joined by
`sep`, in the same order across every column.

## Examples

``` r
if (FALSE) { # \dontrun{
uniprot2details <- get_uniprot_details(c('O76024', 'Q03135', 'A0A5G2QPJ4'))
collapse_uniprot_details_multi_accession(
  uniprot2details, c('O76024', 'Q03135; A0A5G2QPJ4'))
} # }
```
