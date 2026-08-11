# Get UniProt annotation details, with UniParc fallback for retired accessions

Queries
[`uniprotREST::uniprot_map()`](https://csdaw.github.io/uniprotREST/reference/uniprot_map.html)
for gene/protein name annotations for a vector of UniProt accessions,
following the same `uniprot_map()` usage as
[`get_go_terms`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_go_terms.md).
Because MaxQuant searches are run against a UniProt snapshot that may be
weeks to years out of date by the time annotation runs, some accessions
may since have been deleted from UniProtKB; these come back from
`uniprot_map()` with blank `Gene.Names` and
`Protein.names == "deleted"`. For any such accession,
[`get_uniparc_fallback`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniparc_fallback.md)
is used to recover the last-known annotation from UniParc.

## Usage

``` r
get_uniprot_details(accessions, verbosity = 0)
```

## Arguments

- accessions:

  `character vector` UniProt accessions.

- verbosity:

  `integer` Verbosity level for uniprotREST::uniprot_map. Some
  accessions come back from `uniprot_map()` as multiple rows (they have
  subsequently been 'demerged'); these are collapsed to one row per
  accession, with `;`-separated values where the demerged rows differed.

## Value

`data.frame` with one row per accession and columns `UniprotID`,
`Entry.Name`, `Reviewed`, `Protein.names`, `Gene.Names`, `Organism`,
`Length`, `Gene.Names.First`, `Annotation.Source`. `Annotation.Source`
is `"UniProtKB (live)"` for accessions resolved directly from UniProtKB,
or `"UniParc (retired UniProtKB entry)"` for accessions recovered via
[`get_uniparc_fallback`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniparc_fallback.md).
Accessions that are deleted and have no recoverable UniParc annotation
keep their original (blank/`"deleted"`) values from `uniprot_map()`.

## Examples

``` r
if (FALSE) { # \dontrun{
get_uniprot_details(c('O76024', 'Q03135', 'A0A5G2QPJ4'))
} # }
```
