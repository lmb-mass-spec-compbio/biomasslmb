# Recover annotations for retired UniProtKB accessions via UniParc

MaxQuant searches are run against a FASTA snapshot of UniProt taken at
search time. By the time the annotation step runs, some accessions have
since been deleted from UniProtKB (most commonly because the underlying
Ensembl gene model was withdrawn or superseded), so a normal
[`uniprotREST::uniprot_map()`](https://csdaw.github.io/uniprotREST/reference/uniprot_map.html)
call (as used by
[`get_go_terms`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_go_terms.md))
returns blank `Gene.Names`/`Protein.names` for them even though the
protein may be a real, well-supported hit.

For each accession, this function queries the live UniProtKB entry. If
it is inactive, it follows the entry's `uniParcId` to UniParc and
returns the most recently-updated UniParc cross-reference that has both
a gene name and a protein name (some cross-references, e.g. Ensembl,
carry a gene name but no protein name). This two-step lookup (rather
than bulk ID mapping to UniParc) is used because accessions are
versioned internally (e.g. `A0A5G2QPJ4.1` vs `.2`), and bulk mapping to
UniParc can return multiple UniParc entries for one accession with no
clean way to disambiguate; resolving via the entry's own `uniParcId` is
unambiguous.

Intended to be called only with the subset of accessions that came back
empty/`"deleted"` from a prior `uniprot_map()` call, not a full
accession list.

The two lookups this requires (UniProtKB entry, then UniParc entry) are
each done for all `accessions` in parallel via
[`httr2::req_perform_parallel()`](https://httr2.r-lib.org/reference/req_perform_parallel.html),
rather than one accession at a time, since this is the dominant runtime
cost of
[`get_uniprot_details`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
for inputs with many retired accessions.

## Usage

``` r
get_uniparc_fallback(accessions)
```

## Arguments

- accessions:

  `character vector` UniProt accessions to resolve via UniParc.

## Value

`data.frame` with columns `UniprotID`, `Gene.Names`, `Gene.Names.First`,
`Protein.names`, `Annotation.Source`. Accessions that are not inactive,
have no `uniParcId`, or have no UniParc cross-reference with both a gene
name and a protein name are omitted from the result.

## Examples

``` r
if (FALSE) { # \dontrun{
get_uniparc_fallback("A0A5G2QPJ4")
} # }
```
