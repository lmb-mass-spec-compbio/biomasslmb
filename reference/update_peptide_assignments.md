# Update PSM level protein/master protein assignments using peptide-level assignments

When multiple search engines are used in Proteome Discoverer for TMT
proteomics the PSM to protein and master protein assignments may not be
consistent between the search engines. We can resolve this using the
peptide level output from Proteome Discoverer, which has consistent
assignments

## Usage

``` r
update_peptide_assignments(
  obj,
  pep_inf,
  verbose = FALSE,
  master_protein_col = "Master.Protein.Accessions"
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing PSM-level output from Proteome
  Discoverer.

- pep_inf:

  `string`. Filepath for Peptide-level output

- verbose:

  `boolean`. Default is FALSE; Don't output tallies of PSMs/proteins.

- master_protein_col:

  `string`. Name of column containing master proteins. (Only used when
  verbose=TRUE)

## Value

Returns a `SummarisedExperiment` with the standardised assignments

## Details

**Leucine/Isoleucine** Leucine/Isoleucine are isobaric and only generate
distinct fragments with higher collision energies which are not
typically used in proteomics analysis. Thus, L and I are usually
indistinguishable.

The spectrum search engines deal with this differently in terms of how
they report the multiple possible sequence matches and their
assignments.

The peptide level Proteome Discoverer output selects one possible
sequence and uses the union of reported assignments from the PSM level
for the peptide level assignments. Thus, when we use the peptide level
assignments, we will lose PSM level data for the redundant sequences
which were not selected. This is not an issue but will result in the
function output being shorter than the input.
