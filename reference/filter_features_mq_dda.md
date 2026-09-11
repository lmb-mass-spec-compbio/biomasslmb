# Filter MaxQuant DDA output

This function filters the output .txt files (peptide groups or PSMs)
from Proteome Discoverer for DDA, based on various criteria:

1.  Remove hits to the decoy database

2.  Remove features without a master protein

3.  Remove features which are not proteotypic (i.e. Unique..Proteins. is
    yes)

4.  Remove features matching a contaminant protein

5.  Remove features matching any protein associated with a contaminant
    protein (see below)

6.  Remove features without quantification values

## Usage

``` r
filter_features_mq_dda(
  obj,
  master_protein_col = "Leading.razor.protein",
  protein_col = "Proteins",
  unique_master = FALSE,
  proteotypic = FALSE,
  filter_contaminant = TRUE,
  contaminant_proteins = NULL,
  filter_associated_contaminant = TRUE,
  remove_no_quant = TRUE
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing output from MaxQuant. Use
  [`readQFeatures`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html)
  to read in .txt file

- master_protein_col:

  `string`. Name of column containing master proteins.

- protein_col:

  `string`. Name of column containing all protein matches.

- unique_master:

  `logical`. Not available for MaxQuant output, where
  `Leading.razor.protein` always holds a single accession, so there is
  nothing to filter. `TRUE` raises an error pointing at `proteotypic`.

- proteotypic:

  `logical`. Filter out features whose peptide sequence is found in more
  than one protein.

- filter_contaminant:

  `logical`. Filter out features which match a contaminant protein.

- contaminant_proteins:

  `character vector`. The protein IDs form the contaminant proteins

- filter_associated_contaminant:

  `logical`. Filter out features which match a contaminant associated
  protein.

- remove_no_quant:

  `logical`. Remove features with no quantification

## Value

Returns a `SummarisedExperiment` with the filtered MaxQuant output.

## Details

MaxQuant assigns every feature a single razor protein, so unlike the
other search engines it never reports a tie between protein groups and
there is nothing for `unique_master` to filter. The ambiguity is instead
recorded per peptide, in `Unique..Proteins.`, and `proteotypic = TRUE`
is the equivalent filter. See
[`vignette("gotcha_peptide_to_protein")`](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/gotcha_peptide_to_protein.md).

**Associated contaminant proteins** are proteins which have at least one
feature shared with a contaminant protein. It has been observed that the
contaminant fasta files often do not contain all possible contaminant
proteins e.g. some features can be assigned to a keratin which is not in
the provided contaminant database.

In the example below, using `filter_associated_contaminant = TRUE` will
filter out f2 and f3 in addition to f1, regardless of the value in the
Leading.razor.protein column.

    feature  Proteins                          Leading.razor.protein
    f1       protein1, protein2, contaminant,  protein1
    f2       protein1, protein3                protein3
    f3       protein2                          protein2
