# Filter Proteome Discoverer DDA output

This function filters the output .txt files (peptide groups or PSMs)
from Proteome Discoverer for DDA, based on various criteria:

1.  Remove features without a master protein

2.  Remove features without a unique master protein (i.e.
    Number.of.Protein.Groups == 1)

3.  Remove features which are not proteotypic (i.e. Number.of.Proteins
    == 1)

4.  Remove features matching a contaminant protein

5.  Remove features matching any protein associated with a contaminant
    protein (see below)

6.  Remove features without quantification values

## Usage

``` r
filter_features_pd_dda(
  obj,
  master_protein_col = "Master.Protein.Accessions",
  protein_col = "Protein.Accessions",
  unique_master = TRUE,
  proteotypic = FALSE,
  filter_contaminant = TRUE,
  contaminant_proteins = NULL,
  crap_proteins = NULL,
  filter_associated_contaminant = TRUE,
  remove_no_quant = TRUE,
  cont_string = "Cont_"
)
```

## Arguments

- obj:

  `SummarisedExperiment` containing output from Proteome Discoverer. Use
  [`readQFeatures`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html)
  to read in .txt file

- master_protein_col:

  `string`. Name of column containing master proteins.

- protein_col:

  `string`. Name of column containing all protein matches.

- unique_master:

  `logical`. Filter out features where the master protein column does
  not resolve to a single protein accession.

- proteotypic:

  `logical`. Filter out features whose peptide sequence is found in more
  than one protein.

- filter_contaminant:

  `logical`. Filter out features which match a contaminant protein.

- contaminant_proteins:

  `character vector`. The protein IDs form the contaminant proteins

- crap_proteins:

  `character vector`. Same as contaminant_proteins. Available for
  backwards compatibility. Default is NULL. If both contaminant_proteins
  and crap_proteins are set, an error is thrown.

- filter_associated_contaminant:

  `logical`. Filter out features which match a contaminant associated
  protein.

- remove_no_quant:

  `logical`. Remove features with no quantification

- cont_string:

  `string`. string to search for contaminants

## Value

Returns a `SummarisedExperiment` with the filtered Proteome Discoverer
output.

## Details

`unique_master` and `proteotypic` are different filters. `unique_master`
asks whether Proteome Discoverer resolved the feature to a single
protein accession, so it removes features which are ambiguous between
protein groups. `proteotypic` asks whether the peptide sequence occurs
in only one protein in the database, so it also removes features whose
protein group holds several indistinguishable proteins, even where a
master was assigned. See
[`vignette("gotcha_peptide_to_protein")`](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/gotcha_peptide_to_protein.md).

**Associated contaminant proteins** are proteins which have at least one
feature shared with a contaminant protein. It has been observed that the
contaminant fasta files often do not contain all possible contaminant
proteins e.g. some features can be assigned to a keratin which is not in
the provided contaminant database.

In the example below, using `filter_associated_contaminant = TRUE` will
filter out f2 and f3 in addition to f1, regardless of the value in the
Master.Protein.Accession column.

    feature  Protein.Accessions         Master.Protein.Accessions
    f1       protein1, protein2, contaminant,  protein1,
    f2       protein1, protein3         protein3,
    f3       protein2                   protein2

## Examples

``` r
# load PD PSM-level output
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_clock,
  colData = tmt_clock_design,
  quantCols = rownames(tmt_clock_design),
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

# extract the accessions from the contaminant FASTA, in both the prefixed
# and bare forms, since the search may not have renamed its entries
contaminant_fasta <- system.file(
  "extdata", "0602_Universal_Contaminants.fasta.gz", package = "biomasslmb")

contaminant_accessions <- get_contaminant_fasta_accessions(contaminant_fasta)
contaminant_accessions <- c(contaminant_accessions,
                            sub("^Cont_", "", contaminant_accessions))

# remove contaminants and PSMs without a unique master protein
psms_filtered <- filter_features_pd_dda(
  obj = tmt_qf[["psms_raw"]],
  contaminant_proteins = contaminant_accessions,
  filter_contaminant = TRUE,
  filter_associated_contaminant = TRUE,
  unique_master = TRUE)
#> Filtering data...
#> 11281 features found from 790 master proteins => Input
#> 762 contaminant proteins supplied
#> 4949 proteins identified as 'contaminant associated'
#> 9785 features found from 746 master proteins => contaminant features removed
#> 8956 features found from 727 master proteins => associated contaminant features removed
#> 8956 features found from 727 master proteins => PD-labelled 'Contaminants' removed
#> 8953 features found from 726 master proteins => features without a master protein removed
#> 8822 features found from 648 master proteins => features with non-unique master proteins removed
#> 6312 features found from 606 master proteins => features without quantification removed

c(before = nrow(tmt_qf[["psms_raw"]]), after = nrow(psms_filtered))
#> before  after 
#>  11281   6312 
```
