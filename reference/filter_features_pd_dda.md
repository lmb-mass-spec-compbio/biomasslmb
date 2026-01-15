# Filter Proteome Discoverer DDA output

This function filters the output .txt files (peptide groups or PSMs)
from Proteome Discoverer for DDA, based on various criteria:

1.  Remove features without a master protein

2.  Remove features without a unique master protein (i.e.
    Number.of.Protein.Groups == 1)

3.  Remove features matching a contaminant protein

4.  Remove features matching any protein associated with a contaminant
    protein (see below)

5.  Remove features without quantification values

## Usage

``` r
filter_features_pd_dda(
  obj,
  master_protein_col = "Master.Protein.Accessions",
  protein_col = "Protein.Accessions",
  unique_master = TRUE,
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
  [`readQFeatures`](https://rdrr.io/pkg/QFeatures/man/readQFeatures.html)
  to read in .txt file

- master_protein_col:

  `string`. Name of column containing master proteins.

- protein_col:

  `string`. Name of column containing all protein matches.

- unique_master:

  `logical`. Filter out features without a unique master protein.

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
if (FALSE) { # \dontrun{

#### PSMs.txt example ####
# load PD PSMs.txt output
tmt_qf <- readQFeatures(assayData = psm_tmt_total,
 quantCols = 36:45,
 name = "psms_raw")

# extract the UniProt accessions from the contaminant FASTA headers
contaminant_accessions <- get_contaminant_fasta_accessions(contaminant_fasta_inf)

# filter the PSMs
psm2 <- filter_features_pd_dda(
  obj = tmt_qf[['psms_raw']],
  master_protein_col = "Master.Protein.Accessions",
  protein_col = "Protein.Accessions",
  unique_master = TRUE,
  TMT = TRUE,
  filter_contaminant = TRUE,
  contaminant_proteins = contaminant_accessions,
  filter_associated_contaminant = TRUE
)


} # }
```
