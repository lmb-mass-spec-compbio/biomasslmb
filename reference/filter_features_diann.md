# Filter DIA-NN output

This function filters the precursor .txt file from DIA-NN, based on
various criteria:

1.  Remove features without a master protein (Protein.Group column)

2.  Remove features without a unique master protein

3.  Remove features matching a contaminant protein

4.  Remove features matching any protein associated with a contaminant
    protein (see below)

5.  Remove features without quantification values

## Usage

``` r
filter_features_diann(
  obj,
  master_protein_col = "Protein.Group",
  protein_col = "Protein.Ids",
  unique_master = TRUE,
  filter_contaminant = TRUE,
  contaminant_proteins = NULL,
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
