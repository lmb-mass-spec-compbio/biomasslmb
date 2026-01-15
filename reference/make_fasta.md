# Make a FASTA using UniProt accessions

Given a vector of UniProt accessions, this function will:

1.  Download the sequences

2.  Print the current UniProt release (you should put this in the FASTA
    file name)

3.  (Optional) add cRAP numbers to the FASTA headers for each sequence

4.  Save the sequences in a FASTA file

## Usage

``` r
make_fasta(
  accessions,
  file,
  is_crap = FALSE,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- accessions:

  `character vector`, the UniProt accessions to use

- file:

  `character`, filepath to save the fasta to e.g. `"crap.fasta"`

- is_crap:

  `logical`, Is the output going to be a cRAP database? If `TRUE`
  cRAP001, cRAP002, etc. is appended to the sequence headers in the
  FASTA file. Default is `FALSE`

- overwrite:

  `logical`, if the FASTA file already exists should it be overwritten?
  Default is `FALSE`

- verbose:

  `logical`, should the function send messages to the console? Default
  is `TRUE`

## Value

Returns a FASTA file saved to disk at the specified file path.

## Examples

``` r
# specify some UniProt accessions
crap <- get_ccp_crap()
#> Error in get_ccp_crap(): could not find function "get_ccp_crap"

if (FALSE) { # \dontrun{
make_fasta(crap, "2021-01_cRAP.fasta")
} # }
```
