# Append sequences to end of a FASTA

This function is used to add sequences from a FASTA file (file1) onto
the end of another FASTA file (file2).

If file2 is a cRAP database, then you can optionally add cRAP numbers to
the headers of file1, starting at the desired number e.g. cRAP127.

## Usage

``` r
append_fasta(file1, file2, is_crap = FALSE, crap_start = 1)
```

## Arguments

- file1:

  `character`, file path of FASTA to append

- file2:

  `character`, file path of FASTA to append to

- is_crap:

  `logical`, should cRAP numbers e.g. cRAP001, cRAP002, etc. be added to
  sequence headers of file1? Default is `FALSE`

- crap_start:

  `numeric`, what number should the cRAP00x start at? Default is 1.

## Value

Overwrites FASTA file2 with some more sequences added to the end.

## Examples

``` r
# Add some commercial protease sequences onto the end of CCP cRAP FASTA
if (FALSE) { # \dontrun{
append_fasta(
  file1 = "commercial-proteases.fasta",
  file2 = "2021-01_CCP_cRAP.fasta",
  is_crap = TRUE,
  crap_start = 128
)
} # }
```
