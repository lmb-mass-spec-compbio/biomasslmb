# Check the current UniProt release

UniProt releases are published approximately every 8 weeks. This
function checks what the current UniProt release is.

## Usage

``` r
check_uniprot_release()
```

## Value

Returns a `character`, the current release number in the format YYYY_XX
where YYYY is the calendar year and XX a 2-digit number that is
incremented for each release of a given year, e.g. 2010_01, 2010_02,
etc.

## Examples

``` r
# print release number to console
check_uniprot_release()
#> [1] "2025_04"

# save release number and use in e.g. a file name
rls <- check_uniprot_release()

paste0("folder/filename_", rls, ".fasta")
#> [1] "folder/filename_2025_04.fasta"
```
