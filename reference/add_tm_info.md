# Add transmembrane domain details

Parse the Transmembrane column from the output of
query_protein_tm_topology to add columns decscribing the transmembrane
domains.

## Usage

``` r
add_tm_info(query_result)
```

## Arguments

- query_resul:

  `data.frame` output of query_protein_tm_topology

## Value

`data.frame` object.

## Examples

``` r
uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
query_result <- query_protein_tm_topology(uniprotIDs)
#> Waiting 4s for retry backoff ■■■■■■■■■                       
#> Waiting 4s for retry backoff ■■■■■■■■■■■■■■■■■■■             
#> Waiting 4s for retry backoff ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  
query_result_parsed <- add_tm_info(query_result)
```
