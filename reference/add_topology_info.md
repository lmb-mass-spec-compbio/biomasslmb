# Add topology details

Parse the Topological.domain column from the output of
query_protein_tm_topology to add columns decscribing the topology

## Usage

``` r
add_topology_info(result)
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
query_result_parsed <- add_topology_info(query_result)
```
