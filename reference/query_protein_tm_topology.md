# Query UniProt to determine protein transmembrane-domains and topology

Given a set of UniProt IDs, this function returns the information
regarding transmembrane domains and topology.

## Usage

``` r
query_protein_tm_topology(uniprotIDs, verbosity = 0)
```

## Arguments

- uniprotIDs:

  `character vector` Uniprot IDs

- verbosity:

  `integer` Verbosity level for uniprotREST::uniprot_map

## Value

`data.frame` object.

## Examples

``` r
uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
query_protein_tm_topology(uniprotIDs)
#>   UniprotID Length                                           Intramembrane
#> 1    O76024    890                                                        
#> 2    Q03135    178 INTRAMEM 105..125; /note=Helical; /evidence=ECO:0000255
#> 3    Q96T23   1441                                                        
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Transmembrane
#> 1 TRANSMEM 314..334; /note=Helical; /evidence=ECO:0000255; TRANSMEM 340..360; /note=Helical; /evidence=ECO:0000255; TRANSMEM 402..422; /note=Helical; /evidence=ECO:0000255; TRANSMEM 427..447; /note=Helical; /evidence=ECO:0000255; TRANSMEM 465..485; /note=Helical; /evidence=ECO:0000255; TRANSMEM 496..516; /note=Helical; /evidence=ECO:0000255; TRANSMEM 529..549; /note=Helical; /evidence=ECO:0000255; TRANSMEM 563..583; /note=Helical; /evidence=ECO:0000255; TRANSMEM 589..609; /note=Helical; /evidence=ECO:0000255; TRANSMEM 632..652; /note=Helical; /evidence=ECO:0000255; TRANSMEM 870..890; /note=Helical; /evidence=ECO:0000255
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#>                                                                                                       Topological.domain
#> 1                                                                TOPO_DOM 653..869; /note=Lumenal; /evidence=ECO:0000255
#> 2 TOPO_DOM 2..104; /note=Cytoplasmic; /evidence=ECO:0000255; TOPO_DOM 126..178; /note=Cytoplasmic; /evidence=ECO:0000255
#> 3                                                                                                                       
```
