# Obtain protein transmembrane-domains and topology information

Given a set of UniProt IDs, this function queries UniProt to obtain the
information regarding transmembrane domains (TMDs) and topology and then
parses this information to define additional columns describining the
TMDs and topology

## Usage

``` r
get_protein_tm_topology(uniprotIDs, verbosity = 0)
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
get_protein_tm_topology(uniprotIDs)
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
#>   n_tms                                    tm_start
#> 1    11 314;340;402;427;465;496;529;563;589;632;870
#> 2     0                                        <NA>
#> 3     0                                        <NA>
#>                                        tm_end                        tm_length
#> 1 334;360;422;447;485;516;549;583;609;652;890 20;20;20;20;20;20;20;20;20;20;20
#> 2                                        <NA>                             <NA>
#> 3                                        <NA>                             <NA>
#>                                                                                  tm_types
#> 1 Helical;Helical;Helical;Helical;Helical;Helical;Helical;Helical;Helical;Helical;Helical
#> 2                                                                                    <NA>
#> 3                                                                                    <NA>
#>   tm_types_set                topology      n_term      c_term topology_start
#> 1      Helical                 Lumenal     Lumenal     Lumenal            653
#> 2         <NA> Cytoplasmic;Cytoplasmic Cytoplasmic Cytoplasmic          2;126
#> 3         <NA>                    <NA>        <NA>        <NA>           <NA>
#>   topology_end topology_length
#> 1          869             216
#> 2      104;178          102;52
#> 3         <NA>            <NA>
```
