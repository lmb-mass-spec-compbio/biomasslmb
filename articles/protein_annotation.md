# Protein annotation: UniProt details, GO terms and topology

Search engine output identifies a protein by its UniProt accession and,
at best, whatever description string was in the FASTA. Almost everything
downstream needs more than that: gene symbols to label a volcano plot,
GO terms to run an enrichment test, transmembrane topology to ask
whether the proteins gained in a membrane fraction are the ones you
would expect. That annotation has to be fetched from UniProt.

It is worth fetching it as its own step, before any QC or filtering. The
annotations depend only on which proteins were *detected*, not on how
the data are subsequently filtered, normalised or summarised to protein
level. Retrieving them once and caching the result means the rest of the
pipeline can be re-run — after changing a filtering threshold, say —
without querying UniProt again.

The queries below are shown as you would run them, but the results
displayed were retrieved from UniProt release 2026_02 and cached, so
that building this vignette does not depend on UniProt being reachable.
Running the same calls yourself will return the current release’s
annotations, which may differ.

## Load required packages

``` r

library(QFeatures)
library(biomasslmb)
library(dplyr)

lfq_qf <- biomasslmb::lfq_qf
```

## Collecting the accessions to annotate

We use the peptide-level Proteome Discoverer (PD) output from the
[LFQ-DDA
vignette](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/LFQ_DDA_Peptide_QC_Summarisation.md)
— reading the raw search output rather than the processed `QFeatures`
object, since the point of annotating first is that this step does not
depend on the processing.

``` r

pep_inf <- system.file(
  "extdata", "lfq_dda_pd_PeptideGroups.txt",
  package = "biomasslmb"
)

infdf <- read.delim(pep_inf)
```

PD reports a `Master.Protein.Accessions` value for each peptide. Where a
peptide cannot be assigned unambiguously to one protein — most often
because it is shared between near-identical paralogues or isoforms — PD
reports several accessions in a single field, joined by `"; "`. We
therefore need two sets of identifiers: the protein groups exactly as
the search engine reports them, which is what we will eventually join
back onto the data, and the individual accessions, which is what UniProt
can actually be queried with.

``` r

protein_groups <- unique(infdf$Master.Protein.Accessions)

accessions <- protein_groups %>%
  strsplit('; ') %>%
  unlist() %>%
  unique()

length(protein_groups)
#> [1] 670
length(accessions)
#> [1] 694
```

The second number is the larger, and the difference is the number of
extra accessions contributed by multi-protein groups.

Other search engines report the same thing under a different column name
and, importantly, a different separator: MaxQuant uses
`Leading.razor.protein` with `";"` (no space), and DIA-NN uses
`Protein.Group`, also with `";"`. The separator matters, because it is
what you split on here and what you pass to
[`collapse_uniprot_details_multi_accession()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/collapse_uniprot_details_multi_accession.md)
below.

## UniProt protein details

[`get_uniprot_details()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
retrieves the entry name, protein name(s), gene name(s), organism and
sequence length for each accession. It also derives `Gene.Names.First`,
the first of the possibly space-separated gene names UniProt returns,
which is the primary symbol and the one you want for plot labels.

``` r

uniprot2details <- biomasslmb::get_uniprot_details(accessions)
```

``` r

head(as.data.frame(uniprot2details))
#>    UniprotID       Entry.Name   Reviewed
#> 1 A0A075B6R9      KVD24_HUMAN   reviewed
#> 2 A0A075B6S2      KVD29_HUMAN   reviewed
#> 3 A0A075B6S6      KVD30_HUMAN   reviewed
#> 4 A0A075B6Z2 A0A075B6Z2_HUMAN unreviewed
#> 5     A2A3N6      PIPSL_HUMAN   reviewed
#> 6     E9PAV3      NACAM_HUMAN   reviewed
#>                                                                                                          Protein.names
#> 1                                                          Probable non-functional immunoglobulin kappa variable 2D-24
#> 2                                                                                  Immunoglobulin kappa variable 2D-29
#> 3                                                                                  Immunoglobulin kappa variable 2D-30
#> 4                                                                                     T cell receptor alpha joining 56
#> 5                                                              Putative PIP5K1A and PSMD4-like protein (PIP5K1A-PSMD4)
#> 6 Nascent polypeptide-associated complex subunit alpha, muscle-specific form (Alpha-NAC, muscle-specific form) (skNAC)
#>      Gene.Names             Organism Length Gene.Names.First Annotation.Source
#> 1     IGKV2D-24 Homo sapiens (Human)    120        IGKV2D-24  UniProtKB (live)
#> 2     IGKV2D-29 Homo sapiens (Human)    120        IGKV2D-29  UniProtKB (live)
#> 3     IGKV2D-30 Homo sapiens (Human)    120        IGKV2D-30  UniProtKB (live)
#> 4        TRAJ56 Homo sapiens (Human)     21           TRAJ56  UniProtKB (live)
#> 5 PIPSL PSMD4P2 Homo sapiens (Human)    862            PIPSL  UniProtKB (live)
#> 6          NACA Homo sapiens (Human)   2078             NACA  UniProtKB (live)
```

Two complications are handled for you.

**Demerged accessions.** An accession that has since been split into
several current entries — because the original was found to represent
more than one distinct protein — returns several rows. These are
collapsed into one row per queried accession, with the differing fields
joined by `;`.

**Retired accessions.** A search is run against a FASTA snapshot of
UniProt taken at search time, which may be months or years old by the
time you annotate. Accessions deleted from UniProtKB since then come
back from a plain ID mapping with a blank gene name and
`Protein.names == "deleted"`, even though the protein may be a perfectly
real hit. For those,
[`get_uniprot_details()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
follows the entry to UniParc and recovers the last known annotation. The
`Annotation.Source` column records which route each row took:

``` r

table(uniprot2details$Annotation.Source)
#> 
#> UniProtKB (live) 
#>              663
```

None of the accessions in this dataset are retired, so all of them
resolved directly. Where one has been, the fallback looks like this —
`A0A5G2QPJ4` is a pig entry that has since been deleted, and its gene
and protein names are recovered from UniParc:

``` r

uniparc_example <- biomasslmb::get_uniprot_details(c('O76024', 'A0A5G2QPJ4'))
```

``` r

uniparc_example %>%
  select(UniprotID, Gene.Names.First, Protein.names, Annotation.Source) %>%
  as.data.frame()
#>    UniprotID Gene.Names.First            Protein.names
#> 1 A0A5G2QPJ4            CADM1 Cell adhesion molecule 1
#> 2     O76024             WFS1                Wolframin
#>                   Annotation.Source
#> 1 UniParc (retired UniProtKB entry)
#> 2                  UniProtKB (live)
```

### Checking the result is complete

UniProt’s ID mapping service occasionally returns a truncated cached
response, in which a few hundred accessions silently lose their
annotation with no error.
[`get_uniprot_details()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
guards against this with `check_complete = 'error'` (the default), which
compares the number of accessions submitted against the number mapped
plus the number UniProt reports as failed, and raises an error if they
disagree.

That check passing does not mean every accession got an annotation —
accessions UniProt explicitly reports as failed are simply absent from
the result. Some are missing here:

``` r

unannotated <- setdiff(accessions, uniprot2details$UniprotID)

length(unannotated)
#> [1] 31
head(unannotated)
#> [1] "Cont_Q86YZ3" "Cont_P02533" "Cont_P13645" "Cont_P08779" "Cont_P60712"
#> [6] "Cont_Q2KJD0"
```

These are the contaminant proteins. The search database prefixes
contaminant accessions with `Cont_`, which is not a UniProt accession,
so UniProt cannot map them. That is harmless as long as you are aware of
it, because contaminants are removed during QC anyway — but it is the
reason the annotation table has fewer rows than the accession vector you
passed in, and it is worth checking that the accessions that dropped out
are ones you expected to lose.

## Annotating protein groups

`uniprot2details` has one row per individual accession, so it cannot be
joined directly against data keyed by PD’s multi-accession
`Master.Protein.Accessions` strings.
[`collapse_uniprot_details_multi_accession()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/collapse_uniprot_details_multi_accession.md)
builds that table: it splits each group into its accessions, looks each
one up, and re-collapses to one row per group, with every constituent
accession contributing one value to every column so that the parts stay
positionally aligned across columns.

``` r

uniprot2groups <- biomasslmb::collapse_uniprot_details_multi_accession(
  uniprot2details, protein_groups, sep = '; '
)

uniprot2groups %>%
  filter(grepl('; ', UniprotID)) %>%
  select(UniprotID, Gene.Names.First) %>%
  head() %>%
  as.data.frame()
#>                     UniprotID     Gene.Names.First
#> 1      A0A075B6R9; A0A075B6S2 IGKV2D-24; IGKV2D-29
#> 2         Cont_P01966; P69905                HBA1;
#> 3         Cont_P02769; P02768                  ALB
#> 4         Cont_P60712; P62736                ACTA2
#> 5 Cont_P60712; P62736; Q562R1        ACTA2; ACTBL2
#> 6         Cont_P60712; Q6S8J3                POTEE
```

Pass the separator your search engine uses: `sep = ';'` for MaxQuant and
DIA-NN.

A group is only returned if every one of its accessions was annotated,
so the contaminant groups drop out here too:

``` r

length(protein_groups)
#> [1] 670
nrow(uniprot2groups)
#> [1] 631
```

## Adding the annotations to a `QFeatures` object

The annotations belong in the assay’s `rowData`, joined on
`Master.Protein.Accessions` — the column the protein-level row names are
taken from. A left join keeps every row of the assay, in its existing
order, and leaves an `NA` wherever a group has no annotation.
`rowData()` returns a `DataFrame` rather than a `data.frame`, so it
needs converting before `dplyr` will accept it, and select down to the
columns you want so that the join does not bring the whole annotation
table across.

``` r

protein_se <- lfq_qf[['protein']]

rowData(protein_se) <- rowData(protein_se) %>%
  as.data.frame() %>%
  left_join(
    uniprot2groups %>% select(UniprotID, Gene.Names.First, Protein.names),
    by = c('Master.Protein.Accessions' = 'UniprotID')
  )

lfq_qf[['protein']] <- protein_se

rowData(lfq_qf[['protein']])[1:5, c('Gene.Names.First', 'Protein.names')] %>%
  as.data.frame()
#>        Gene.Names.First
#> O00159            MYO1C
#> O00422            SAP18
#> O00425          IGF2BP3
#> O14744            PRMT5
#> O15027           SEC16A
#>                                                                                                                                                                                                                                                                                                          Protein.names
#> O00159                                                                                                                                                                                                                                                      Unconventional myosin-Ic (Myosin I beta) (MMI-beta) (MMIb)
#> O00422                                                                                                                                            Histone deacetylase complex subunit SAP18 (18 kDa Sin3-associated polypeptide) (2HOR0202) (Cell growth-inhibiting gene 38 protein) (Sin3-associated polypeptide p18)
#> O00425                                                                                                 Insulin-like growth factor 2 mRNA-binding protein 3 (IGF2 mRNA-binding protein 3) (IMP-3) (IGF-II mRNA-binding protein 3) (KH domain-containing protein overexpressed in cancer) (hKOC) (VICKZ family member 3)
#> O14744 Protein arginine N-methyltransferase 5 (PRMT5) (EC 2.1.1.320) (72 kDa ICln-binding protein) (Histone-arginine N-methyltransferase PRMT5) (Jak-binding protein 1) (Shk1 kinase-binding protein 1 homolog) (SKB1 homolog) (SKB1Hs) [Cleaved into: Protein arginine N-methyltransferase 5, N-terminally processed]
#> O15027                                                                                                                                                                                                                                                       Protein transport protein Sec16A (SEC16 homolog A) (p250)
```

Reach for
[`left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
rather than [`merge()`](https://rdrr.io/r/base/merge.html) here:
[`merge()`](https://rdrr.io/r/base/merge.html) sorts its output by the
key and drops unmatched rows unless told otherwise, either of which puts
`rowData` out of step with the rows of the assay.

Every group matched here, since the contaminants were removed from
`lfq_qf` during QC:

``` r

sum(!rowData(protein_se)$Master.Protein.Accessions %in% uniprot2groups$UniprotID)
#> [1] 0
```

DIA-NN output is the exception to all of this: it carries `Genes` and
`Protein.Names` columns of its own, taken from the search FASTA, and the
[exploration and statistical
testing](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/exploration_and_statistical_testing.md)
vignette uses them directly. They are convenient, but they are frozen at
whatever release the FASTA came from, and they do not include GO terms
or topology, so a UniProt query is still needed for anything beyond
labelling.

## GO terms

[`get_go_terms()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_go_terms.md)
returns one row per protein per GO term.

``` r

go_res <- biomasslmb::get_go_terms(accessions)
```

``` r

head(as.data.frame(go_res))
#>   UNIPROTKB                           GO_desc      GO.ID
#> 1    O14744                         chromatin GO:0000785
#> 2    O14744                         cytoplasm GO:0005737
#> 3    O14744                           cytosol GO:0005829
#> 4    O14744                   Golgi apparatus GO:0005794
#> 5    O14744 histone methyltransferase complex GO:0035097
#> 6    O14744                       methylosome GO:0034709
```

GO terms form a hierarchy running from specific to general, and a
protein annotated with a specific term is implicitly also described by
every broader term above it. Enrichment tools that do not walk the
hierarchy themselves will therefore miss enrichment of a broad term
simply because the proteins were only ever annotated with its more
specific children. `expand_terms = TRUE` adds all ancestor terms
explicitly, which is what you want when the result is destined for GO
over-representation testing.

``` r

go_res_all <- biomasslmb::get_go_terms(accessions, expand_terms = TRUE)
```

``` r

head(as.data.frame(go_res_all))
#>   UNIPROTKB      GO.ID                            TERM ONTOLOGY
#> 1    P16403 GO:0000018 regulation of DNA recombination       BP
#> 2    Q9Y230 GO:0000018 regulation of DNA recombination       BP
#> 3    P43246 GO:0000018 regulation of DNA recombination       BP
#> 4    P61160 GO:0000018 regulation of DNA recombination       BP
#> 5    Q9UNS1 GO:0000018 regulation of DNA recombination       BP
#> 6    P52292 GO:0000018 regulation of DNA recombination       BP
```

Expansion is much the slower call, and it changes the shape of the
output as well as its size: the unexpanded result has a `GO_desc`
column, while the expanded result has `TERM` and `ONTOLOGY` columns
looked up from `GO.db`. The `UNIPROTKB` and `GO.ID` columns are common
to both, and are the two columns that the enrichment functions take as
their `gene2cat` mapping.

``` r

nrow(go_res)
#> [1] 12714
nrow(go_res_all)
#> [1] 79410
```

Not every protein has GO annotation — unreviewed entries frequently have
none — so the mapping covers fewer proteins than you queried:

``` r

length(unique(go_res_all$UNIPROTKB))
#> [1] 658
```

Note that GO terms are keyed by individual accession, not by protein
group, so a multi-accession group needs a decision about which of its
accessions’ terms to use before enrichment testing.

## Transmembrane domains and topology

[`get_protein_tm_topology()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_protein_tm_topology.md)
retrieves UniProt’s transmembrane and topological domain features and
parses them into columns describing each protein’s membrane topology:
the number of transmembrane segments, their start and end positions, and
which side of the membrane the N- and C-termini and each loop are on.

``` r

tm_topology <- biomasslmb::get_protein_tm_topology(
  c('P48962', 'P51881', 'P62827')
)
```

``` r

tm_topology %>%
  select(UniprotID, Length, n_tms, tm_start, tm_end, n_term, c_term) %>%
  as.data.frame()
#>   UniprotID Length n_tms             tm_start                tm_end
#> 1    P48962    298     6 8;75;110;179;211;274 37;99;130;199;231;291
#> 2    P51881    298     6 8;75;110;179;211;274 37;99;130;199;231;291
#> 3    P62827    216     0                 <NA>                  <NA>
#>                        n_term                      c_term
#> 1 Mitochondrial intermembrane Mitochondrial intermembrane
#> 2 Mitochondrial intermembrane Mitochondrial intermembrane
#> 3                        <NA>                        <NA>
```

The two ADP/ATP translocases are six-pass mitochondrial inner membrane
carriers; the third protein is soluble, so its topology columns are
`NA`. Proteins whose segments are all beta-stranded are also returned as
`NA`, since the start/end parsing assumes alpha-helical segments.

## Caching the annotations

Save the results so the rest of the analysis can read them rather than
re-querying, and record which UniProt release they came from — releases
are published roughly every eight weeks, and annotations retrieved
months apart are not necessarily identical.

``` r

biomasslmb::check_uniprot_release()

saveRDS(uniprot2details, 'uniprot2details.rds')
saveRDS(uniprot2groups, 'uniprot2groups.rds')
saveRDS(go_res_all, 'go_res_all.rds')
```

## Things that catch people out

**The separator differs between search engines.** PD uses `"; "`,
MaxQuant and DIA-NN use `";"`. Splitting on the wrong one leaves the
multi-accession groups intact, so they are passed to UniProt as single
unrecognised identifiers and come back unannotated.

**Contaminant accessions carry a prefix.** `Cont_P13646` is not a
UniProt accession and will not map. This is the same prefix mismatch
that silently leaves contaminants in the data if you filter on
unprefixed accessions.

**A complete-looking result is not a complete result.** `check_complete`
catches a truncated response, but accessions UniProt reports as failed
are absent from the output with no warning. Compare the accessions in
against the accessions out.

**Search engine gene names are frozen.** The gene and protein names in
the search output come from the FASTA, and so reflect UniProt as it was
at search time. Where they disagree with a fresh query, the fresh query
is the current answer — which is also why it is worth recording the
release the annotations came from.

## Where to go next

With annotations in hand, the acquisition-specific QC and summarisation
vignette for your data is the next step — see [Getting
started](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/biomasslmb.md)
for which one applies. The GO mapping produced here is the input to
functional enrichment testing, performed after differential abundance
testing in [data exploration and statistical
testing](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/exploration_and_statistical_testing.md).

## Getting help

If your experiment does not fit what this article assumes, or a result
looks wrong, please get in touch with Tom Smith (<tsmith@mrclmb.ac.uk>)
rather than guessing — it is far easier to help while an analysis is in
progress than to unpick a decision afterwards, and easier still before
the samples are run. [Getting
started](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/biomasslmb.md)
sets out which article applies to which experiment.
