# Pitfall: contaminants and protein-level FDR

Every sample carries protein that did not come from the experiment:
keratin from skin and dust, trypsin from the digest, albumin from serum
in the culture medium. These are searched against deliberately — a
contaminant database gives their spectra somewhere correct to go,
instead of leaving them to be force-matched onto real proteins — and
then have to be removed before anything is quantified.

The removal is where things go quietly wrong, because the QC vignettes
apply *three* overlapping defences against contaminants and, on a
typical export, any one of them is sufficient. When one is broken you do
not find out, because the other two cover for it. This article is about
what each one actually does, and what an export looks like when the
cover is not there.

It also covers protein-level FDR, which is a different kind of filter
applied at the same stage and for a related reason: not everything the
search engine reports is something you should quantify.

## Load required packages

``` r

library(QFeatures)
library(biomasslmb)
library(dplyr)
```

``` r

pep_inf <- system.file("extdata", "lfq_dda_pd_PeptideGroups.txt",
                       package = "biomasslmb")
infdf <- read.delim(pep_inf)

abundance_cols_ix <- grep('^Abundance', colnames(infdf))
exp_design <- data.frame(
  quantCols = colnames(infdf)[abundance_cols_ix],
  Condition = rep(c('WT', 'Mutant'), each = 3),
  Replicate = rep(1:3, times = 2))

read_peptides <- function(df) {
  qf <- readQFeatures(assayData = df, quantCols = abundance_cols_ix,
                      colData = exp_design, name = 'peptides')
  qf <- sync_coldata(qf, 'peptides')
  qf[['peptides']]
}

peptides <- read_peptides(infdf)
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.
nrow(peptides)
#> [1] 3544
```

## The three defences

**The search engine’s own flag.** Proteome Discoverer writes a
`Contaminant` column, MaxQuant a `Potential.contaminant` column, marking
anything matched to the contaminant database.
[`filter_features_pd_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_pd_dda.md)
and
[`filter_features_mq_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_mq_dda.md)
use it automatically.

**The accession prefix.** Contaminant databases are normally built with
the entries renamed, so a contaminant accession reads `Cont_P02538` in a
PD search or `CON__P02538` in a MaxQuant one. The filtering functions
grep for that prefix — `cont_string`, which defaults to `Cont_` for the
PD functions and `CON__` for the MaxQuant ones.

**An explicit accession list.**
[`get_contaminant_fasta_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_contaminant_fasta_accessions.md)
parses the contaminant FASTA and returns the accessions it holds;
[`get_maxquant_cont_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_maxquant_cont_accessions.md)
does the same for MaxQuant’s bundled database. That list is passed as
`contaminant_proteins`.

``` r

contaminant_accessions <- get_contaminant_fasta_accessions(system.file(
  "extdata", "0602_Universal_Contaminants.fasta.gz", package = "biomasslmb"))

maxquant_accessions <- get_maxquant_cont_accessions()
#> No contaminants fasta supplied, defaulting to contaminants.fasta.gz from biomasslmb

c(fasta = length(contaminant_accessions), maxquant = length(maxquant_accessions))
#>    fasta maxquant 
#>      381      246
head(contaminant_accessions, 3)
#> [1] "Cont_P00722" "Cont_P09870" "Cont_P30879"
head(maxquant_accessions, 3)
#> [1] "CON__P00761" "CON__Q32MB2" "CON__P19013"
```

Note that the two helpers return accessions in **different formats**,
matching the database each search used. A list from one is useless
against output from the other, and nothing checks that you passed the
right one.

## On a normal export, the list does nothing

``` r

filter_contaminants <- function(obj, accessions) {
  nrow(filter_features_pd_dda(
    obj,
    contaminant_proteins = accessions,
    filter_contaminant = TRUE,
    filter_associated_contaminant = TRUE,
    unique_master = FALSE,
    remove_no_quant = TRUE))
}

data.frame(
  accession_list = c('correct (from the search database)',
                     'MaxQuant list — wrong format entirely'),
  peptides_retained = c(
    filter_contaminants(peptides, c(contaminant_accessions,
                                    sub('^Cont_', '', contaminant_accessions))),
    filter_contaminants(peptides, maxquant_accessions)))
#>                          accession_list peptides_retained
#> 1    correct (from the search database)              2325
#> 2 MaxQuant list — wrong format entirely              2325
```

Identical. Supplying a completely inapplicable contaminant list changes
nothing, because the `Cont_` prefix grep and PD’s own `Contaminant`
column between them have already caught every contaminant in the file.
If your pipeline has a broken contaminant list, this is what it looks
like: exactly like a working one.

## When the cover is not there

The prefix grep and the engine’s flag both depend on the search having
been set up with a renamed contaminant database. That is the normal case
but not the only one — appending a contaminant FASTA to the target
database without renaming its entries is an easy thing to do, and then
contaminants are reported under their ordinary UniProt accessions with
nothing to distinguish them.

Below, the same data is altered to look like that search: prefixes
stripped and the flag cleared. Nothing else changes.

``` r

unprefixed <- infdf
unprefixed$Master.Protein.Accessions <- gsub(
  'Cont_', '', unprefixed$Master.Protein.Accessions)
unprefixed$Protein.Accessions <- gsub('Cont_', '', unprefixed$Protein.Accessions)
unprefixed$Contaminant <- 'False'

unprefixed_peptides <- read_peptides(unprefixed)
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.
```

``` r

data.frame(
  accession_list = c('correct list, including bare accessions',
                     'correct list, Cont_-prefixed forms only',
                     'MaxQuant list — wrong format entirely'),
  peptides_retained = c(
    filter_contaminants(unprefixed_peptides,
                        c(contaminant_accessions,
                          sub('^Cont_', '', contaminant_accessions))),
    filter_contaminants(unprefixed_peptides, contaminant_accessions),
    filter_contaminants(unprefixed_peptides, maxquant_accessions)))
#>                            accession_list peptides_retained
#> 1 correct list, including bare accessions              2415
#> 2 correct list, Cont_-prefixed forms only              2467
#> 3   MaxQuant list — wrong format entirely              2467
```

Now the list is the only thing standing between you and the
contaminants, and getting its format wrong leaves around fifty keratin
and trypsin peptides in the data, silently.

This is why the QC vignettes pass both forms:

``` r

contaminant_accessions <- c(contaminant_accessions,
                            sub('^Cont_', '', contaminant_accessions))
```

It costs nothing when the prefixes are present and saves you when they
are not. The general principle is worth stating plainly: **check what
your contaminant accessions actually look like in your own data before
trusting the filter**, rather than assuming the helper you called
returns the format your search produced.

``` r

master_proteins <- unique(unlist(strsplit(
  infdf$Master.Protein.Accessions, '; ')))

c(prefixed_in_this_export = length(grep('^Cont_', master_proteins)),
  matching_the_fasta_list = length(intersect(master_proteins,
                                             contaminant_accessions)))
#> prefixed_in_this_export matching_the_fasta_list 
#>                      31                      31
```

## Associated contaminants

A subtler problem. A peptide can be shared between a contaminant and a
real protein — keratins have relatives in the human proteome, and
trypsin is a protease like many others. If that peptide is assigned to
the real protein, the contaminant’s signal is quantified under a name
that looks legitimate.

`filter_associated_contaminant = TRUE` handles this by collecting every
protein that appears alongside a contaminant on any peptide, and
removing features assigned to those too.

``` r

with_associated <- filter_features_pd_dda(
  peptides, contaminant_proteins = contaminant_accessions,
  filter_contaminant = TRUE, filter_associated_contaminant = TRUE,
  unique_master = FALSE, remove_no_quant = TRUE)

without_associated <- filter_features_pd_dda(
  peptides, contaminant_proteins = contaminant_accessions,
  filter_contaminant = TRUE, filter_associated_contaminant = FALSE,
  unique_master = FALSE, remove_no_quant = TRUE)

c(with_associated = nrow(with_associated),
  without_associated = nrow(without_associated),
  extra_removed = nrow(without_associated) - nrow(with_associated))
#>    with_associated without_associated      extra_removed 
#>               2325               2336                 11
```

A small number, but the peptides it removes are the ones most likely to
mislead, because they are attached to real proteins with plausible
names. It is worth knowing that this is deliberately aggressive: a real
protein that happens to share one peptide with a keratin is discarded
entirely, not just for that peptide. If a protein you care about
disappears during filtering, this is the step to check first.

## Protein-level FDR

A separate filter applied at the same stage. Search engines control the
false discovery rate at the level of peptide-spectrum matches, but a
protein identified from a handful of borderline matches can still be a
false identification even when every one of its PSMs passed. Proteome
Discoverer therefore also reports a protein-level FDR confidence, in the
protein export rather than the peptide one.

[`filter_by_protein_fdr()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_by_protein_fdr.md)
reads that file and drops peptides belonging to proteins that did not
reach high confidence.

``` r

protein_inf <- system.file("extdata", "lfq_dda_pd_Proteins.txt.gz",
                           package = "biomasslmb")

proteins <- read.delim(gzfile(protein_inf))
table(proteins$Protein.FDR.Confidence.Combined)
#> 
#> High 
#>  583
```

Every protein here is high confidence, so no peptide can be removed for
failing the FDR threshold. Watch what happens anyway.

``` r

fdr_filtered <- filter_by_protein_fdr(with_associated,
                                      protein_fdr_filename = gzfile(protein_inf))

c(before = nrow(with_associated), after = nrow(fdr_filtered))
#> before  after 
#>   2325   2159
```

[`filter_by_protein_fdr()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_by_protein_fdr.md)
does two things, and only one of them is in its name. It joins the
peptide table to the protein table on the master protein accession, and
drops any peptide whose master protein has **no row to join to** —
regardless of confidence.

``` r

master_protein <- rowData(with_associated)$Master.Protein.Accessions

c(no_matching_protein_row = sum(!master_protein %in% proteins$Accession),
  of_those_a_protein_group = sum(!master_protein %in% proteins$Accession &
                                   grepl('; ', master_protein)),
  of_those_simply_absent = sum(!master_protein %in% proteins$Accession &
                                 !grepl('; ', master_protein)))
#>  no_matching_protein_row of_those_a_protein_group   of_those_simply_absent 
#>                      166                      143                       23
```

Most are peptides whose master protein is a group of several accessions:
the protein table has one row per protein, so a `Q9Y294; O00566` master
matches nothing and the peptide goes. The remainder are proteins the
peptide table names but the protein table does not, because Proteome
Discoverer grouped them under a different representative.

Neither is unreasonable behaviour — a peptide whose protein is absent
from the protein-level output cannot have a protein-level FDR — but it
means the step removes far more than its name suggests, and how much
depends on what you did beforehand. Filtering to unique master proteins
first separates the two effects:

``` r

unique_master_only <- filter_features_pd_dda(
  peptides, contaminant_proteins = contaminant_accessions,
  filter_contaminant = TRUE, filter_associated_contaminant = TRUE,
  unique_master = TRUE, remove_no_quant = TRUE)

c(before = nrow(unique_master_only),
  after = nrow(filter_by_protein_fdr(
    unique_master_only, protein_fdr_filename = gzfile(protein_inf))))
#> before  after 
#>   2182   2159
```

**Apply it after unique-master filtering**, so that what it removes is
what it is meant to remove.

As for the FDR filtering itself: in most Proteome Discoverer exports the
protein table contains nothing but `High`, so the step does nothing,
which is why it is easy to omit without noticing. Omitting it is still
the wrong habit, for the same reason as everything else here — the
export where it would have mattered looks exactly like the export where
it did not. Check the distribution above on your own data; `Medium` and
`Low` entries mean the filter is doing real work, and skipping it means
quantifying proteins the search engine was not confident it had
identified at all.

`retain_proteins` exists for the case where a protein of specific
interest — a bait, a tagged construct — falls below the threshold and
you want to keep it knowingly rather than lose it silently.

## Summary

- Contaminant removal has three overlapping defences: the search
  engine’s flag, the accession prefix, and an explicit accession list.
  On a normal export any one suffices, so a broken list is invisible.
- [`get_contaminant_fasta_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_contaminant_fasta_accessions.md)
  and
  [`get_maxquant_cont_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_maxquant_cont_accessions.md)
  return **differently formatted** accessions. Passing both the prefixed
  and bare forms costs nothing and covers the case where the search
  database was not renamed — which is the case where the list is the
  only defence you have.
- `filter_associated_contaminant = TRUE` also removes real proteins that
  share a peptide with a contaminant. This is intentional and
  aggressive; check it first if a protein you expected has vanished.
- Protein-level FDR is a separate control from PSM-level FDR and lives
  in a separate file.
  [`filter_by_protein_fdr()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_by_protein_fdr.md)
  applies it, and usually removes nothing on that count — which is not a
  reason to omit it.
- That function also silently drops every peptide whose master protein
  has no row in the protein table, including all peptides with a
  multi-protein master. Run it **after** `unique_master = TRUE` so the
  two effects are not confused.

## Getting help

If your experiment does not fit what this article assumes, or a result
looks wrong, please get in touch with Tom Smith (<tsmith@mrclmb.ac.uk>)
rather than guessing — it is far easier to help while an analysis is in
progress than to unpick a decision afterwards, and easier still before
the samples are run. [Getting
started](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/biomasslmb.md)
sets out which article applies to which experiment.

``` r

sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] dplyr_1.2.1                 biomasslmb_0.1.0           
#>  [3] QFeatures_1.20.0            MultiAssayExperiment_1.36.2
#>  [5] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [7] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#>  [9] IRanges_2.44.0              S4Vectors_0.48.1           
#> [11] BiocGenerics_0.56.0         generics_0.1.4             
#> [13] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1        farver_2.1.2            blob_1.3.0             
#>  [4] Biostrings_2.78.0       S7_0.2.2                fastmap_1.2.0          
#>  [7] lazyeval_0.2.3          XML_3.99-0.24           digest_0.6.39          
#> [10] lifecycle_1.0.5         cluster_2.1.8.2         ProtGenerics_1.42.0    
#> [13] survival_3.8-6          KEGGREST_1.50.0         RSQLite_3.53.3         
#> [16] magrittr_2.0.5          genefilter_1.92.0       compiler_4.5.3         
#> [19] rlang_1.3.0             sass_0.4.10             tools_4.5.3            
#> [22] igraph_2.3.3            yaml_2.3.12             corrplot_0.95          
#> [25] knitr_1.52              S4Arrays_1.10.1         htmlwidgets_1.6.4      
#> [28] bit_4.6.0               DelayedArray_0.36.1     plyr_1.8.9             
#> [31] RColorBrewer_1.1-3      abind_1.4-8             withr_3.0.3            
#> [34] purrr_1.2.2             desc_1.4.3              grid_4.5.3             
#> [37] xtable_1.8-8            ggplot2_4.0.3           scales_1.4.0           
#> [40] MASS_7.3-65             cli_3.6.6               rmarkdown_2.32         
#> [43] crayon_1.5.3            ragg_1.5.2              otel_0.2.0             
#> [46] robustbase_0.99-7       httr_1.4.9              reshape2_1.4.5         
#> [49] BiocBaseUtils_1.12.0    DBI_1.3.0               cachem_1.1.0           
#> [52] stringr_1.6.0           splines_4.5.3           AnnotationDbi_1.72.0   
#> [55] AnnotationFilter_1.34.0 XVector_0.50.0          vctrs_0.7.3            
#> [58] Matrix_1.7-4            jsonlite_2.0.0          naniar_1.1.0           
#> [61] visdat_0.6.0            bit64_4.8.6             clue_0.3-68            
#> [64] systemfonts_1.3.2       tidyr_1.3.2             jquerylib_0.1.4        
#> [67] annotate_1.88.0         glue_1.8.1              DEoptimR_1.2-1         
#> [70] pkgdown_2.2.1           uniprotREST_1.0.0       stringi_1.8.9          
#> [73] gtable_0.3.6            tibble_3.3.1            pillar_1.11.1          
#> [76] htmltools_0.5.9         R6_2.6.1                textshaping_1.0.5      
#> [79] evaluate_1.0.5          lattice_0.22-9          backports_1.5.1        
#> [82] png_0.1-9               memoise_2.0.1           bslib_0.12.0           
#> [85] Rcpp_1.1.2              checkmate_2.3.4         SparseArray_1.10.10    
#> [88] xfun_0.60               MsCoreUtils_1.22.1      fs_2.1.0               
#> [91] pkgconfig_2.0.3
```
