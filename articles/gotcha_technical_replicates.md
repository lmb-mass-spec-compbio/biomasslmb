# Pitfall: technical replicates are not biological replicates

A biological replicate is a separate biological sample: a different
culture, a different animal, a different pulldown. A technical replicate
is the same sample measured more than once — the same digest injected
twice, most often.

They answer different questions. Biological replicates tell you how much
the thing you are studying varies; technical replicates tell you how
much your instrument varies. A statistical test comparing two conditions
is a statement about the first, and it needs to know how many biological
replicates there are. If you hand it technical replicates as though they
were biological ones, it will believe you.

This article shows what that costs, and what to do with technical
replicates instead.

## Load required packages

``` r

library(QFeatures)
library(biomasslmb)
library(dplyr)
library(limma)
```

## The data

`lfq_dda_techrep_peptides.txt.gz` is MaxQuant peptide-level output from
a bait pulldown against a tag-only control. Three biological replicates
of each, and every one of those six samples was acquired twice — twelve
runs in total.

The sample names carry that structure, and the rest of the article
relies on it: `Condition_Replicate_Techrep`, so `IP_2_1` and `IP_2_2`
are the two acquisitions of the same pulldown, while `IP_1_1` and
`IP_2_1` are different pulldowns. The design table below is built by
pulling those three fields back out of the column names.

``` r

pep_inf <- system.file("extdata", "lfq_dda_techrep_peptides.txt.gz",
                       package = "biomasslmb")
infdf <- read.delim(gzfile(pep_inf))

intensity_cols_ix <- grep('^(Control|IP)_', colnames(infdf))

exp_design <- data.frame(quantCols = colnames(infdf)[intensity_cols_ix]) %>%
  mutate(
    Condition = sub('_.*', '', quantCols),
    Replicate = as.integer(sub('^[A-Za-z]+_(\\d+)_\\d+$', '\\1', quantCols)),
    Tech_rep  = as.integer(sub('.*_(\\d+)$', '\\1', quantCols)))

knitr::kable(exp_design)
```

| quantCols   | Condition | Replicate | Tech_rep |
|:------------|:----------|----------:|---------:|
| Control_1_1 | Control   |         1 |        1 |
| Control_1_2 | Control   |         1 |        2 |
| Control_2_1 | Control   |         2 |        1 |
| Control_2_2 | Control   |         2 |        2 |
| Control_3_1 | Control   |         3 |        1 |
| Control_3_2 | Control   |         3 |        2 |
| IP_1_1      | IP        |         1 |        1 |
| IP_1_2      | IP        |         1 |        2 |
| IP_2_1      | IP        |         2 |        1 |
| IP_2_2      | IP        |         2 |        2 |
| IP_3_1      | IP        |         3 |        1 |
| IP_3_2      | IP        |         3 |        2 |

`proteotypic = TRUE` keeps only the peptides MaxQuant flagged as
matching a single protein, in its `Unique..Proteins.` column. MaxQuant
assigns every peptide one leading razor protein rather than reporting a
tie, so this is the only protein-inference filter
[`filter_features_mq_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_mq_dda.md)
has to offer.

``` r

qf <- readQFeatures(assayData = infdf, quantCols = intensity_cols_ix,
                    colData = exp_design, name = 'peptides_raw')
qf <- sync_coldata(qf, 'peptides_raw')
qf[['peptides_raw']] <- zeroIsNA(qf[['peptides_raw']])

qf[['peptides']] <- filter_features_mq_dda(
  qf[['peptides_raw']],
  contaminant_proteins = get_maxquant_cont_accessions(),
  filter_contaminant = TRUE, filter_associated_contaminant = TRUE,
  proteotypic = TRUE, remove_no_quant = TRUE)

qf[['peptides']] <- normalize(logTransform(qf[['peptides']], base = 2),
                              method = 'diff.median')
qf <- sync_coldata(qf, 'peptides')

nrow(qf[['peptides']])
#> [1] 2632
```

## How much do they actually agree?

Before anything else, it is worth confirming that these are technical
replicates in behaviour and not just in name. For every pair of samples
within a condition, take the standard deviation of the peptide-level
log2 differences: a small value means the two agree.

``` r

quant <- assay(qf[['peptides']])

pair_spread <- function(a, b) sd(quant[, a] - quant[, b], na.rm = TRUE)

# Technical pairs: same biological replicate r, its two acquisitions 1 and 2
technical <- unlist(lapply(c('Control', 'IP'), function(cond) {
  sapply(1:3, function(r) {
    pair_spread(paste0(cond, '_', r, '_1'), paste0(cond, '_', r, '_2'))
  })
}))

# Biological pairs: different replicates, holding the acquisition t fixed so
# that the two kinds of pair differ only in which factor is being varied
biological <- unlist(lapply(c('Control', 'IP'), function(cond) {
  lapply(1:2, function(t) {
    combn(1:3, 2, function(p) {
      pair_spread(paste0(cond, '_', p[1], '_', t),
                  paste0(cond, '_', p[2], '_', t))
    })
  })
}))

c(technical_pairs = round(median(technical), 3),
  biological_pairs = round(median(biological), 3))
#>  technical_pairs biological_pairs 
#>            0.410            0.848
```

Two acquisitions of the same digest disagree about half as much as two
separate biological samples of the same condition. That is the expected
picture, and it is exactly why the two cannot be interchanged: a
substantial part of the variation between biological replicates is real
biological variation, and no amount of re-injecting will reduce it.

## Summarising to protein

``` r

qf <- aggregateFeatures(qf, i = 'peptides',
                        fcol = 'Leading.razor.protein', name = 'protein',
                        fun = MsCoreUtils::robustSummary, maxit = 10000)
qf <- sync_coldata(qf, 'protein')

protein <- assay(qf[['protein']])
sample_data <- colData(qf[['protein']])

nrow(protein)
#> [1] 393
```

``` r

test_condition <- function(quant, condition) {
  group <- factor(condition, levels = c('Control', 'IP'))
  n_quant <- sapply(levels(group), function(l) {
    rowSums(!is.na(quant[, group == l, drop = FALSE]))
  })
  testable <- apply(n_quant, 1, min) >= 2

  fit <- eBayes(lmFit(quant[testable, ], model.matrix(~ group)))
  tt <- topTable(fit, coef = 2, number = Inf)
  c(tested = sum(testable), significant = sum(tt$adj.P.Val < 0.05))
}
```

## The error: twelve samples that are really six

The tempting thing to do with twelve columns is to test them as twelve
samples.

``` r

test_condition(protein, sample_data$Condition)
#>      tested significant 
#>         296          17
```

Now the correct analysis. The two acquisitions of a biological sample
are averaged first, giving one value per biological replicate, and the
test sees three against three.

``` r

biological_sample <- paste(sample_data$Condition, sample_data$Replicate)

averaged <- sapply(unique(biological_sample), function(s) {
  rowMeans(protein[, biological_sample == s, drop = FALSE], na.rm = TRUE)
})
averaged[is.nan(averaged)] <- NA

test_condition(averaged, sub(' .*', '', colnames(averaged)))
#>      tested significant 
#>         282           7
```

Roughly twice as many proteins reach significance when the technical
replicates are treated as independent.

The extra hits are not evidence of anything. `limma` estimates each
protein’s variance from the scatter of samples within a group, and the
two acquisitions of one digest scatter half as much as two separate
samples would. Including both therefore does two things at once: it
halves the standard error by pretending there are six independent
observations per group, and it pulls the variance estimate down towards
the technical noise rather than the biological noise. Both push p-values
in the same direction, and neither reflects any additional biological
information — the experiment still contains three biological replicates
per condition, whatever the file says.

This is pseudo-replication, and proteomics makes it unusually easy to
commit, because a technical replicate arrives as just another column and
nothing in the file marks it as one.

## Averaging beats picking one

Given that the two acquisitions have to be reduced to one value, an
alternative to averaging is to keep one and discard the other. It is
worth seeing what that costs.

``` r

first_only <- protein[, sample_data$Tech_rep == 1]
second_only <- protein[, sample_data$Tech_rep == 2]

rbind(
  averaged = test_condition(averaged, sub(' .*', '', colnames(averaged))),
  technical_replicate_1 = test_condition(
    first_only, sample_data$Condition[sample_data$Tech_rep == 1]),
  technical_replicate_2 = test_condition(
    second_only, sample_data$Condition[sample_data$Tech_rep == 2]))
#>                       tested significant
#> averaged                 282           7
#> technical_replicate_1    255           3
#> technical_replicate_2    243           0
```

Averaging wins on both counts, and by a wide margin. It tests more
proteins, because a protein missed in one acquisition can be rescued by
the other, and it finds more of them, because averaging two measurements
of the same digest reduces the technical component of the noise without
touching the biological component.

Note also how much the two single-acquisition analyses differ from each
other. Neither is wrong; they are two draws from the same technical
noise, and the gap between them is a fair warning about how much a
result can move on measurement noise alone when the number of biological
replicates is small.

## Modelling the pairing instead

Averaging is not the only correct treatment. `limma` will take all
twelve columns as they are, provided it is told which of them belong
together:
[`duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
estimates a single consensus correlation between observations that share
a block, and [`lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) then
uses it to discount the fact that two acquisitions of one digest are not
independent measurements.

Which proteins can be tested has to be decided on biological samples
rather than on columns, so that is worked out first.

``` r

group <- factor(sample_data$Condition, levels = c('Control', 'IP'))

n_bio <- sapply(levels(group), function(l) {
  cols <- which(group == l)
  apply(protein[, cols, drop = FALSE], 1, function(x) {
    length(unique(biological_sample[cols][!is.na(x)]))
  })
})

testable_bio <- apply(n_bio, 1, min) >= 2

# the same rule applied to columns rather than to biological samples
testable_acq <- apply(sapply(levels(group), function(l) {
  rowSums(!is.na(protein[, group == l, drop = FALSE]))
}), 1, min) >= 2

c(by_biological_sample = sum(testable_bio), by_acquisition = sum(testable_acq))
#> by_biological_sample       by_acquisition 
#>                  282                  296

design <- model.matrix(~ group)

corfit <- duplicateCorrelation(protein[testable_bio, ], design,
                               block = biological_sample)

round(corfit$consensus.correlation, 3)
#> [1] 0.637
```

That correlation is the pairing made explicit, and it is the number the
whole approach turns on: it says how much of a protein’s value is shared
between two acquisitions of the same digest. A value near zero would
mean the repeat acquisitions carry independent information and could be
treated as separate samples; a value near one would mean the second
acquisition adds nothing at all.

``` r

fit <- eBayes(lmFit(protein[testable_bio, ], design,
                    block = biological_sample,
                    correlation = corfit$consensus.correlation))

blocked <- topTable(fit, coef = 2, number = Inf)

rbind(
  averaged = test_condition(averaged, sub(' .*', '', colnames(averaged))),
  blocked = c(tested = sum(testable_bio),
              significant = sum(blocked$adj.P.Val < 0.05)),
  twelve_samples = test_condition(protein, sample_data$Condition))
#>                tested significant
#> averaged          282           7
#> blocked           282           8
#> twelve_samples    296          17
```

**Blocking lands where averaging lands, not where the naive analysis
lands.** On the same set of proteins it finds one more than averaging
and less than half what treating the acquisitions as independent samples
produced. That is the result to take from this section: the two
defensible routes agree with each other, and the gap between them and
the twelve-sample analysis is the pseudo-replication rather than a
difference of method.

Which to use is mostly a matter of what else the design contains.
Averaging is simpler, needs no extra assumption, and leaves an ordinary
matrix that any downstream tool will accept. Blocking keeps the
individual acquisitions visible, which matters when the repeats are
unbalanced — three injections of one sample and one of another — since
averaging silently gives the better-measured sample the same weight as
the other. It also assumes one correlation describes every block, which
is a real assumption and is worth checking against the per-pair spreads
calculated earlier when block sizes differ.

The testability rule deserves the attention it was given above. Counting
quantified acquisitions instead of quantified biological samples would
have admitted 14 further proteins here, each one passing on two values
in a condition that both came from the same digest. The block structure
discounts them correctly in the variance, but no amount of correct
discounting turns one biological sample into two, and a fold change
resting on one is not something to report as though it rested on three.

## What technical replicates are for

None of this makes them a waste of instrument time. They are the only
thing in an experiment that separates measurement variability from
biological variability, which makes them useful for questions the
biological replicates cannot answer:

- **Is a surprising result real or an acquisition artefact?** If a
  protein’s two technical replicates disagree far more than the typical
  pair, its value in that sample is unreliable regardless of what the
  statistics say.
- **Is one sample’s acquisition bad?** The `pair_spread` calculation
  above, per sample rather than pooled, identifies a run that went wrong
  — a poor injection, a fouled column — in a way that a single
  acquisition per sample cannot.
- **Is the instrument the limiting factor?** If technical variability
  approaches biological variability, more biological replicates will not
  help as much as fixing the acquisition.

The rule is simply that they belong in quality control, not in the
degrees of freedom of the test.

## Summary

- A technical replicate measures the instrument; a biological replicate
  measures the biology. Here two acquisitions of one digest disagreed
  about half as much as two biological samples.
- **Testing technical replicates as though they were biological
  replicates roughly doubled the number of significant proteins here.**
  The extra hits come from a halved standard error and a variance
  estimate pulled towards technical noise, not from any additional
  information.
- Average technical replicates to one value per biological sample before
  testing. Do not pick one and discard the rest: averaging tested more
  proteins and found more of them than either single acquisition.
- Blocking on biological sample with
  [`duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
  is the alternative to averaging and agreed with it closely here. It
  suits unbalanced repeats, at the cost of assuming one correlation
  describes every block.
- Decide what is testable from the number of quantified biological
  samples, not the number of quantified columns, whichever route is
  taken.
- Use them for quality control — spotting unreliable measurements, bad
  runs, and whether the instrument or the biology is the limiting source
  of variation.

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
#>  [1] limma_3.66.0                dplyr_1.2.1                
#>  [3] biomasslmb_0.1.0            QFeatures_1.20.0           
#>  [5] MultiAssayExperiment_1.36.2 SummarizedExperiment_1.40.0
#>  [7] Biobase_2.70.0              GenomicRanges_1.62.1       
#>  [9] Seqinfo_1.0.0               IRanges_2.44.0             
#> [11] S4Vectors_0.48.1            BiocGenerics_0.56.0        
#> [13] generics_0.1.4              MatrixGenerics_1.22.0      
#> [15] matrixStats_1.5.0          
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1        farver_2.1.2            blob_1.3.0             
#>  [4] Biostrings_2.78.0       S7_0.2.2                fastmap_1.2.0          
#>  [7] lazyeval_0.2.3          XML_3.99-0.24           digest_0.6.39          
#> [10] lifecycle_1.0.5         cluster_2.1.8.2         ProtGenerics_1.42.0    
#> [13] statmod_1.5.2           survival_3.8-6          KEGGREST_1.50.0        
#> [16] RSQLite_3.53.3          magrittr_2.0.5          genefilter_1.92.0      
#> [19] compiler_4.5.3          rlang_1.3.0             sass_0.4.10            
#> [22] tools_4.5.3             igraph_2.3.3            yaml_2.3.12            
#> [25] corrplot_0.95           knitr_1.52              S4Arrays_1.10.1        
#> [28] htmlwidgets_1.6.4       bit_4.6.0               DelayedArray_0.36.1    
#> [31] plyr_1.8.9              RColorBrewer_1.1-3      abind_1.4-8            
#> [34] withr_3.0.3             purrr_1.2.2             desc_1.4.3             
#> [37] grid_4.5.3              xtable_1.8-8            ggplot2_4.0.3          
#> [40] scales_1.4.0            MASS_7.3-65             cli_3.6.6              
#> [43] rmarkdown_2.32          crayon_1.5.3            ragg_1.5.2             
#> [46] otel_0.2.0              robustbase_0.99-7       httr_1.4.9             
#> [49] reshape2_1.4.5          BiocBaseUtils_1.12.0    DBI_1.3.0              
#> [52] cachem_1.1.0            stringr_1.6.0           splines_4.5.3          
#> [55] AnnotationDbi_1.72.0    AnnotationFilter_1.34.0 XVector_0.50.0         
#> [58] vctrs_0.7.3             Matrix_1.7-4            jsonlite_2.0.0         
#> [61] naniar_1.1.0            visdat_0.6.0            bit64_4.8.6            
#> [64] clue_0.3-68             systemfonts_1.3.2       tidyr_1.3.2            
#> [67] jquerylib_0.1.4         annotate_1.88.0         glue_1.8.1             
#> [70] DEoptimR_1.2-1          pkgdown_2.2.1           uniprotREST_1.0.0      
#> [73] stringi_1.8.9           gtable_0.3.6            tibble_3.3.1           
#> [76] pillar_1.11.1           htmltools_0.5.9         R6_2.6.1               
#> [79] textshaping_1.0.5       evaluate_1.0.5          lattice_0.22-9         
#> [82] backports_1.5.1         png_0.1-9               memoise_2.0.1          
#> [85] bslib_0.12.0            Rcpp_1.1.2              checkmate_2.3.4        
#> [88] SparseArray_1.10.10     xfun_0.60               MsCoreUtils_1.22.1     
#> [91] fs_2.1.0                pkgconfig_2.0.3
```
