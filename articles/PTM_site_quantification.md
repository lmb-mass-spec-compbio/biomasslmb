# PTM designs: site-level quantification and normalisation

A post-translational modification (PTM) experiment asks a different
question from a whole-proteome one. The quantity of interest is not how
much of a protein there is, but how much of it carries a modification at
a particular residue. Everything below follows from that difference.

Because modified peptides are a small fraction of the peptides in a
digest, they have to be enriched before they can be quantified with any
depth. A common design is to TMT label every sample once, pool the
labelled material, and then split the pool: one part is enriched for the
modification, the other is left alone. The two parts are acquired as
separate runs, giving two datasets over the same samples — an
**enriched** one that reports modified peptides, and a **total** one
that reports protein abundance.

The alternative is to enrich each sample separately and quantify
label-free, and the reason not to is that the enrichment is the most
variable step in the workflow. Recovery from an immobilised metal
affinity or titanium dioxide column depends on binding capacity, wash
stringency and elution, and none of those repeat exactly between
preparations. Enriching separately imprints that variability on the
samples individually, where nothing in the data distinguishes it from
biology. Labelling first and pooling moves the enrichment downstream of
the point at which the samples become one tube, so a single enrichment
is applied to all of them and its idiosyncrasies land on every channel
alike.

Two further things follow from pooling. Modified peptides sit low in the
dynamic range and are identified less reproducibly than unmodified ones,
so a label-free DDA experiment loses many sites to missing values in a
subset of samples; within a plex, a site identified once is quantified
in every channel from the same spectrum. And the total fraction is split
from the pool after labelling, which is what later makes it usable as a
normalisation reference: it is the same labelled material, so it carries
the same loading differences as the enriched fraction.

The costs are the usual isobaric ones. The plex sets a ceiling on how
many samples can share an enrichment, and a design that exceeds it needs
bridge channels and the treatment in [multi-plex
TMT](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/TMT_multiplex.md).
Co-isolation compresses fold changes, which SPS MS3 (McAlister et al.
2014) and the filtering below limit, but cannot completely resolve.
Label-free is a reasonable choice where the enrichment can be tightly
controlled, and DIA acquisition answers the missing value argument, but
neither recovers the shared enrichment.

With TMT-label, pool then enrich experiment design, there are three
important considerations which we cover here:

- **Which residue is modified matters**, so the search engine’s
  localisation confidence has to be assessed and the site placed within
  the protein, not just within the peptide.
- **The enriched fraction cannot be normalised against itself.** A
  treatment that changes phosphorylation globally is the signal, and
  median-centring the enriched data would remove it. The total fraction
  supplies the reference instead.
- **A site can change because its protein changed.** Dividing the site
  out by its protein separates the two, and the answer is not always the
  same as the site abundance.

This vignette works through a phosphoproteomics experiment searched with
Proteome Discoverer (PD). The principles apply to any enriched PTM; only
the search-engine columns differ. MaxQuant’s equivalents are
[`filter_maxquant_ptm()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_maxquant_ptm.md)
and
[`add_filter_ptm_pos_rowdata_mq()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_filter_ptm_pos_rowdata_mq.md),
which take the same role as
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
below, and [peptides that cannot be
localised](#peptides-that-cannot-be-localised) works the same data
through both engines side by side.

Everything up to the localisation step is the routine PSM processing
covered in [TMT QC and
summarisation](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/TMT_PSM_QC_Summarisation.md),
so it is run here with less commentary.

## Load required packages

``` r

library(QFeatures)
library(biomasslmb)
library(ggplot2)
library(dplyr)
library(limma)
```

## The experimental design

`psm_tmt_phospho`, `psm_tmt_phospho_total` and `tmt_phospho_design` are
datasets available from the `biomasslmb` package, derived from a real
TMTpro experiment in mouse fibroblasts. A drug treatment is compared
against a vehicle control at two timepoints, with four replicates of
each combination, all labelled in a single plex. `psm_tmt_phospho` is
the PD PSM-level output for the phospho-enriched fraction and
`psm_tmt_phospho_total` for the matched total fraction, both truncated
to a subset of proteins for a manageable vignette.

``` r

knitr::kable(tmt_phospho_design)
```

|              | Condition | Timepoint | Replicate | quantCols    |
|:-------------|:----------|:----------|:----------|:-------------|
| Control_T1_1 | Control   | T1        | 1         | Control_T1_1 |
| Control_T1_2 | Control   | T1        | 2         | Control_T1_2 |
| Control_T1_3 | Control   | T1        | 3         | Control_T1_3 |
| Control_T1_4 | Control   | T1        | 4         | Control_T1_4 |
| Control_T2_1 | Control   | T2        | 1         | Control_T2_1 |
| Control_T2_2 | Control   | T2        | 2         | Control_T2_2 |
| Control_T2_3 | Control   | T2        | 3         | Control_T2_3 |
| Control_T2_4 | Control   | T2        | 4         | Control_T2_4 |
| Treated_T1_1 | Treated   | T1        | 1         | Treated_T1_1 |
| Treated_T1_2 | Treated   | T1        | 2         | Treated_T1_2 |
| Treated_T1_3 | Treated   | T1        | 3         | Treated_T1_3 |
| Treated_T1_4 | Treated   | T1        | 4         | Treated_T1_4 |
| Treated_T2_1 | Treated   | T2        | 1         | Treated_T2_1 |
| Treated_T2_2 | Treated   | T2        | 2         | Treated_T2_2 |
| Treated_T2_3 | Treated   | T2        | 3         | Treated_T2_3 |
| Treated_T2_4 | Treated   | T2        | 4         | Treated_T2_4 |

The two fractions share this design, because they are two acquisitions
of the same labelled pool. The TMT tag that carried a given sample is
the same in both.

## Defining the contaminant proteins

As for any experiment, we need the contaminant accessions to filter
against. This search used the ‘0602_Universal Contaminants’ database
(Frankenfield et al. 2022).

``` r

contaminant_fasta_inf <- system.file(
  "extdata", "0602_Universal_Contaminants.fasta.gz",
  package = "biomasslmb"
)

contaminant_accessions <- get_contaminant_fasta_accessions(contaminant_fasta_inf)
contaminant_accessions <- c(contaminant_accessions,
                            sub('^Cont_', '', contaminant_accessions))
```

## Read in and filter both fractions

The two fractions are read into separate `QFeatures` objects. They
cannot share one, because they hold different features: modified
peptides in one, all peptides in the other. Keeping them apart also
keeps it obvious which assay is which when they are combined at the end.

The filtering applied here is the same for both, and follows the
PSM-level processing in the [TMT PSM QC and
summarisation](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/TMT_PSM_QC_Summarisation.md)
vignette: drop contaminants and PSMs without a unique master protein,
then drop PSMs with high co-isolation, low signal:noise, or a search
engine rank below 1.

``` r

read_and_filter <- function(infdf) {
  qf <- readQFeatures(assayData = infdf,
                      colData = tmt_phospho_design,
                      name = 'psm_raw')
  qf <- sync_coldata(qf, 'psm_raw')

  # A more accurate average S:N ratio value than PD reports
  qf[['psm_raw']] <- update_average_sn(qf[['psm_raw']])

  qf[['psm_filtered']] <- filter_features_pd_dda(
    qf[['psm_raw']],
    contaminant_proteins = contaminant_accessions,
    filter_contaminant = TRUE,
    filter_associated_contaminant = TRUE,
    unique_master = TRUE)

  qf[['psm_quality']] <- filter_TMT_PSMs(
    qf[['psm_filtered']], inter_thresh = 50, sn_thresh = 5)

  filterFeatures(qf, ~ Rank == 1, i = 'psm_quality')
}

phospho <- read_and_filter(psm_tmt_phospho)
total <- read_and_filter(psm_tmt_phospho_total)
```

``` r

data.frame(
  fraction = c('Phospho-enriched', 'Total'),
  raw = c(nrow(phospho[['psm_raw']]), nrow(total[['psm_raw']])),
  filtered = c(nrow(phospho[['psm_quality']]), nrow(total[['psm_quality']])))
#>           fraction   raw filtered
#> 1 Phospho-enriched 13354     6757
#> 2            Total 17045     9119
```

## Localising the modification

A search engine reports which peptide was identified and which residues
carry the modification, but the second of those is much less certain
than the first. A peptide with several serines and one phosphate can
often be explained about equally well by a phosphate on any of them: the
fragment ions that would distinguish the possibilities may not have been
observed. PD’s ptmRS node scores each candidate site, and
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
reads those scores out of the `ptmRS.Best.Site.Probabilities` column.

`threshold` sets the minimum score a site needs. If any site on a
peptide falls below it, the modifications on that peptide are all
discarded rather than partially kept — a peptide whose phosphate could
be on either of two residues cannot be assigned to either. It is not
silent about the region, though, and [peptides that cannot be
localised](#peptides-that-cannot-be-localised) below recovers what it
does say.

``` r

phospho[['psm_localised']] <- parse_PTM_scores_pd(
  phospho[['psm_quality']], threshold = 75)
#> Removed 0 Features where the ptm_col value == `Inconclusive data`
#> Total Features: 6757
#> Total detected PTMFeatures: 6554
#> Features passing filter: 5058
#> Features failing filter: 1496
#> BiPTM/multiPTM Features where some sites fail filter: 349
#> Total detected sites: 10398
#> Sites passing filter: 6701
#> Sites failing filter: 3697
#> monoPTM passing filter: 3619
#> biPTM passing filter: 1237
#> multiPTM passing filter: 202
#> Too many isoforms: 0
```

The log distinguishes features from sites, which matters when peptides
carry more than one modification.
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
annotates rather than filters, writing empty values where nothing
passed, so the failures still need removing.

``` r

phospho[['psm_localised']] <- phospho[['psm_localised']][
  grepl('Phospho', rowData(phospho[['psm_localised']])$ptms), ]

nrow(phospho[['psm_localised']])
#> [1] 5058
```

A threshold of 75 is a common choice and is the one used here; 95 is
stricter and is worth considering when the question rests on a single
site. The cost is one-sided — raising it discards peptides, it does not
add any — so the appropriate value depends on whether you would rather
lose sites or believe some that are misplaced.

## Placing sites within the protein

The positions reported so far are positions within the peptide. Two PSMs
of the same site reached by different missed cleavages will disagree
about it, and no two proteins can be compared on it.
[`add_ptm_positions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ptm_positions.md)
converts them to positions within the protein by re-digesting the
protein sequences and locating each peptide.

That needs the sequences. In a real analysis you would fetch them for
your master proteins with
[`make_fasta()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/make_fasta.md);
a FASTA covering this dataset’s proteins is included with the package so
that building the vignette does not depend on UniProt being reachable.

``` r

proteome_fasta <- tempfile(fileext = '.fasta')
make_fasta(unique(rowData(phospho[['psm_localised']])$Master.Protein.Accessions),
           file = proteome_fasta)
```

``` r

proteome_fasta <- system.file(
  "extdata", "tmt_phospho_proteome.fasta.gz", package = "biomasslmb")

phospho[['psm_localised']] <- add_ptm_positions(
  phospho[['psm_localised']],
  proteome_fasta = proteome_fasta,
  master_protein_col = 'Master.Protein.Accessions')
```

This adds `start` and `end` for the peptide, `ptm_positions_prot` for
the sites within the protein, and `ptm_name`, which combines residue and
position into the conventional form.

``` r

rowData(phospho[['psm_localised']]) %>%
  data.frame() %>%
  select(Sequence, Master.Protein.Accessions, ptm_positions,
         start, ptm_positions_prot, ptm_name) %>%
  head(4)
#>       Sequence Master.Protein.Accessions ptm_positions start ptm_positions_prot
#> 2 SHSSSSSPSPSR                    Q8BTI8         9; 11   903           911; 913
#> 4      GDSDDGR                    Q9CYI0             3    53                 55
#> 5    SASGSSSDR                    Q8BTI8          6; 7  2031         2036; 2037
#> 7     SDGAGGAR                    Q6ZQ58             1   302                302
#>       ptm_name
#> 2   S911; S913
#> 4          S55
#> 5 S2036; S2037
#> 7         S302
```

Two things stop a peptide from being placed. Its master protein may have
no sequence in the FASTA, which happens when an accession has been
withdrawn or merged since the search; and the peptide may occur at more
than one position in its protein, in which case `start` holds several
values and the site cannot be assigned to one of them. Neither is
common, and in a subset this size either count can come out at zero, but
both have to be checked for and removed: a site with no unambiguous
position cannot be summarised or compared.

``` r

site_positions <- rowData(phospho[['psm_localised']])

data.frame(
  reason = c('No sequence for the master protein',
             'Peptide occurs more than once in its protein',
             'Placed'),
  PSMs = c(sum(is.na(site_positions$start)),
           sum(grepl(';', site_positions$start)),
           sum(!is.na(site_positions$start) & !grepl(';', site_positions$start))))
#>                                         reason PSMs
#> 1           No sequence for the master protein    2
#> 2 Peptide occurs more than once in its protein    0
#> 3                                       Placed 5056
```

``` r

phospho[['psm_sites']] <- phospho[['psm_localised']][
  !is.na(site_positions$start) & !grepl(';', site_positions$start), ]
```

## Peptides that cannot be localised

The localisation filter is the only step in this workflow that discards
data on the strength of a probability rather than a measurement, and it
is the most expensive one. Applied to the quality-filtered PSMs of this
dataset at a threshold of 75, it removes a little under a quarter of
them.

``` r

psm_all_ptms <- phospho[['psm_quality']]
candidates <- parse_ptm_candidates_pd(psm_all_ptms)

resolution <- candidates %>%
  group_by(row, n_ptms) %>%
  summarise(n_localised = sum(prob >= 0.75), .groups = 'drop') %>%
  mutate(resolved = n_localised == n_ptms)

c(PSMs = nrow(resolution),
  discarded = sum(!resolution$resolved),
  percent = round(100 * mean(!resolution$resolved), 1))
#>      PSMs discarded   percent 
#>    6554.0    1498.0      22.9
```

**A peptide whose phosphate is tied between two serines is not an
uninformative peptide.** It establishes that the region is
phosphorylated and by how much, and leaves open only which of the two
residues carries the phosphate. The filter throws away the
quantification along with the ambiguity, and what it throws away is not
a random subset: serine-rich and threonine-rich regions are exactly the
ones a spectrum struggles to resolve, so the sites that go missing are
concentrated in the sequences most worth looking at.

The alternative is to keep those peptides by pooling them. Peptides
whose sites all reach the threshold are untouched and keep their own
single-residue identity. Peptides that do not resolve contribute every
candidate residue that is still plausible, and candidate sets that
overlap are merged into one group, which is quantified as a unit. The
cost is stated plainly at the end of this section: a pooled group names
a region rather than a residue, and several things you would want to do
with a site cannot be done with one.

### Reading every candidate residue

[`parse_ptm_candidates_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_pd.md)
reads `ptmRS.Phospho.Site.Probabilities`, **not** the
`ptmRS.Best.Site.Probabilities` column that
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
uses. The distinction is easy to miss and fails silently. Best Site
reports only the winning isoform’s sites, so the alternatives a group
would be built from are absent from it, and pointing the parser at that
column yields groups containing only the sites ptmRS already preferred.

``` r

example_row <- resolution$row[!resolution$resolved][1]

rowData(psm_all_ptms)[example_row, c('ptmRS.Best.Site.Probabilities',
                                     'ptmRS.Phospho.Site.Probabilities')] %>%
  data.frame() %>%
  t()
#>                                  17                                                                           
#> ptmRS.Best.Site.Probabilities    "S3(Phospho): 46.88; S5(Phospho): 46.88; S9(Phospho): 100; S11(Phospho): 100"
#> ptmRS.Phospho.Site.Probabilities "S(1): 6.2; S(3): 46.9; S(5): 46.9; S(9): 100.0; S(11): 100.0"
```

The parser returns one row per candidate residue, with the position
within the peptide, the residue, the probability and the number of
modifications on the peptide.
[`parse_ptm_candidates_mq()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_mq.md)
returns the same table from MaxQuant’s `Phospho..STY..Probabilities`
column, and everything after this point is shared between the two.

``` r

head(candidates, 4)
#>   row pep_pos residue prob n_ptms
#> 1   1       1       S    0      2
#> 2   1       3       S    0      2
#> 3   1       4       S    0      2
#> 4   1       5       S    0      2
```

Probabilities are on a 0 to 1 scale in both parsers, so a threshold
means the same thing whichever search engine produced the data. ptmRS
reports percentages and
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
takes its `threshold` as one — 75 above, 0.75 below — which is the one
place in this workflow where the same quantity is expressed two ways.

### The probabilities are a budget, not a score

Each modification has to sit on some residue, so the probabilities
reported for a peptide sum to the number of modifications on it. That
holds in this dataset to within the rounding ptmRS applies.

``` r

candidates %>%
  group_by(row, n_ptms) %>%
  summarise(total = sum(prob), .groups = 'drop') %>%
  summarise(min = round(min(total / n_ptms), 3),
            max = round(max(total / n_ptms), 3))
#> # A tibble: 1 × 2
#>     min   max
#>   <dbl> <dbl>
#> 1 0.997  1.00
```

**This is what makes the pooling quantitative rather than cosmetic.** A
candidate’s probability is its share of one modification, so the summed
probability of a set of candidates is the chance that the set contains
the true site, and the probability left outside the set is the chance
that it does not. A group is therefore not a vague claim about a region;
it is a statement with a coverage figure attached.

### Grouping the candidates

Grouping happens in protein coordinates, so the peptides need placing
first, as in the previous section — this time before the localisation
filter has been applied rather than after.

``` r

psm_all_ptms <- add_peptide_positions_from_cleavage(
  psm_all_ptms, proteome_fasta, master_protein_col = 'Master.Protein.Accessions')

psm_all_ptms <- add_ambiguous_ptm_group_rowdata(
  psm_all_ptms, candidates,
  master_protein_col = 'Master.Protein.Accessions',
  min_prob = 0.75)
```

Candidate residues are nodes, and two residues are joined when they are
candidates on the same peptide. The groups are the connected components
of that graph, so a “S32 or S33” peptide and a “S33 or S37” peptide
share a node and become one group over all three residues. The graph is
built separately for each protein and each modification count, which
stops a singly phosphorylated peptide from merging with a doubly
phosphorylated one covering the same residues: they overlap in sequence
but are different molecular species.

``` r

grouped <- data.frame(rowData(psm_all_ptms))

c(resolved_PSMs = sum(grouped$ptm_group_resolved, na.rm = TRUE),
  pooled_PSMs = sum(!grouped$ptm_group_resolved %in% TRUE &
                      !is.na(grouped$ptm_group_id)),
  unassigned = sum(is.na(grouped$ptm_group_id)))
#> resolved_PSMs   pooled_PSMs    unassigned 
#>          5054          1498           205
```

The PSMs the filter would have kept are unchanged — each still reports
its own residue — and the PSMs it would have discarded are now
distributed over a much smaller number of quantifiable groups.

``` r

grouped %>%
  filter(!is.na(ptm_group_id)) %>%
  group_by(pooled = !ptm_group_resolved) %>%
  summarise(PSMs = n(), features = n_distinct(ptm_group_id))
#> # A tibble: 2 × 3
#>   pooled  PSMs features
#>   <lgl>  <int>    <int>
#> 1 FALSE   5054     1351
#> 2 TRUE    1498      535
```

A pooled group’s size is set by how much genuine ambiguity there is, not
by how many PSMs fell into it.

``` r

grouped %>%
  filter(!ptm_group_resolved %in% TRUE, !is.na(ptm_group_id)) %>%
  distinct(ptm_group_id, ptm_group_n_candidates) %>%
  pull(ptm_group_n_candidates) %>%
  table()
#> .
#>   2   3   4   5   6   7   8   9  10  11 
#> 142 148 117  51  31  17  15   9   3   2
```

### The same grouping from MaxQuant output

Everything from
[`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)
onwards is the same for either search engine, because the parser has
already reduced the output to one candidate table. Only the step that
produces that table differs, and `psm_tmt_phospho_mq` — the
phospho-enriched fraction of the TMT18plex whose total fraction is
`psm_tmt_factorial` — is the MaxQuant equivalent of the data used above.

The PSM filtering is the MaxQuant workflow’s, covered in [TMT QC and
summarisation](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/TMT_PSM_QC_Summarisation.md):
reporter intensities of exactly zero become `NA`, contaminant accessions
come from
[`get_maxquant_cont_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_maxquant_cont_accessions.md),
and the decoy hits MaxQuant leaves in the export are removed.

``` r

mq_qf <- QFeatures::readQFeatures(assayData = psm_tmt_phospho_mq,
                               colData = tmt_factorial_design,
                               quantCols = rownames(tmt_factorial_design),
                               name = 'psm_raw')
mq_qf <- sync_coldata(mq_qf, 'psm_raw')
mq_qf[['psm_raw']] <- QFeatures::zeroIsNA(mq_qf[['psm_raw']])

mq_qf[['psm_quality']] <- filter_features_mq_dda(
  mq_qf[['psm_raw']],
  contaminant_proteins = get_maxquant_cont_accessions(),
  filter_contaminant = TRUE,
  filter_associated_contaminant = TRUE)

mq_psms <- mq_qf[['psm_quality']][
  rowData(mq_qf[['psm_quality']])$Phospho..STY..Probabilities != '', ]
```

MaxQuant writes the probabilities inline in the peptide sequence rather
than as a separate list, so a residue’s position is implied by how much
sequence precedes its value.
[`parse_ptm_candidates_mq()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_mq.md)
reads that format and returns the identical five columns.

``` r

rowData(mq_psms)$Phospho..STY..Probabilities[1]
#> [1] "AAAAAAS(0.003)AGS(0.144)S(0.814)AS(0.039)SGNQPPQELGLGELLEEFSR"

mq_candidates <- parse_ptm_candidates_mq(mq_psms)
head(mq_candidates, 3)
#>   row pep_pos residue  prob n_ptms
#> 1   1       7       S 0.003      1
#> 2   1      10       S 0.144      1
#> 3   1      11       S 0.814      1
```

The probabilities are already on a 0 to 1 scale here, and they obey the
same budget — MaxQuant’s rounding is finer than ptmRS’s, so the totals
sit closer to one.

``` r

mq_proteome_fasta <- system.file(
  "extdata", "tmt_phospho_mq_proteome.fasta.gz", package = "biomasslmb")

mq_psms <- add_peptide_positions_from_cleavage(mq_psms, mq_proteome_fasta)

mq_psms <- add_ambiguous_ptm_group_rowdata(mq_psms, mq_candidates, min_prob = 0.75)
```

That call is identical to the Proteome Discoverer one except for
`master_protein_col`, which defaults to MaxQuant’s
`Leading.razor.protein`. Supporting a third search engine would mean
writing a third parser and changing nothing else.

**The two engines disagree about how much is ambiguous, not about what
to do with it.** Compared at the same threshold, ptmRS leaves a larger
share of PSMs unlocalised and spreads each one over roughly twice as
many candidate residues, so it produces fewer but not smaller groups
from a similar number of PSMs.

``` r

engine_summary <- function(candidates, engine, min_prob = 0.75) {
  unresolved <- candidates %>%
    group_by(row, n_ptms) %>%
    summarise(resolved = sum(prob >= min_prob) == n_ptms[1], .groups = 'drop')

  data.frame(
    engine = engine,
    PSMs = nrow(unresolved),
    pc_unresolved = round(100 * mean(!unresolved$resolved), 1),
    candidates_per_PSM = round(nrow(candidates) / n_distinct(candidates$row), 2),
    pc_candidates_at_zero = round(100 * mean(candidates$prob == 0), 1),
    lowest_prob = min(candidates$prob))
}

rbind(engine_summary(mq_candidates, 'MaxQuant'),
      engine_summary(candidates, 'ptmRS'))
#>     engine PSMs pc_unresolved candidates_per_PSM pc_candidates_at_zero
#> 1 MaxQuant 6790          18.3               2.14                   0.0
#> 2    ptmRS 6554          22.9               4.16                  34.7
#>   lowest_prob
#> 1       0.001
#> 2       0.000
```

The last two columns are the practical difference, and they change what
`min_candidate_prob` does. MaxQuant never writes a probability below
0.001, so its lowest values are a reporting floor and a threshold of
0.02 sits an order of magnitude above it. ptmRS writes exact zeros and
writes a great many of them — a third of all its candidates here — so on
Proteome Discoverer data any threshold above zero clears that block out
before the particular value chosen starts to matter.

A threshold of 0.75 is used above so that the two engines are compared
on equal terms. MaxQuant pipelines commonly use 0.501 instead, which is
the `min_prob` default, and on this dataset it halves the share of PSMs
that fail to localise, from 18.3% to 9.2%.

### Choosing min_candidate_prob

`min_candidate_prob` is the probability below which a residue on an
unresolved peptide is treated as ruled out rather than kept as a
candidate. Because the probabilities are a budget, it is a coverage
guarantee: the mass it discards is the chance the resulting group
excludes the real site.
[`summarise_ptm_groups()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/summarise_ptm_groups.md)
reports both sides of that trade at a given value.

``` r

lapply(c(0.01, 0.02, 0.05, 0.1, 0.2), function(threshold) {
  summarise_ptm_groups(psm_all_ptms, candidates, threshold, min_prob = 0.75,
                       master_protein_col = 'Master.Protein.Accessions')
}) %>%
  bind_rows() %>%
  select(min_candidate_prob, pooled_groups, median_size, pc_size_5_plus,
         p95_span, pc_miss_over_5)
#>   min_candidate_prob pooled_groups median_size pc_size_5_plus p95_span
#> 1               0.01           535           3           23.9       18
#> 2               0.02           535           3           23.9       18
#> 3               0.05           535           3           23.9       18
#> 4               0.10           535           3           23.9       18
#> 5               0.20           535           3           23.9       18
#>   pc_miss_over_5
#> 1           0.00
#> 2           0.33
#> 3           4.54
#> 4          18.29
#> 5          29.91
```

**Raising the threshold buys very little resolution and pays for it
steeply in coverage.** The number of groups barely moves, because the
same peptides are pooled either way and only the size of the group they
land in changes; the median group holds three candidates at every value
in the table. What does move is the miss rate, which is the argument for
keeping the threshold low.

The default of `min_candidate_prob = 0.02` is not taken from this
subset, which is too small to settle it — 410 proteins, and the sweep
above rests on 1498 pooled PSMs. It comes from applying the same
criterion to two complete datasets: the largest threshold at which under
1% of pooled peptides have more than a 5% chance of excluding the true
site. That gives 0.03 on a whole-proteome MaxQuant phosphoproteome and
0.02 on a complete Proteome Discoverer one, and the stricter of the two
is the default. At 0.02, 0.10% and 0.62% of pooled peptides respectively
exceed a 5% miss rate.

### The span guard

`max_group_span` caps how far apart the first and last candidate in a
group may be. It is a distance in residues, not a count of sites, and it
exists to stop a chain of overlapping missed-cleavage peptides fusing
genuinely unrelated regions of a protein into one group.

``` r

grouped %>%
  filter(!ptm_group_resolved %in% TRUE, !is.na(ptm_group_id)) %>%
  distinct(ptm_group_id, ptm_group_members) %>%
  pull(ptm_group_members) %>%
  strsplit(';') %>%
  sapply(function(pos) diff(range(as.numeric(pos)))) %>%
  max()
#> [1] 33
```

The default of 50 is insurance rather than a tuning parameter. The
widest group that forms unaided covers 33 residues here, 39 in the
complete MaxQuant dataset and 43 in the complete Proteome Discoverer
one, so on all three it never acts. Lowering it to 30 or below starts to
cut groups, and a peptide whose candidates end up on both sides of a cut
has no single group and is left unassigned rather than being forced into
one.

### What a pooled group cannot do

A pooled group is a feature you can quantify and test, and
`ptm_group_id` is usable directly as an `fcol` for
[`aggregateFeatures()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
in place of the `site` identifier built earlier. Four things do not
carry over.

A pooled group has no single residue, so it has no surrounding sequence.
[`add_site_sequence()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_site_sequence.md)
and the `site_seq` column it produces are undefined for one, and kinase
and motif enrichment cannot consume pooled groups as they stand.

A pooled group must never be labelled as though it were a residue.
`ptm_group_id` keeps the candidates visible, and
`ptm_group_n_candidates` is the column to carry through to any label
drawn elsewhere — a group reported as `S32` when it is one of four
candidates is a claim the data does not support.

``` r

group_labels <- grouped %>%
  filter(!is.na(ptm_group_id)) %>%
  distinct(ptm_group_id, ptm_group_members, ptm_group_n_candidates)

rbind(head(filter(group_labels, ptm_group_n_candidates == 1), 1),
      head(filter(group_labels, ptm_group_n_candidates > 1), 1))
#>                     ptm_group_id   ptm_group_members ptm_group_n_candidates
#> 2              Q8BTI8_n2_911;913             911;913                      1
#> 17 Q8BL97_n3_252;254;256;260;262 252;254;256;260;262                      5
```

`ptm_group_n_candidates` counts possible assignments rather than
residues, which is why the localised peptide above reports two positions
and a count of 1: it carries two phosphates, both placed, and there is
one way to assign them. A count above 1 marks a pooled group, and its
members are alternatives to one another rather than sites that co-occur.
The distinction is the whole content of the column, and it is what any
downstream label has to preserve.

Pooling sums intensity across residues, so a group containing two sites
that move in opposite directions reports their average, which may be no
change at all. The filter’s answer to this case is to report nothing;
pooling’s is to report something attenuated, and neither is the
site-level answer.

And the residues in a pooled group may also be quantified as resolved
sites in their own right, from peptides that did localise. The same
modification event then contributes to two features, which matters for
multiple testing and for any enrichment analysis that treats features as
independent.

``` r

resolved_sites <- grouped %>%
  filter(ptm_group_resolved %in% TRUE) %>%
  pull(ptm_group_members) %>%
  unique()

pooled_members <- grouped %>%
  filter(!ptm_group_resolved %in% TRUE, !is.na(ptm_group_id)) %>%
  distinct(ptm_group_id, ptm_group_members)

c(pooled_groups = nrow(pooled_members),
  overlapping_a_resolved_site = sum(sapply(
    strsplit(pooled_members$ptm_group_members, ';'),
    function(pos) any(pos %in% resolved_sites))))
#>               pooled_groups overlapping_a_resolved_site 
#>                         535                         484
```

**That overlap is the rule rather than the exception.** Nearly every
pooled group in this dataset contains at least one residue that is also
measured as a localised site, which is unsurprising — a region gets
pooled because some spectra could not place the modification there, not
because none ever could. The consequence is that pooled and localised
features are not two disjoint sets of measurements, and the number of
independent tests is smaller than the number of features. How much
smaller is not quantified: the overlap is easy to count, as above, but
its effect on FDR depends on the test applied downstream, and this
vignette does not attempt a correction for it.

Whether to group at all follows from the question. An analysis that
rests on individual residues — a motif, a kinase substrate, a site to
mutate — needs the localised sites and nothing else. An analysis asking
which proteins and which regions respond to a treatment is the one that
pays most for the filter, and it is the one grouping is for.

Finally, the identifier the PSMs will be summarised over. A site is only
meaningful with its protein attached — `S614` names a residue in some
protein, not a measurable thing — so the accession and the site name are
pasted together.

``` r

rowData(phospho[['psm_sites']])$site <- paste(
  rowData(phospho[['psm_sites']])$Master.Protein.Accessions,
  gsub(' ', '', rowData(phospho[['psm_sites']])$ptm_name),
  sep = '_')

head(unique(rowData(phospho[['psm_sites']])$site))
#> [1] "Q8BTI8_S911;S913"   "Q9CYI0_S55"         "Q8BTI8_S2036;S2037"
#> [4] "Q6ZQ58_S302"        "Q8BTI8_S1457;S1458" "Q8BTI8_S268;S270"
```

## Summarising to sites and proteins

With an identifier per fraction, both are summarised the same way as any
TMT experiment: drop the PSMs with missing values, sum the rest, and
log-transform. The enriched fraction is summarised over the site
identifier, the total fraction over the master protein.

``` r

phospho[['psm_complete']] <- filterNA(phospho[['psm_sites']], 0)
total[['psm_complete']] <- filterNA(total[['psm_quality']], 0)

phospho <- aggregateFeatures(phospho, i = 'psm_complete', fcol = 'site',
                             name = 'site_raw', fun = base::colSums)
total <- aggregateFeatures(total, i = 'psm_complete',
                           fcol = 'Master.Protein.Accessions',
                           name = 'protein_raw', fun = base::colSums)

phospho[['site_raw']] <- logTransform(phospho[['site_raw']], base = 2)
total[['protein_raw']] <- logTransform(total[['protein_raw']], base = 2)

c(sites = nrow(phospho[['site_raw']]), proteins = nrow(total[['protein_raw']]))
#>    sites proteins 
#>     1348      419
```

A peptide carrying two phosphates becomes one feature named for both,
not two features of one site each, because that is what was measured:
the quantification belongs to the doubly-modified form of that peptide,
and there is no way to divide it between the two residues. Peptides
carrying one of the sites alone, if they were observed, form their own
separate feature.

``` r

table(sites_per_feature = lengths(
  strsplit(rowData(phospho[['site_raw']])$ptm_name, '; ')))
#> sites_per_feature
#>   1   2   3   4 
#> 994 287  65   2
```

This means the same residue can appear in more than one feature, and
that a feature’s abundance is not the abundance of any single site.
Treating each feature as an independent observation of its residues
would double-count them.

## Normalising the enriched fraction

Unequal loading between channels has to be corrected in a PTM experiment
just as in any other. The usual correction — centring every sample on a
common median, as `QFeatures::normalize(method = 'diff.median')` does —
assumes that most features are unchanged between samples. In the total
fraction that is reasonable. In the enriched fraction it is not: a
treatment that changes phosphorylation globally will shift the median of
the enriched data, and centring it away removes exactly the effect the
experiment was run to detect.

The total fraction gives the reference instead. It comes from the same
labelled pool, so it carries the same loading differences, and its
median is a legitimate estimate of them.
[`get_medians()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_medians.md)
extracts the per-sample medians, and
[`center_normalise_to_ref()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/center_normalise_to_ref.md)
applies them to whichever assay you pass.

``` r

total_medians <- get_medians(total[['protein_raw']])

total[['protein']] <- center_normalise_to_ref(
  total[['protein_raw']], total_medians, on_log_scale = TRUE)

phospho[['site']] <- center_normalise_to_ref(
  phospho[['site_raw']], total_medians, on_log_scale = TRUE)
```

[`center_normalise_to_ref()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/center_normalise_to_ref.md)
with the default `center_to_zero = FALSE` subtracts each median only
after the medians have themselves been centred on their mean, so it
applies the *relative* differences between samples and leaves the
overall level of the data alone. That is why the same reference can be
applied to two assays on quite different scales.

``` r

plot_quant(phospho[['site_raw']], log2transform = FALSE, method = 'density') +
  theme_biomasslmb(base_size = 9) +
  theme(legend.key.height = unit(3, 'mm')) +
  labs(x = 'Site abundance (log2)', title = 'Before normalisation')

plot_quant(phospho[['site']], log2transform = FALSE, method = 'density') +
  theme_biomasslmb(base_size = 9) +
  theme(legend.key.height = unit(3, 'mm')) +
  labs(x = 'Site abundance (log2)', title = 'After normalisation')
```

![Site-level abundance distributions, before and after
normalisation](PTM_site_quantification_files/figure-html/unnamed-chunk-33-1.png)![Site-level
abundance distributions, before and after
normalisation](PTM_site_quantification_files/figure-html/unnamed-chunk-33-2.png)

Site-level abundance distributions, before and after normalisation

## Site abundance and occupancy

A site’s abundance can change for two reasons: the fraction of the
protein that is modified has changed, or the amount of the protein has
changed. Only the first is a change in signalling; the second would show
up just as well in the total fraction.

Subtracting the protein’s log abundance from the site’s separates them.
The result is often called occupancy or stoichiometry, though neither
name is quite right — an enrichment step does not preserve absolute
proportions, so this is a relative measure, comparable between samples
but not interpretable as “x% of the protein is modified”.

``` r

site_protein <- rowData(phospho[['site']])$Master.Protein.Accessions
has_protein <- site_protein %in% rownames(total[['protein']])

occupancy <- assay(phospho[['site']])[has_protein, ] -
  assay(total[['protein']])[site_protein[has_protein], ]

c(sites = nrow(phospho[['site']]), with_a_protein_estimate = sum(has_protein))
#>                   sites with_a_protein_estimate 
#>                    1348                    1321
```

Not every site gets one. A protein can be quantified in the enriched
fraction and not the total, which happens most for low-abundance
proteins, where the enrichment is the only reason they were seen at all
— so the sites that lose their protein estimate are not a random subset.

## Testing

Both matrices can go into `limma` in the usual way. The design here has
a treatment and a timepoint, and the question is about the treatment, so
the timepoint enters as a blocking factor.

``` r

condition <- factor(phospho[['site']]$Condition,
                    levels = c('Control', 'Treated'))
timepoint <- factor(phospho[['site']]$Timepoint)
design <- model.matrix(~ timepoint + condition)

test_treatment <- function(quant) {
  fit <- eBayes(lmFit(quant, design))
  topTable(fit, coef = 'conditionTreated', number = Inf, sort.by = 'none')
}

site_res <- test_treatment(assay(phospho[['site']]))
protein_res <- test_treatment(assay(total[['protein']]))
occupancy_res <- test_treatment(occupancy)

data.frame(
  level = c('Site abundance', 'Occupancy'),
  tested = c(nrow(site_res), nrow(occupancy_res)),
  significant = c(sum(site_res$adj.P.Val < 0.05),
                  sum(occupancy_res$adj.P.Val < 0.05)))
#>            level tested significant
#> 1 Site abundance   1348           6
#> 2      Occupancy   1321           6
```

Those are the two answers you might report. The protein-level test is
not a third one — it is a different question over a different set of
features — but it is what tells you how far to trust the agreement
between these two. Only 1 of the 419 proteins quantified in the total
fraction changes detectably with treatment, so for most sites there is
nothing for the correction to remove, and the two site-level answers
mostly agree.

``` r

comparison <- data.frame(
  site = rownames(occupancy),
  protein = site_protein[has_protein],
  site_logFC = site_res$logFC[has_protein],
  site_padj = site_res$adj.P.Val[has_protein],
  occupancy_logFC = occupancy_res$logFC,
  occupancy_padj = occupancy_res$adj.P.Val) %>%
  mutate(
    protein_logFC = protein_res[protein, 'logFC'],
    protein_padj = protein_res[protein, 'adj.P.Val'],
    called = case_when(
      site_padj < 0.05 & occupancy_padj < 0.05 ~ 'both',
      site_padj < 0.05 ~ 'site only',
      occupancy_padj < 0.05 ~ 'occupancy only',
      TRUE ~ 'neither'))

table(comparison$called)
#> 
#>           both        neither occupancy only      site only 
#>              4           1313              2              2
```

The interesting cases are the disagreements.

``` r

comparison %>%
  filter(called != 'neither', called != 'both') %>%
  select(site, called, site_logFC, site_padj,
         occupancy_logFC, occupancy_padj, protein_logFC, protein_padj) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 2)))
#>           site         called site_logFC site_padj occupancy_logFC
#> 1  P15105_S343      site only       0.73     0.001            0.29
#> 2  Q3B7Z2_S377 occupancy only       0.27     0.090            0.31
#> 3 Q8K4S1_S1118      site only       0.37     0.044            0.41
#> 4  Q9D032_S347 occupancy only       0.20     0.250            0.27
#>   occupancy_padj protein_logFC protein_padj
#> 1          0.180         0.440      0.00057
#> 2          0.035        -0.039      0.95000
#> 3          0.059        -0.039      0.99000
#> 4          0.035        -0.063      0.87000
```

``` r

ggplot(comparison, aes(site_logFC, occupancy_logFC)) +
  geom_hline(yintercept = 0, colour = 'grey80') +
  geom_vline(xintercept = 0, colour = 'grey80') +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey50') +
  geom_point(aes(colour = called), size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = c(get_cat_palette(3), 'grey70'),
                      breaks = c('both', 'occupancy only', 'site only', 'neither'),
                      name = 'Significant at') +
  theme_biomasslmb(base_size = 9) +
  labs(x = 'Site abundance (log2 fold change)',
       y = 'Occupancy (log2 fold change)')
```

![Site-level against occupancy-level fold change. Points off the
diagonal are sites whose protein also
moved.](PTM_site_quantification_files/figure-html/unnamed-chunk-38-1.png)

Site-level against occupancy-level fold change. Points off the diagonal
are sites whose protein also moved.

A site called at one level and not the other is not a contradiction; the
two are answering different questions. A site that moves only in the
abundance test, while its protein moves by a similar amount, is a
protein-level change showing through — the modified fraction is constant
and more of the protein is present. A site that moves only in the
occupancy test is one where protein-level variation was adding noise,
and removing it made a modest change detectable. With so few
protein-level changes in this subset, the first of those cases appears
only once here; where a treatment moves protein abundances more widely,
it is correspondingly commoner.

Which of the two to report depends on the question. Occupancy is the
right answer for “did signalling through this site change”, and it is
what most phosphoproteomics is asking. Site abundance is the right
answer for “how much of this modified form is present”, which is the
relevant quantity when the modified form is itself the effector.
Reporting both, as above, costs nothing and makes the distinction
visible; reporting site abundance alone and describing it as a change in
phosphorylation is the mistake this section exists to prevent.

Note also that occupancy is only available for sites whose protein was
quantified in the total fraction, so restricting to it discards sites —
and, as noted above, not a random subset of them.

## Testing both levels in one model

Subtracting the protein from the site and testing the result treats the
protein estimate as if it were exact. It is not — it is a measurement
with its own uncertainty, made on the same sample — and the subtraction
throws that uncertainty away before `limma` ever sees it.

The alternative is to give `limma` both measurements and let it do the
subtraction as part of the model. Stack the protein and site abundances
into one matrix, one row per site, with a `type` factor saying which
kind of measurement each column holds.

``` r

site_mat <- assay(phospho[['site']])[has_protein, ]
protein_mat <- assay(total[['protein']])[site_protein[has_protein], ]

joint <- cbind(protein_mat, site_mat)
colnames(joint) <- paste(rep(c('protein', 'site'), each = ncol(site_mat)),
                         colnames(site_mat), sep = '.')

joint_sample <- rep(colnames(site_mat), 2)
joint_type <- factor(rep(c('protein', 'site'), each = ncol(site_mat)),
                     levels = c('protein', 'site'))
joint_condition <- factor(rep(phospho[['site']]$Condition, 2),
                          levels = c('Control', 'Treated'))
joint_timepoint <- factor(rep(phospho[['site']]$Timepoint, 2))

dim(joint)
#> [1] 1321   32
```

Every term that can act differently on the two measurement types is then
interacted with `type`. The treatment effect on the protein is
`conditionTreated`; the extra effect on the site, over and above the
protein’s, is `conditionTreated:typesite`. That interaction is the
quantity occupancy was constructed to isolate.

``` r

joint_design <- model.matrix(~ joint_timepoint + joint_condition * joint_type +
                               joint_timepoint:joint_type)

colnames(joint_design)
#> [1] "(Intercept)"                          
#> [2] "joint_timepointT2"                    
#> [3] "joint_conditionTreated"               
#> [4] "joint_typesite"                       
#> [5] "joint_conditionTreated:joint_typesite"
#> [6] "joint_timepointT2:joint_typesite"
```

The same sample supplies one protein column and one site column, so
those two are not independent observations. `duplicateCorrelation`
estimates a single consensus correlation within blocks, and `lmFit` uses
it to weight the fit. [Technical
replicates](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/gotcha_technical_replicates.md)
covers the same machinery applied to repeated acquisitions of one
sample.

``` r

corfit <- duplicateCorrelation(joint, joint_design, block = joint_sample)

round(corfit$consensus.correlation, 3)
#> [1] 0.233
```

``` r

joint_fit <- eBayes(lmFit(joint, joint_design, block = joint_sample,
                          correlation = corfit$consensus.correlation))

joint_res <- topTable(joint_fit, coef = 'joint_conditionTreated:joint_typesite',
                      number = Inf, sort.by = 'none')

table(joint = joint_res$adj.P.Val < 0.05,
      occupancy = occupancy_res$adj.P.Val < 0.05)
#>        occupancy
#> joint   FALSE TRUE
#>   FALSE  1313    1
#>   TRUE      2    5
```

**The two routes estimate the same fold changes and disagree only on the
uncertainty.** The interaction coefficient and the occupancy fold change
are the same number to machine precision, which they have to be — both
are the difference between the site’s response and the protein’s.

``` r

max(abs(joint_res$logFC - occupancy_res$logFC))
#> [1] 3.497203e-15
```

What differs is the standard error each is divided by, and that is
enough to move a few sites across the threshold in both directions. The
gain is modest here because the pairing is weak in this dataset — a
consensus correlation of 0.233 — but it is the block term that delivers
it, and dropping it gives back most of the difference.

``` r

unblocked <- eBayes(lmFit(joint, joint_design))

c(blocked = sum(joint_res$adj.P.Val < 0.05),
  unblocked = sum(topTable(unblocked, coef = 'joint_conditionTreated:joint_typesite',
                           number = Inf)$adj.P.Val < 0.05),
  occupancy = sum(occupancy_res$adj.P.Val < 0.05))
#>   blocked unblocked occupancy 
#>         7         4         6
```

The reason to reach for this rather than the subtraction is not the
handful of sites it moves. It is that the model can express designs the
subtraction cannot. A term interacted with `type` is a term allowed to
act differently on protein and site, so a covariate that affects protein
abundance and site occupancy in different ways can be adjusted for as
such, and a design with more than two groups yields one interaction
coefficient per group from a single fit rather than a subtraction
repeated per contrast.

The cost is the assumption `duplicateCorrelation` makes: one correlation
describes every block. Where some samples pair their two fractions much
more tightly than others — different enrichment batches, say — that
single number describes none of them well.

## Summary

A PTM experiment adds three steps to the whole-proteome workflow, and
each of them is a place to get the answer wrong:

- **Localisation.** The search engine’s site assignment is a
  probabilistic claim, scored separately from the peptide
  identification.
  [`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
  reads the scores and discards peptides whose sites cannot be placed
  confidently. Positions then have to be converted from peptide
  coordinates to protein coordinates with
  [`add_ptm_positions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ptm_positions.md),
  dropping the peptides that occur at more than one position in their
  protein. Discarding is not the only option:
  [`parse_ptm_candidates_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_pd.md)
  and
  [`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)
  keep the unlocalised peptides by pooling their candidate residues into
  a group with a stated chance of containing the true site, which buys
  back a quarter of the PSMs here at the price of a feature that names a
  region rather than a residue.
- **Normalisation.** The enriched fraction’s own median is not a valid
  reference, because a global change in modification is signal rather
  than a loading artefact.
  [`get_medians()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_medians.md)
  on the total fraction and
  [`center_normalise_to_ref()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/center_normalise_to_ref.md)
  on the enriched one keeps the loading correction and leaves the
  biology.
- **Occupancy.** A change in site abundance is a change in modification
  only once the protein’s own change has been divided out. Dividing it
  out by subtraction and testing the result is the direct route; fitting
  protein and site abundances jointly, with a `type` interaction and
  `duplicateCorrelation` on the sample, reaches the same fold changes
  with an uncertainty that accounts for the pairing, and extends to
  designs the subtraction cannot express.

Two of the steps in this vignette are shared with any TMT experiment and
are covered in more detail elsewhere: PSM filtering and summarisation in
[TMT QC and
summarisation](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/TMT_PSM_QC_Summarisation.md),
and testing in [Data exploration and statistical
testing](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/exploration_and_statistical_testing.md).
If the design spans more than one plex, both fractions need bringing
onto a common scale before any of the above, which [multi-plex
TMT](https://lmb-mass-spec-compbio.github.io/biomasslmb/articles/TMT_multiplex.md)
covers — the bridge correction is applied to the enriched and total
assays separately, using each one’s own bridge channels.

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
#>  [3] ggplot2_4.0.3               biomasslmb_0.1.0           
#>  [5] QFeatures_1.20.0            MultiAssayExperiment_1.36.2
#>  [7] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [9] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#> [11] IRanges_2.44.0              S4Vectors_0.48.1           
#> [13] BiocGenerics_0.56.0         generics_0.1.4             
#> [15] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.3.0               rlang_1.3.0             magrittr_2.0.5         
#>  [4] clue_0.3-68             cleaver_1.48.0          otel_0.2.0             
#>  [7] compiler_4.5.3          RSQLite_3.53.3          png_0.1-9              
#> [10] systemfonts_1.3.2       vctrs_0.7.3             reshape2_1.4.5         
#> [13] stringr_1.6.0           ProtGenerics_1.42.0     pkgconfig_2.0.3        
#> [16] crayon_1.5.3            fastmap_1.2.0           backports_1.5.1        
#> [19] XVector_0.50.0          labeling_0.4.3          rmarkdown_2.32         
#> [22] visdat_0.6.0            ragg_1.5.2              purrr_1.2.2            
#> [25] bit_4.6.0               xfun_0.61               cachem_1.1.0           
#> [28] jsonlite_2.0.0          blob_1.3.0              DelayedArray_0.36.1    
#> [31] cluster_2.1.8.2         R6_2.6.1                bslib_0.12.0           
#> [34] stringi_1.8.9           RColorBrewer_1.1-3      genefilter_1.92.0      
#> [37] jquerylib_0.1.4         Rcpp_1.1.2              knitr_1.52             
#> [40] BiocBaseUtils_1.12.0    Matrix_1.7-4            splines_4.5.3          
#> [43] igraph_2.3.3            tidyselect_1.2.1        abind_1.4-8            
#> [46] yaml_2.3.12             lattice_0.22-9          tibble_3.3.1           
#> [49] plyr_1.8.9              withr_3.0.3             KEGGREST_1.50.0        
#> [52] S7_0.2.2                evaluate_1.0.5          uniprotREST_1.0.0      
#> [55] desc_1.4.3              survival_3.8-6          Biostrings_2.78.0      
#> [58] pillar_1.11.1           corrplot_0.95           checkmate_2.3.4        
#> [61] scales_1.4.0            xtable_1.8-8            glue_1.8.1             
#> [64] lazyeval_0.2.3          tools_4.5.3             robustbase_0.99-7      
#> [67] annotate_1.88.0         fs_2.1.0                XML_3.99-0.24          
#> [70] grid_4.5.3              tidyr_1.3.2             MsCoreUtils_1.22.1     
#> [73] AnnotationDbi_1.72.0    naniar_1.1.0            cli_3.6.6              
#> [76] textshaping_1.0.5       S4Arrays_1.10.1         AnnotationFilter_1.34.0
#> [79] gtable_0.3.6            DEoptimR_1.2-1          sass_0.4.10            
#> [82] digest_0.6.39           SparseArray_1.10.10     htmlwidgets_1.6.4      
#> [85] farver_2.1.2            memoise_2.0.1           htmltools_0.5.9        
#> [88] pkgdown_2.2.1           lifecycle_1.0.5         httr_1.4.9             
#> [91] statmod_1.5.2           bit64_4.8.6             MASS_7.3-65
```

Frankenfield, Ashley M., Jiawei Ni, Mustafa Ahmed, and Ling Hao. 2022.
“Protein Contaminants Matter: Building Universal Protein Contaminant
Libraries for DDA and DIA Proteomics.” *Journal of Proteome Research* 21
(9): 2104–13. <https://doi.org/10.1021/acs.jproteome.2c00145>.

McAlister, Graeme C., David P. Nusinow, Mark P. Jedrychowski, et al.
2014. “MultiNotch MS3 Enables Accurate, Sensitive, and Multiplexed
Detection of Differential Expression Across Cancer Cell Line Proteomes.”
*Analytical Chemistry* 86 (14): 7150–58.
<https://doi.org/10.1021/ac502040v>.
