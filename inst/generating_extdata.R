library(here)
library(dplyr)
# file from smb://istore.lmb.internal/mass_spec_share/Monday_Data/Cat/Cat_research/HYE/forTomSmith/DIA_Eclipse/DIANN
x <- read.delim('~/Downloads/DIANN/report.tsv')
set.seed(42)
random_500_proteins <- x %>% pull(Protein.Group) %>% unique() %>% sample(500)
x %>%
  filter(Protein.Group %in% random_500_proteins) %>%
  write.table(here('inst/extdata/diann_report.tsv'), sep='\t', row.names=FALSE, quote=FALSE)

# file from smb://istore.lmb.internal/mass_spec_share/Monday_Data/Cat/Cat_research/HYE/forTomSmith/TMT_SPSMS3
x <- read.delim('~/Downloads/TMT_SPSMS3/TMT_HYE_SPSMS3_PSMs.txt')
set.seed(42)
random_500_proteins <- x %>% pull(Master.Protein.Accessions) %>% unique() %>% sample(500)
x %>%
  filter(Master.Protein.Accessions %in% random_500_proteins) %>%
  write.table(here('inst/extdata/tmt_pd_PSMs.tsv'), sep='\t', row.names=FALSE, quote=FALSE)


# PD/Sequest HT PeptideGroups export from a real DDA LFQ TurboID (proximity
# biotinylation) experiment, searched against the Hao lab's '0602_Universal
# Contaminants' database. Collapsed here to one bait's biotin-labelled vs
# unlabelled control samples (2 conditions x 3 replicates), keeping just the
# columns the LFQ-DDA vignette uses, then truncated to 500 proteins.
x <- read.delim('~/Downloads/DDA/PeptideGroups.txt', check.names = FALSE)
keep_cols <- c(
  'Checked', 'Confidence', 'Annotated Sequence', 'Modifications', 'Contaminant',
  'Number of Protein Groups', 'Number of Proteins', 'Number of PSMs',
  'Master Protein Accessions', 'Protein Accessions',
  'Number of Missed Cleavages', 'Theo MHplus in Da',
  'Abundance F2 Sample BCAM biotin', 'Abundance F4 Sample BCAM biotin', 'Abundance F6 Sample BCAM biotin',
  'Abundance F1 Sample BCAM control', 'Abundance F3 Sample BCAM control', 'Abundance F5 Sample BCAM control',
  'Search Engine Rank by Search Engine Sequest HT',
  'Delta M in ppm by Search Engine Sequest HT',
  'RT in min by Search Engine Sequest HT'
)
x <- x[, keep_cols]
set.seed(42)
random_500_proteins <- x %>% pull(`Master Protein Accessions`) %>% unique() %>% sample(500)
x %>%
  filter(`Master Protein Accessions` %in% random_500_proteins) %>%
  write.table(here('inst/extdata/lfq_dda_pd_PeptideGroups.txt'), sep='\t', row.names=FALSE, quote=FALSE)


# DIA-NN report.parquet and sample metadata for a public plasma proteomics
# case series comparing Mpox, COVID-19 and healthy controls (Wang et al.
# 2022, EMBO Mol Med). Searched with a contaminants FASTA included directly
# in the search database, so contaminant proteins are identifiable via a
# `Cont_` prefix carried through into Protein.Ids/Protein.Group. Used in full
# (no subsetting/anonymisation needed) for the LFQ-DIA vignette; also used,
# with accompanying teaching material, in the Proteomics Data Analysis Using
# R course (bioinformatics-core-shared-training/Proteomics_stats_R).
file.copy(
  '~/git_repos/training/course_expression_proteomics/materials/dia/data/monkeypox_plasma_proteomes.parquet',
  here('inst/extdata/monkeypox_plasma_proteomes.parquet'), overwrite = TRUE)
file.copy(
  '~/git_repos/training/course_expression_proteomics/materials/dia/data/monkeypox_metadata.xlsx',
  here('inst/extdata/monkeypox_metadata.xlsx'), overwrite = TRUE)
