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


# file from smb://istore.lmb.internal/mass_spec_share/Monday_Data/Cat/Cat_research/HYE/forTomSmith/DDA
x <- read.delim('~/Downloads/DDA/HFX_DDA_230min_HYE_PeptideGroups.txt')
set.seed(42)
random_500_proteins <- x %>% pull(Master.Protein.Accessions) %>% unique() %>% sample(500)
x %>%
  filter(Master.Protein.Accessions %in% random_500_proteins) %>%
  write.table(here('inst/extdata/lfq_dda_pd_PeptideGroups.txt'), sep='\t', row.names=FALSE, quote=FALSE)
