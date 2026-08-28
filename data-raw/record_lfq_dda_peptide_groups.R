# Records inst/extdata/lfq_dda_pd_PeptideGroups.txt, the Proteome Discoverer
# peptide-level export used by the LFQ-DDA QC and summarisation vignette.
#
# Source: 2026_05 project (real MRC-LMB proteomics data), subsetted to a
# teaching-sized dataset. A whole-proteome LFQ-DDA comparison between a wild
# type human cell line and a point mutant of the same line, three replicates
# each, searched with CHIMERYS. Condition labels are generic, since the
# experiment is unpublished.
#
# The abundance columns keep the file IDs Proteome Discoverer assigns
# (Abundance.F7.Sample ...), because relabelling them from a design table is
# the first thing the vignette does. F7-F12 are samples 1-6 of
# raw/Sample_Info.xlsx, i.e. the three wild type replicates followed by the
# three mutant replicates.

library(dplyr)

set.seed(42)

pep_inf <- "~/git_repos/projects/2026/2026_05_AMRC5_WT_vs_Mutant/raw/PD/2733983520_Roberta_20R_PeptideGroups.txt"

infdf <- read.delim(pep_inf)

# ---- Columns to retain ---------------------------------------------------

keep_cols <- c(
  "Checked", "Confidence", "Sequence", "Annotated.Sequence", "Modifications",
  "Contaminant", "Number.of.Protein.Groups", "Number.of.Proteins",
  "Number.of.PSMs", "Master.Protein.Accessions", "Protein.Accessions",
  "Number.of.Missed.Cleavages", "Theo.MHplus.in.Da",
  grep("^Abundance\\.F", colnames(infdf), value = TRUE),
  "PSM.Search.Engine.Rank.by.Search.Engine.CHIMERYS",
  "PSM.Delta.M.in.ppm.by.Search.Engine.CHIMERYS",
  "PSM.RT.in.min.by.Search.Engine.CHIMERYS")

stopifnot(all(keep_cols %in% colnames(infdf)))

infdf <- infdf[, keep_cols]

# ---- Subsample to a teaching-sized dataset -------------------------------

# 1) A sample of contaminant peptides, so the contaminant filtering step has
#    something real to remove
contaminant_proteins <- unique(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "True"])
keep_contaminants <- sample(contaminant_proteins,
                            min(20, length(contaminant_proteins)))

# 2) A random sample of proteins with a unique master protein assignment,
#    keeping every peptide for each so per-protein aggregation is realistic
unique_master_proteins <- unique(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "False" &
                                    infdf$Number.of.Protein.Groups == 1])
keep_proteins <- sample(unique_master_proteins, 500)

# 3) A slice of peptides without a unique master protein, so that filter also
#    removes something real
non_unique_ix <- which(infdf$Number.of.Protein.Groups > 1)
keep_non_unique_ix <- sample(non_unique_ix, min(300, length(non_unique_ix)))

peptide_groups <- infdf %>%
  filter(Master.Protein.Accessions %in% c(keep_contaminants, keep_proteins) |
           row_number() %in% keep_non_unique_ix) %>%
  distinct()

cat(sprintf("lfq_dda_pd_PeptideGroups.txt: %d peptides, %d master proteins\n",
            nrow(peptide_groups),
            length(unique(peptide_groups$Master.Protein.Accessions))))

write.table(peptide_groups, "inst/extdata/lfq_dda_pd_PeptideGroups.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
