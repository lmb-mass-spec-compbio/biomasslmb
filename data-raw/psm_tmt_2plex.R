# Build `psm_tmt_2plex` and `tmt_2plex_design`: a real, small PD PSM-level
# dataset spread over two TMTpro 18plex plexes, each carrying a pooled bridge
# channel. Used in Part C of the TMT_PSM_QC_Summarisation vignette, which
# covers the multi-plex path: per-plex processing, joining the plexes and
# bridge correction. The single-plex datasets (`psm_tmt_clock`,
# `psm_tmt_per2_mq`) cannot demonstrate any of that.
#
# Source: 2026_11 project (real MRC-LMB proteomics data), subsetted to a
# teaching-sized dataset. The design is a 3 x 3 factorial: three cell lines
# (a wild type and two knockouts) crossed with three drug treatment
# timepoints, three replicates each, plus one pooled bridge channel per plex.
# Genotype and timepoint labels are generic, since the experiment is
# unpublished.

library(dplyr)

set.seed(42)

raw_dir <- "~/git_repos/projects/2026/2026_11_mTOR_KOs_Irene_Diaz/raw"

psm_inf <- file.path(raw_dir, "2671070900_Irene_TMT_PSMs.txt.gz")
design_inf <- file.path(raw_dir, "2671070900_Irene_TMT_200326.xlsx")

# ---- Experimental design -------------------------------------------------

labelling <- openxlsx::read.xlsx(design_inf, sheet = 2)
samples <- openxlsx::read.xlsx(design_inf, sheet = 3) %>%
  mutate(sample_n = as.character(sample_n))

tmt_2plex_design <- labelling %>%
  # channels with no sample_n were not labelled in that plex
  filter(!is.na(sample_n)) %>%
  left_join(samples, by = "sample_n") %>%
  tidyr::separate(sample, into = c("Genotype", "Time", "Replicate"),
                  sep = " ", fill = "right") %>%
  mutate(
    Genotype = recode(Genotype, "C" = "WT", "PDCD4" = "KO1", "SERBP1" = "KO2"),
    Time = recode(Time, "DMSO" = "0h", "1DAY" = "24h"),
    is_bridge = grepl("^Pool", sample_n),
    Genotype = ifelse(is_bridge, "Bridge", Genotype),
    Time = ifelse(is_bridge, NA, Time),
    sample = ifelse(is_bridge,
                    paste0("Bridge_", Plex),
                    paste(Genotype, Time, Replicate, sep = "_"))
  ) %>%
  # readQFeatures() splits the PSMs on the `runCol` column of the input and
  # matches channels to runs via a colData column of the same name
  mutate(runCol = Plex, Plex = factor(Plex)) %>%
  select(Sample = sample, quantCols = tag, runCol, Plex,
         Genotype, Time, Replicate) %>%
  arrange(Plex, Genotype, Time, Replicate)

# ---- Raw PSMs ------------------------------------------------------------

infdf <- read.delim(psm_inf)

# Proteome Discoverer file indices 6 and 7 are the two plexes
infdf$Plex <- ifelse(grepl("^F6\\.", infdf$File.ID), "1", "2")

abundance_cols <- grep("^Abundance\\.", colnames(infdf))
colnames(infdf)[abundance_cols] <- sub("^Abundance\\.", "",
                                       colnames(infdf)[abundance_cols])

stopifnot(all(tmt_2plex_design$quantCols %in% colnames(infdf)))

# ---- Subsample to a teaching-sized dataset -------------------------------

# 1) Proteins with at least two PSMs in *both* plexes, so that every retained
#    protein has something to contribute to the plex merge and the bridge
#    correction. All PSMs are kept for each selected protein.
in_both_plexes <- infdf %>%
  filter(Contaminant == "False", Number.of.Proteins == 1) %>%
  count(Master.Protein.Accessions, Plex) %>%
  filter(n >= 2) %>%
  count(Master.Protein.Accessions) %>%
  filter(n == 2)

keep_proteins <- sample(in_both_plexes$Master.Protein.Accessions, 300)

# 2) A sample of contaminant proteins, so the contaminant filtering step has
#    something real to remove
contaminant_proteins <- unique(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "True"])
keep_contaminants <- sample(contaminant_proteins, 10)

# 3) A slice of PSMs without a unique master protein, so that filter has
#    something real to remove too
non_unique_ix <- which(infdf$Number.of.Proteins > 1)
keep_non_unique_ix <- sample(non_unique_ix, 300)

psm_tmt_2plex <- infdf %>%
  filter(Master.Protein.Accessions %in% c(keep_proteins, keep_contaminants) |
           row_number() %in% keep_non_unique_ix) %>%
  distinct()

cat(sprintf("psm_tmt_2plex: %d PSMs, %d master proteins\n",
            nrow(psm_tmt_2plex),
            length(unique(psm_tmt_2plex$Master.Protein.Accessions))))
print(table(psm_tmt_2plex$Plex))

usethis::use_data(psm_tmt_2plex, overwrite = TRUE)
usethis::use_data(tmt_2plex_design, overwrite = TRUE)

# ---- Matching peptide-level output ---------------------------------------

# update_peptide_assignments() takes the peptide-level export, which carries a
# single consistent set of protein assignments per sequence. Only the
# assignment columns are needed, and only for the sequences retained above.
pep_inf <- file.path(raw_dir, "2671070900_Irene_TMT_PeptideGroups.txt.gz")

peptide_groups <- read.delim(pep_inf) %>%
  filter(Sequence %in% psm_tmt_2plex$Sequence) %>%
  select(Sequence, Protein.Accessions, Master.Protein.Accessions,
         Positions.in.Master.Proteins)

out_inf <- gzfile("inst/extdata/tmt_2plex_PeptideGroups.txt.gz", "w")
write.table(peptide_groups, out_inf, sep = "\t", quote = FALSE,
            row.names = FALSE)
close(out_inf)

cat(sprintf("tmt_2plex_PeptideGroups.txt.gz: %d peptides\n",
            nrow(peptide_groups)))
