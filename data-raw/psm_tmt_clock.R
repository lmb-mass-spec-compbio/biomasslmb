# Build `psm_tmt_clock` and `tmt_clock_design`: a real, small PD TMT12plex
# PSM-level dataset with a genuine 2-condition (Control vs Mutant) design,
# used in the TMT_PSM_QC_Summarisation vignette. Unlike the legacy
# `psm_tmt_total` (inherited from camprotR, no experimental design attached),
# this dataset lets the vignette demonstrate design-aware steps: PCA by
# condition, missingness-by-condition, differential abundance testing.
#
# Source: project 2026_04 (real MRC-LMB proteomics data), subsetted and
# reduced to a teaching-sized dataset. See data-raw/source_paths.R for why the
# input paths are globbed rather than written out.

library(dplyr)

set.seed(42)

source("data-raw/source_paths.R")

# ---- Experimental design -----------------------------------------------

exp_design <- openxlsx::read.xlsx(
  source_path("2026/2026_04_*/raw/labelling_scheme.xlsx")) %>%
  tidyr::separate(sample, into = c(NA, "Condition", "Replicate"), sep = "", remove = FALSE) %>%
  mutate(Condition = recode(Condition, "C" = "Control", "M" = "Mutant"))

# ---- Raw PSMs ------------------------------------------------------------

infdf <- read.delim(source_path("2026/2026_04_*/raw/*_TMT_12labels_PSMs.txt"))

abundance_cols <- which(grepl("^Abundance\\.", colnames(infdf)))

# map TMT tag -> sample name (e.g. "126" -> "C1") and rename abundance columns
tag2sample <- setNames(exp_design$sample, exp_design$tag)
colnames(infdf)[abundance_cols] <- gsub("^Abundance\\.", "", colnames(infdf)[abundance_cols])
colnames(infdf)[abundance_cols] <- tag2sample[colnames(infdf)[abundance_cols]]

# ---- Subsample to a teaching-sized dataset -------------------------------

# 1) A sample of contaminant PSMs (PD-flagged), so the contaminant-filtering
#    step in the vignette has something real to remove
contaminant_proteins <- unique(infdf$Master.Protein.Accessions[infdf$Contaminant == "True"])
keep_contaminant_proteins <- sample(contaminant_proteins, min(20, length(contaminant_proteins)))

# 2) The proteins most differentially abundant between Control and Mutant in
#    the FULL dataset (identified via the same PD-filter -> sum-summarise ->
#    limma::treat pipeline used in the vignettes, run once on the complete
#    125k-PSM file). These are real quantification values for real proteins,
#    not fabricated - included deliberately so the teaching subsample still
#    shows genuine biological signal (e.g. in PCA/volcano plots) rather than
#    only a random, likely-null background. NOTE: because the vignettes only
#    test this ~1000-protein subsample rather than the full ~5000-protein
#    dataset, the multiple-testing burden - and hence the FDR-adjusted
#    p-values - will differ from an analysis of the complete experiment.
top_hit_proteins <- c(
  "Q9Z0V7", "P01837", "E9Q4P1", "Q64471", "P56387", "P03980", "Q8BWU8",
  "Q78YY6", "Q8K070", "P24472", "Q3V3R1", "Q91W61", "Q8CEE7", "O89086",
  "Q9WUC3", "S4R1W1")

# 3) A random sample of further proteins with a unique master protein
#    assignment, keeping *all* PSMs for each selected protein so per-protein
#    aggregation behaves realistically
unique_master_proteins <- unique(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "False" & infdf$Number.of.Proteins == 1]
)
unique_master_proteins <- setdiff(unique_master_proteins, top_hit_proteins)
keep_unique_proteins <- c(top_hit_proteins, sample(unique_master_proteins, 600))

# 4) A slice of PSMs without a unique master protein, so the
#    unique-master-protein filtering step also removes something real
non_unique_ix <- which(infdf$Number.of.Proteins > 1)
keep_non_unique_ix <- sample(non_unique_ix, min(300, length(non_unique_ix)))

psm_tmt_clock <- infdf %>%
  filter(
    Master.Protein.Accessions %in% c(keep_contaminant_proteins, keep_unique_proteins) |
      row_number() %in% keep_non_unique_ix
  ) %>%
  distinct()

cat(sprintf("psm_tmt_clock: %d PSMs, %d unique Master.Protein.Accessions\n",
            nrow(psm_tmt_clock), length(unique(psm_tmt_clock$Master.Protein.Accessions))))

# ---- Design table for the retained samples -------------------------------

tmt_clock_design <- exp_design %>%
  select(sample, tag, Condition, Replicate) %>%
  mutate(quantCols = sample) %>%
  tibble::column_to_rownames("sample")

usethis::use_data(psm_tmt_clock, overwrite = TRUE)
usethis::use_data(tmt_clock_design, overwrite = TRUE)
