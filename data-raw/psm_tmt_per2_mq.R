# Build `psm_tmt_per2_mq` and `tmt_per2_mq_design`: a real, small MaxQuant
# TMT18plex PSM-level (evidence.txt) dataset with a genuine 2-condition
# (bait IP vs control IP) design, used as the MaxQuant worked example in
# the TMT_PSM_QC_Summarisation vignette, alongside the PD-based
# `psm_tmt_clock` example.
#
# Source: 2025_23_PER2_interactors_stoichiometry project (real MRC-LMB
# proteomics data). The full experiment includes a "Booster" channel and
# several unused/empty TMT channels which are dropped here to leave a clean
# 2-condition x 6-replicate design. The original condition labels (`P2H`
# for the bait pulldown, `HALO` for the control pulldown) and the
# `Raw.file` column (which embeds the analyst's name and the bait protein)
# are anonymised below to generic `IP`/`Control` labels and placeholder
# run names, since - unlike the protein accessions, which are just real
# quantification data - these two columns identify the specific project
# and person rather than describing the (deliberately generic) teaching
# experiment.

library(dplyr)

set.seed(42)

raw_dir <- "~/git_repos/projects/2025/2025_23_PER2_interactors_stoichiometry/raw"

# ---- Experimental design -----------------------------------------------

exp_design <- openxlsx::read.xlsx(
  file.path(raw_dir, "2103938239_Andrei_TMT12_210825.xlsx"), sheet = 2) %>%
  filter(!is.na(Sample)) %>%
  tidyr::separate(Sample, into = c("Condition", "Replicate"), remove = FALSE) %>%
  tibble::column_to_rownames("Sample")

# ---- Raw PSMs (evidence.txt) ---------------------------------------------

infdf <- read.delim(file.path(raw_dir, "evidence.txt"))

# All Reporter-intensity-related columns (corrected/plain/count) - drop the
# plain and count variants entirely, and restrict the corrected columns to
# just the real bait/control samples (drop Booster + empty channels)
all_reporter_cols <- which(grepl("^Reporter\\.intensity", colnames(infdf)))
abundance_cols <- which(grepl("^Reporter\\.intensity\\.corrected\\.", colnames(infdf)))

tag2sample <- setNames(rownames(exp_design), exp_design$Tag_n)
sample_names <- tag2sample[gsub("^Reporter\\.intensity\\.corrected\\.", "",
                                colnames(infdf)[abundance_cols])]

keep_samples <- intersect(sample_names, rownames(exp_design)[exp_design$Condition %in% c("P2H", "HALO")])
keep_abundance_cols <- abundance_cols[match(keep_samples, sample_names)]

infdf <- infdf[, c(setdiff(seq_len(ncol(infdf)), all_reporter_cols), keep_abundance_cols)]

# The abundance columns are now the last length(keep_samples) columns, in
# the order of keep_samples (via the match() above)
n <- ncol(infdf)
colnames(infdf)[(n - length(keep_samples) + 1):n] <- keep_samples

# ---- Subsample to a teaching-sized dataset -------------------------------

# The search database included custom entries for the tagged bait construct
# itself (a fusion protein, not a real UniProt accession), so PSMs
# identifying the bait are labelled with these pseudo-accessions rather
# than a UniProt ID. Since that would directly identify the gene and the
# tag used, these are excluded from the sampling pool entirely (not merely
# relabelled downstream).
bait_construct_ids <- c("PER2MNG", "PER2HALO", "HALO001")
is_bait_construct <- function(x) {
  Reduce(`|`, lapply(bait_construct_ids, function(id) grepl(id, x, fixed = TRUE)))
}
infdf <- infdf[
  !is_bait_construct(infdf$Proteins) &
    !is_bait_construct(infdf$Leading.proteins) &
    !is_bait_construct(infdf$Leading.razor.protein), ]

# 1) A sample of contaminant PSMs, so the contaminant-filtering step in the
#    vignette has something real to remove
contaminant_proteins <- unique(infdf$Leading.razor.protein[infdf$Potential.contaminant == "+"])
keep_contaminant_proteins <- sample(contaminant_proteins, min(20, length(contaminant_proteins)))

# 2) A random sample of proteins with a unique master protein assignment
unique_master_proteins <- unique(
  infdf$Leading.razor.protein[
    infdf$Potential.contaminant == "" & infdf$Reverse == "" &
      !grepl(";", infdf$Leading.proteins)])
keep_unique_proteins <- sample(unique_master_proteins, 600)

# 3) A slice of PSMs without a unique master protein
non_unique_ix <- which(grepl(";", infdf$Leading.proteins))
keep_non_unique_ix <- sample(non_unique_ix, min(200, length(non_unique_ix)))

psm_tmt_per2_mq <- infdf %>%
  filter(
    Leading.razor.protein %in% c(keep_contaminant_proteins, keep_unique_proteins) |
      row_number() %in% keep_non_unique_ix
  ) %>%
  filter(Reverse == "") %>%
  distinct()

cat(sprintf("psm_tmt_per2_mq: %d PSMs, %d unique Leading.razor.protein\n",
            nrow(psm_tmt_per2_mq), length(unique(psm_tmt_per2_mq$Leading.razor.protein))))

# ---- Design table for the retained samples -------------------------------

tmt_per2_mq_design <- exp_design[keep_samples, c("Condition", "Replicate")]

# ---- Anonymise: generic condition labels and sample names -----------------

condition_map <- c(P2H = "IP", HALO = "Control")
tmt_per2_mq_design$Condition <- condition_map[tmt_per2_mq_design$Condition]
new_sample_names <- paste(tmt_per2_mq_design$Condition, tmt_per2_mq_design$Replicate, sep = "_")

sample_rename <- setNames(new_sample_names, keep_samples)
colnames(psm_tmt_per2_mq)[match(keep_samples, colnames(psm_tmt_per2_mq))] <- sample_rename[keep_samples]

rownames(tmt_per2_mq_design) <- new_sample_names
tmt_per2_mq_design$quantCols <- rownames(tmt_per2_mq_design)

# ---- Anonymise: Raw.file embeds the analyst's name and the bait protein ---

run_map <- setNames(
  sprintf("run_%02d", seq_along(unique(psm_tmt_per2_mq$Raw.file))),
  unique(psm_tmt_per2_mq$Raw.file))
psm_tmt_per2_mq$Raw.file <- run_map[psm_tmt_per2_mq$Raw.file]

usethis::use_data(psm_tmt_per2_mq, overwrite = TRUE)
usethis::use_data(tmt_per2_mq_design, overwrite = TRUE)
