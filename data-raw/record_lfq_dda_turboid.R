# Records inst/extdata/lfq_dda_pd_turboid_PeptideGroups.txt, the Proteome
# Discoverer peptide-level export used by Part B of
# vignettes/interactome_designs.Rmd and referenced from R/qf_data.R.
#
# Source: 2025_11 project (real MRC-LMB proteomics data), subsetted to a
# teaching-sized dataset. A TurboID proximity-labelling experiment in mouse
# embryonic stem cells, searched with Sequest HT. The vignette uses it for the
# presence/absence structure that LFQ interactome data has and TMT does not.
#
# The export covers two baits, six samples each. Only the first bait's arm is
# kept -- three biotin-treated samples against three untreated controls. One
# bait is enough for what the vignette shows, and the second would only add
# columns the vignette immediately drops.
#
# What is deliberately not carried over:
#   the bait names -- the source export names both baits in its abundance
#     column headers, and the source filename names a bait and a collaborator.
#     The experiment is unpublished and the bait identifies the project, so the
#     abundance columns are reduced to the file ID and the condition. Nothing
#     in the vignette needs the bait. The name is derived below rather than
#     written out, so this script does not record what the data drops. See
#     data-raw/source_paths.R for why the input path is globbed.
#
# Column names keep Proteome Discoverer's spacing (check.names = FALSE),
# because the vignette reads the file the way the instrument software wrote it.

set.seed(42)

source("data-raw/source_paths.R")

n_proteins <- 500

pep_inf <- source_path("2025/2025_11_*/raw/*_PeptideGroups.txt")

infdf <- read.delim(pep_inf, check.names = FALSE)

# ---- Pick one bait's six samples -----------------------------------------

abundance_cols <- grep("^Abundance F[0-9]+ Sample ",
                       colnames(infdf), value = TRUE)

file_no <- as.integer(sub("^Abundance F([0-9]+) .*$", "\\1", abundance_cols))
bait <- sub("^Abundance F[0-9]+ Sample (.*) (biotin|control)$", "\\1",
            abundance_cols)
condition <- sub("^.* ", "", abundance_cols)

# The first bait in acquisition order, identified by position rather than name.
keep_bait <- bait == bait[which.min(file_no)]
abundance_cols <- abundance_cols[keep_bait]
file_no <- file_no[keep_bait]
condition <- condition[keep_bait]

stopifnot(length(abundance_cols) == 6, sum(condition == "biotin") == 3)

# Biotin-treated samples first, then the controls, each in acquisition order.
col_order <- order(condition != "biotin", file_no)
abundance_cols <- abundance_cols[col_order]

# ---- Columns to retain ---------------------------------------------------

keep_cols <- c(
  "Checked", "Confidence", "Annotated Sequence", "Modifications", "Contaminant",
  "Number of Protein Groups", "Number of Proteins", "Number of PSMs",
  "Master Protein Accessions", "Protein Accessions",
  "Number of Missed Cleavages", "Theo MHplus in Da",
  abundance_cols,
  "Search Engine Rank by Search Engine Sequest HT",
  "Delta M in ppm by Search Engine Sequest HT",
  "RT in min by Search Engine Sequest HT")

stopifnot(all(keep_cols %in% colnames(infdf)))

infdf <- infdf[, keep_cols]

# Drop the bait from the abundance headers, keeping the file ID and condition.
colnames(infdf)[match(abundance_cols, colnames(infdf))] <- paste(
  sub("^(Abundance F[0-9]+ Sample) .*$", "\\1", abundance_cols),
  condition[col_order])

# ---- Subsample to a teaching-sized dataset -------------------------------

keep_proteins <- sample(unique(infdf[["Master Protein Accessions"]]),
                        n_proteins)

peptide_groups <- infdf[infdf[["Master Protein Accessions"]] %in%
                          keep_proteins, ]

cat(sprintf(
  "lfq_dda_pd_turboid_PeptideGroups.txt: %d peptides, %d master proteins\n",
  nrow(peptide_groups),
  length(unique(peptide_groups[["Master Protein Accessions"]]))))

write.table(peptide_groups, "inst/extdata/lfq_dda_pd_turboid_PeptideGroups.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
