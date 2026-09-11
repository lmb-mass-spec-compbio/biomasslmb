# Build `psm_tmt_phospho`, `psm_tmt_phospho_total` and `tmt_phospho_design`,
# plus `inst/extdata/tmt_phospho_proteome.fasta.gz`: a matched pair of real PD
# PSM-level datasets from a phosphoproteomics experiment, used in the
# PTM_site_quantification vignette.
#
# The experiment TMT-labelled each sample once, pooled, then split the pool
# into a phospho-enriched fraction and a total (non-enriched) fraction, which
# were acquired as separate runs. That is what makes the pair useful: the
# vignette needs both to demonstrate normalising the enriched data using the
# medians of the total data, and to express site abundance relative to the
# protein it sits on.
#
# Source: project 2025_17, second plex. Mouse fibroblasts, a vehicle control
# against a drug treatment, sampled at two timepoints with four replicates
# each. The second plex is used because the first has too few differential
# sites at this subsample size for the vignette to have anything to show at the
# end. Condition and timepoint labels are generic, since the experiment is
# unpublished. See data-raw/source_paths.R for why the input paths are globbed
# rather than written out.

library(dplyr)
library(tidyr)

devtools::load_all()

set.seed(42)

source("data-raw/source_paths.R")

total_psm_inf <- source_path("2025/2025_17_*_fibroblasts/raw/*_TMT_TP-*_PSMs.txt.gz")
phos_psm_inf <- source_path("2025/2025_17_*_fibroblasts/raw/*_phos_*_TMT1_PSMs.txt.gz")
design_inf <- source_path("2025/2025_17_*_fibroblasts/raw/experimental_design.xlsx")

# ---- Experimental design -------------------------------------------------

# The `Condition` column packs the treatment, the two timepoint components and
# the replicate into one string, e.g. "<treatment> <t1> + <t2> v<n>".
split_design <- openxlsx::read.xlsx(design_inf, sheet = 1) %>%
  filter(Plex == 2, !is.na(Condition), Condition != "is") %>%
  mutate(Condition = gsub(" \\+ | \\+", "+", Condition)) %>%
  separate(Condition, into = c("Condition", "arm", "Replicate"), sep = " ") %>%
  separate(arm, into = c("Timepoint", "post_treatment"), sep = "\\+")

# The two conditions sort with the treatment before the vehicle control, and
# the timepoints sort chronologically. Both orderings are properties of this
# design sheet, so assert the shape rather than trusting it silently.
source_conditions <- sort(unique(split_design$Condition))
source_timepoints <- sort(unique(split_design$Timepoint))

stopifnot(length(source_conditions) == 2, length(source_timepoints) == 2)

condition_map <- setNames(c("Treated", "Control"), source_conditions)
timepoint_map <- setNames(c("T1", "T2"), source_timepoints)

plex_design <- split_design %>%
  mutate(
    Condition = unname(condition_map[Condition]),
    Timepoint = unname(timepoint_map[Timepoint]),
    Replicate = sub("^v", "", Replicate),
    Sample = paste(Condition, Timepoint, Replicate, sep = "_")) %>%
  arrange(Condition, Timepoint, Replicate)

# The TMT tags the retained samples were labelled with, in sample order
tags <- plex_design$label

tmt_phospho_design <- plex_design %>%
  select(Sample, Condition, Timepoint, Replicate) %>%
  tibble::column_to_rownames("Sample")

tmt_phospho_design$quantCols <- rownames(tmt_phospho_design)

# ---- Raw PSMs ------------------------------------------------------------

# Each fraction was acquired over its own pair of PD file indices, one per
# plex, the second index being the second plex.
read_plex_2 <- function(inf, tags, sample_names) {
  infdf <- read.delim(inf)

  file_prefixes <- sort(unique(sub("\\..*", "", infdf$File.ID)))
  stopifnot(length(file_prefixes) == 2)
  infdf <- infdf[sub("\\..*", "", infdf$File.ID) == file_prefixes[2], ]

  colnames(infdf) <- sub("^Abundance\\.", "", colnames(infdf))
  stopifnot(all(tags %in% colnames(infdf)))

  # Drop the bridge channel and the unlabelled channel, then name the retained
  # tags for their samples
  all_tag_cols <- grep("^1(2[6-9]|3[0-5])(N|C)?$", colnames(infdf))
  infdf <- infdf[, c(setdiff(seq_len(ncol(infdf)), all_tag_cols),
                     match(tags, colnames(infdf)))]
  colnames(infdf)[(ncol(infdf) - length(tags) + 1):ncol(infdf)] <- sample_names

  infdf
}

total_df <- read_plex_2(total_psm_inf, tags, rownames(tmt_phospho_design))
phos_df <- read_plex_2(phos_psm_inf, tags, rownames(tmt_phospho_design))

# ---- Subsample to a teaching-sized dataset -------------------------------

# Keep proteins seen in both fractions, with enough PSMs in each that the site
# and protein estimates the vignette compares are both supported by more than a
# single observation. All PSMs are kept for each selected protein.
localisable <- phos_df %>%
  filter(Contaminant == "False", Number.of.Proteins == 1,
         !ptmRS.Best.Site.Probabilities %in% c("", "Inconclusive data"),
         grepl("Phospho", ptmRS.Best.Site.Probabilities)) %>%
  count(Master.Protein.Accessions, name = "n_phospho")

quantified <- total_df %>%
  filter(Contaminant == "False", Number.of.Proteins == 1) %>%
  count(Master.Protein.Accessions, name = "n_total")

shared <- merge(localisable, quantified) %>%
  filter(n_phospho >= 2, n_total >= 3)

keep_proteins <- sample(shared$Master.Protein.Accessions, 400)

subsample <- function(infdf) {
  # A sample of contaminant proteins and a slice of PSMs without a unique
  # master protein, so the routine filtering steps have something real to remove
  contaminant_proteins <- unique(
    infdf$Master.Protein.Accessions[infdf$Contaminant == "True"])
  keep_contaminants <- sample(contaminant_proteins,
                              min(10, length(contaminant_proteins)))

  non_unique_ix <- which(infdf$Number.of.Proteins > 1)
  keep_non_unique_ix <- sample(non_unique_ix, min(200, length(non_unique_ix)))

  infdf %>%
    filter(Master.Protein.Accessions %in% c(keep_proteins, keep_contaminants) |
             row_number() %in% keep_non_unique_ix) %>%
    distinct()
}

psm_tmt_phospho <- subsample(phos_df)
psm_tmt_phospho_total <- subsample(total_df)

cat(sprintf("psm_tmt_phospho: %d PSMs, %d master proteins\n",
            nrow(psm_tmt_phospho),
            length(unique(psm_tmt_phospho$Master.Protein.Accessions))))
cat(sprintf("psm_tmt_phospho_total: %d PSMs, %d master proteins\n",
            nrow(psm_tmt_phospho_total),
            length(unique(psm_tmt_phospho_total$Master.Protein.Accessions))))

usethis::use_data(psm_tmt_phospho, overwrite = TRUE)
usethis::use_data(psm_tmt_phospho_total, overwrite = TRUE)
usethis::use_data(tmt_phospho_design, overwrite = TRUE)

# ---- Proteome fasta for locating sites within proteins --------------------

# add_ptm_positions() re-digests the protein sequences to place each peptide,
# and hence each site, within its protein. Only the master proteins of the
# retained PSMs are needed, and shipping them means the vignette does not have
# to query UniProt while building.
fasta_accessions <- sort(unique(c(psm_tmt_phospho$Master.Protein.Accessions,
                                 psm_tmt_phospho_total$Master.Protein.Accessions)))
fasta_accessions <- fasta_accessions[
  fasta_accessions != "" & !grepl(";", fasta_accessions)]

fasta_tmp <- tempfile(fileext = ".fasta")
make_fasta(fasta_accessions, file = fasta_tmp)

fasta_out <- "inst/extdata/tmt_phospho_proteome.fasta.gz"
proteome <- Biostrings::readAAStringSet(fasta_tmp)
Biostrings::writeXStringSet(proteome, fasta_out, compress = TRUE)
file.remove(fasta_tmp)

cat(sprintf("tmt_phospho_proteome.fasta.gz: %d sequences of %d accessions\n",
            length(proteome), length(fasta_accessions)))
