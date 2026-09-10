# Build `psm_tmt_factorial` and `tmt_factorial_design`: a real MaxQuant
# TMT18plex PSM-level (evidence.txt) dataset from a whole-proteome experiment
# with a crossed 3 x 2 design, used as the MaxQuant worked example in the
# TMT_PSM_QC_Summarisation vignette and as the factorial worked example in
# exploration_and_statistical_testing.
#
# The package's other MaxQuant TMT dataset, `psm_tmt_per2_mq`, is a bait
# pulldown, so it cannot stand in for a whole-proteome comparison. This one
# can, and its design carries a second thing no other packaged dataset has: a
# nuisance factor (cell line) crossed with the factor of interest (treatment),
# which is what makes it possible to show the cost of ignoring one.
#
# Source: project 2025_30. The experiment ran a total and a phospho-enriched
# fraction on the same plex, split here on the `Raw.file` name; only the total
# fraction is kept.
#
# Anonymisation, following the same reasoning as data-raw/psm_tmt_per2_mq.R:
# the quantification is just numbers, but the design labels, the `Raw.file`
# column and a set of custom search-database entries all identify the specific
# project, person and tagged construct rather than describing the (deliberately
# generic) teaching experiment. All three are replaced below. Label maps are
# derived from the structure of the design rather than written out as literals,
# so this public script never records the original names. See
# data-raw/source_paths.R for why the input paths are globbed.

library(dplyr)

set.seed(42)

source("data-raw/source_paths.R")

# ---- Experimental design -------------------------------------------------

exp_design <- openxlsx::read.xlsx(
  source_path("2025/2025_30_*/raw/*_TMT18_*.xlsx"),
  startRow = 7, colNames = FALSE) %>%
  setNames(c("sample_name", "sample_conditions", "reporter_n", "tag",
             "reporter_mass")) %>%
  tidyr::separate(sample_conditions,
                  into = c("Genotype", "Treatment", "Replicate"), sep = " ")

stopifnot(nrow(exp_design) == 18)

# ---- Derive the anonymised label maps ------------------------------------

# Three cell lines: two carry a single tagged allele, the third carries both
# and is written as a pair separated by "/". That structure is what names the
# levels here, so the original line names are never written down.
source_genotypes <- unique(exp_design$Genotype)
stopifnot(length(source_genotypes) == 3)

double_line <- grep("/", source_genotypes, value = TRUE)
single_lines <- sort(setdiff(source_genotypes, double_line))
stopifnot(length(double_line) == 1, length(single_lines) == 2)

genotype_map <- setNames(c("Line_A", "Line_B", "Line_AB"),
                         c(single_lines, double_line))

# Two treatments, a vehicle and a compound. The vehicle sorts first.
source_treatments <- sort(unique(exp_design$Treatment))
stopifnot(length(source_treatments) == 2)

treatment_map <- setNames(c("Control", "Treated"), source_treatments)

exp_design$Genotype <- unname(genotype_map[exp_design$Genotype])
exp_design$Treatment <- unname(treatment_map[exp_design$Treatment])

# ---- Raw PSMs (evidence.txt) ---------------------------------------------

infdf <- read.delim(source_path("2025/2025_30_*/raw/evidence.txt.gz"))

# The plex was run as a total and a phospho-enriched fraction. Only the total
# fraction belongs here; the phospho fraction is a different vignette's
# subject, and `psm_tmt_phospho` already serves it.
infdf <- infdf[!grepl("phos", infdf$Raw.file), ]

# All Reporter-intensity-related columns (corrected/plain/count) - drop the
# plain and count variants entirely and keep the corrected ones, which are the
# ones any analysis should use.
all_reporter_cols <- which(grepl("^Reporter\\.intensity", colnames(infdf)))
abundance_cols <- which(grepl("^Reporter\\.intensity\\.corrected\\.",
                              colnames(infdf)))

reporter_n <- as.integer(gsub("^Reporter\\.intensity\\.corrected\\.", "",
                              colnames(infdf)[abundance_cols]))
stopifnot(setequal(reporter_n, exp_design$reporter_n))

exp_design <- exp_design[match(reporter_n, exp_design$reporter_n), ]

infdf <- infdf[, c(setdiff(seq_len(ncol(infdf)), all_reporter_cols),
                   abundance_cols)]

sample_names <- with(exp_design,
                     paste(Genotype, Treatment, Replicate, sep = "_"))

n <- ncol(infdf)
colnames(infdf)[(n - length(sample_names) + 1):n] <- sample_names

# ---- Exclude the tagged construct ----------------------------------------

# The search database included custom entries for the tagged constructs
# themselves (fusion proteins, not real UniProt accessions), so PSMs
# identifying them are labelled with pseudo-accessions naming the gene and the
# tag. Those are excluded entirely rather than relabelled, since the tagged
# gene is the one thing the treatment was designed to act on and naming it
# would identify the experiment.
construct_ids <- c("PER2MNG", "PER2HALO", "HALO001")
is_construct <- function(x) {
  Reduce(`|`, lapply(construct_ids, function(id) grepl(id, x, fixed = TRUE)))
}
infdf <- infdf[
  !is_construct(infdf$Proteins) &
    !is_construct(infdf$Leading.proteins) &
    !is_construct(infdf$Leading.razor.protein), ]

# ---- Subsample to a teaching-sized dataset -------------------------------

# 1) A sample of contaminant PSMs, so the contaminant-filtering step has
#    something real to remove
contaminant_proteins <- unique(
  infdf$Leading.razor.protein[infdf$Potential.contaminant == "+"])
keep_contaminant_proteins <- sample(contaminant_proteins,
                                    min(20, length(contaminant_proteins)))

# 2) A random sample of proteins with a unique master protein assignment.
#    Random rather than chosen for effect size: the vignette compares hit
#    counts between models, and those counts are only worth reporting if the
#    proteins behind them were not picked for the answer they give.
unique_master_proteins <- unique(
  infdf$Leading.razor.protein[
    infdf$Potential.contaminant == "" & infdf$Reverse == "" &
      !grepl(";", infdf$Leading.proteins)])
keep_unique_proteins <- sample(unique_master_proteins,
                               min(1500, length(unique_master_proteins)))

# 3) A slice of PSMs without a unique master protein
non_unique_ix <- which(grepl(";", infdf$Leading.proteins))
keep_non_unique_ix <- sample(non_unique_ix, min(200, length(non_unique_ix)))

# 4) A slice of decoy hits, so the decoy filter is doing real work too
decoy_ix <- which(infdf$Reverse != "")

psm_tmt_factorial <- infdf %>%
  filter(
    Leading.razor.protein %in% c(keep_contaminant_proteins,
                                 keep_unique_proteins) |
      row_number() %in% c(keep_non_unique_ix, decoy_ix)
  ) %>%
  distinct()

cat(sprintf("psm_tmt_factorial: %d PSMs, %d unique Leading.razor.protein\n",
            nrow(psm_tmt_factorial),
            length(unique(psm_tmt_factorial$Leading.razor.protein))))

# ---- Trim columns that cannot be used here --------------------------------

# `evidence.txt` cross-references MaxQuant's other output tables by row ID.
# None of those tables are packaged, so these columns are pointers to nothing,
# and between them they account for most of the object's size.
dead_reference_cols <- c("MS.MS.IDs", "MS.MS.scan.numbers", "MS3.scan.numbers",
                         "MS.MS.scan.number", "Protein.group.IDs", "Peptide.ID",
                         "Mod..peptide.ID", "Best.MS.MS", "id",
                         "Oxidation..M..site.IDs", "Phospho..STY..site.IDs")

psm_tmt_factorial <- psm_tmt_factorial[
  , setdiff(colnames(psm_tmt_factorial), dead_reference_cols)]

# Anything constant across every retained PSM carries no information. On this
# subset that removes the Phospho columns (this is the total fraction, searched
# without a phospho modification) and the match-between-runs columns (it was
# not enabled), rather than anything a QC step could have used.
constant_cols <- names(psm_tmt_factorial)[
  sapply(psm_tmt_factorial, function(x) length(unique(x[!is.na(x)])) <= 1)]

cat("Dropping constant columns:", paste(constant_cols, collapse = ", "), "\n")

psm_tmt_factorial <- psm_tmt_factorial[
  , setdiff(colnames(psm_tmt_factorial), constant_cols)]

stopifnot(all(c("Reverse", "Potential.contaminant", "Leading.razor.protein",
                "Proteins", "PIF") %in% colnames(psm_tmt_factorial)))

# ---- Anonymise: Raw.file embeds the project ID, analyst and date ----------

source_runs <- sort(unique(psm_tmt_factorial$Raw.file))
run_map <- setNames(sprintf("run_%02d", seq_along(source_runs)), source_runs)
psm_tmt_factorial$Raw.file <- unname(run_map[psm_tmt_factorial$Raw.file])

# ---- Design table ---------------------------------------------------------

tmt_factorial_design <- data.frame(
  Genotype = exp_design$Genotype,
  Treatment = exp_design$Treatment,
  Replicate = exp_design$Replicate,
  row.names = sample_names)

tmt_factorial_design$quantCols <- rownames(tmt_factorial_design)

stopifnot(all(table(tmt_factorial_design$Genotype,
                    tmt_factorial_design$Treatment) == 3))

usethis::use_data(psm_tmt_factorial, overwrite = TRUE)
usethis::use_data(tmt_factorial_design, overwrite = TRUE)
