# Build `psm_tmt_phospho_mq` plus `inst/extdata/tmt_phospho_mq_proteome.fasta.gz`:
# the phospho-enriched fraction of the same TMT18plex that `psm_tmt_factorial`
# takes its total fraction from, used as the MaxQuant worked example in the
# PTM_site_quantification vignette.
#
# The package's PD phospho example, `psm_tmt_phospho`, cannot stand in for it.
# The two search engines report localisation confidence in different columns
# and in different formats - MaxQuant writes probabilities inline in the
# peptide sequence, ptmRS writes them as a separate list of percentages - so a
# vignette that shows only one leaves the other's parser undemonstrated.
#
# Source: project 2025_30, the `_phos_` raw files. `data-raw/psm_tmt_factorial.R`
# takes the `_prot_` files from the same evidence.txt and the same plex, so the
# two datasets share `tmt_factorial_design` and no new design object is built
# here. Anonymisation follows that script exactly, and for the same reasons:
# the quantification is just numbers, but the design labels, the `Raw.file`
# column and the custom search-database entries identify the specific project,
# person and tagged construct. See data-raw/source_paths.R for why the input
# paths are globbed rather than written out.

library(dplyr)

set.seed(42)

source("data-raw/source_paths.R")

# ---- Experimental design -------------------------------------------------

# Identical to data-raw/psm_tmt_factorial.R: the same plex, split into two
# fractions. The labels are re-derived here rather than read back from the
# packaged `tmt_factorial_design`, so that this script does not depend on the
# order in which the two are built.
exp_design <- openxlsx::read.xlsx(
  source_path("2025/2025_30_*/raw/*_TMT18_*.xlsx"),
  startRow = 7, colNames = FALSE) %>%
  setNames(c("sample_name", "sample_conditions", "reporter_n", "tag",
             "reporter_mass")) %>%
  tidyr::separate(sample_conditions,
                  into = c("Genotype", "Treatment", "Replicate"), sep = " ")

stopifnot(nrow(exp_design) == 18)

source_genotypes <- unique(exp_design$Genotype)
stopifnot(length(source_genotypes) == 3)

double_line <- grep("/", source_genotypes, value = TRUE)
single_lines <- sort(setdiff(source_genotypes, double_line))
stopifnot(length(double_line) == 1, length(single_lines) == 2)

genotype_map <- setNames(c("Line_A", "Line_B", "Line_AB"),
                         c(single_lines, double_line))

source_treatments <- sort(unique(exp_design$Treatment))
stopifnot(length(source_treatments) == 2)

treatment_map <- setNames(c("Control", "Treated"), source_treatments)

exp_design$Genotype <- unname(genotype_map[exp_design$Genotype])
exp_design$Treatment <- unname(treatment_map[exp_design$Treatment])

# ---- Raw PSMs (evidence.txt) ---------------------------------------------

infdf <- read.delim(source_path("2025/2025_30_*/raw/evidence.txt.gz"))

# The complement of the split in data-raw/psm_tmt_factorial.R: the enriched
# fraction rather than the total.
infdf <- infdf[grepl("phos", infdf$Raw.file), ]

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

# As in data-raw/psm_tmt_factorial.R: the search database included custom
# entries for the tagged constructs, and the pseudo-accessions naming them are
# excluded entirely rather than relabelled.
construct_ids <- c("PER2MNG", "PER2HALO", "HALO001")
is_construct <- function(x) {
  Reduce(`|`, lapply(construct_ids, function(id) grepl(id, x, fixed = TRUE)))
}
infdf <- infdf[
  !is_construct(infdf$Proteins) &
    !is_construct(infdf$Leading.proteins) &
    !is_construct(infdf$Leading.razor.protein), ]

# ---- Subsample to a teaching-sized dataset -------------------------------

has_phospho <- !is.na(infdf$Phospho..STY..Probabilities) &
  infdf$Phospho..STY..Probabilities != ""

# 1) Proteins carrying more than one localisable phosphopeptide PSM. The
#    requirement is for the protein to have something to show at all, not for
#    it to show any particular thing: within that pool the draw is random, so
#    the proportion of sites that fail to localise is the dataset's own and not
#    an artefact of the selection.
phospho_proteins <- infdf %>%
  filter(has_phospho, Potential.contaminant == "", Reverse == "",
         !grepl(";", Leading.proteins)) %>%
  count(Leading.razor.protein) %>%
  filter(n >= 2)

keep_phospho_proteins <- sample(phospho_proteins$Leading.razor.protein,
                                min(700, nrow(phospho_proteins)))

# 2) A sample of contaminant PSMs, so the contaminant-filtering step has
#    something real to remove
contaminant_proteins <- unique(
  infdf$Leading.razor.protein[infdf$Potential.contaminant == "+"])
keep_contaminant_proteins <- sample(contaminant_proteins,
                                    min(20, length(contaminant_proteins)))

# 3) A slice of PSMs without a unique master protein
non_unique_ix <- which(grepl(";", infdf$Leading.proteins))
keep_non_unique_ix <- sample(non_unique_ix, min(200, length(non_unique_ix)))

# 4) A slice of decoy hits, so the decoy filter is doing real work too
decoy_ix <- which(infdf$Reverse != "")
keep_decoy_ix <- sample(decoy_ix, min(200, length(decoy_ix)))

psm_tmt_phospho_mq <- infdf %>%
  filter(
    Leading.razor.protein %in% c(keep_phospho_proteins,
                                 keep_contaminant_proteins) |
      row_number() %in% c(keep_non_unique_ix, keep_decoy_ix)
  ) %>%
  distinct()

# ---- Trim columns that cannot be used here --------------------------------

# `evidence.txt` cross-references MaxQuant's other output tables by row ID.
# None of those tables are packaged, so these columns are pointers to nothing.
dead_reference_cols <- c("MS.MS.IDs", "MS.MS.scan.numbers", "MS3.scan.numbers",
                         "MS.MS.scan.number", "Protein.group.IDs", "Peptide.ID",
                         "Mod..peptide.ID", "Best.MS.MS", "id",
                         "Oxidation..M..site.IDs", "Phospho..STY..site.IDs")

psm_tmt_phospho_mq <- psm_tmt_phospho_mq[
  , setdiff(colnames(psm_tmt_phospho_mq), dead_reference_cols)]

constant_cols <- names(psm_tmt_phospho_mq)[
  sapply(psm_tmt_phospho_mq, function(x) length(unique(x[!is.na(x)])) <= 1)]

cat("Dropping constant columns:", paste(constant_cols, collapse = ", "), "\n")

psm_tmt_phospho_mq <- psm_tmt_phospho_mq[
  , setdiff(colnames(psm_tmt_phospho_mq), constant_cols)]

# The localisation columns are the point of this dataset, so assert they
# survived the trimming rather than discovering their absence in the vignette.
stopifnot(all(c("Phospho..STY..Probabilities", "Phospho..STY.", "Sequence",
                "Modified.sequence", "Reverse", "Potential.contaminant",
                "Leading.razor.protein", "Proteins", "PIF") %in%
                colnames(psm_tmt_phospho_mq)))

# ---- Anonymise: Raw.file embeds the project ID, analyst and date ----------

source_runs <- sort(unique(psm_tmt_phospho_mq$Raw.file))
run_map <- setNames(sprintf("run_%02d", seq_along(source_runs)), source_runs)
psm_tmt_phospho_mq$Raw.file <- unname(run_map[psm_tmt_phospho_mq$Raw.file])

stopifnot(!any(is_construct(unlist(
  psm_tmt_phospho_mq[, c("Proteins", "Leading.proteins",
                         "Leading.razor.protein")]))))

cat(sprintf("psm_tmt_phospho_mq: %d PSMs, %d unique Leading.razor.protein\n",
            nrow(psm_tmt_phospho_mq),
            length(unique(psm_tmt_phospho_mq$Leading.razor.protein))))
cat(sprintf("  with localisation probabilities: %d\n",
            sum(psm_tmt_phospho_mq$Phospho..STY..Probabilities != "")))

usethis::use_data(psm_tmt_phospho_mq, overwrite = TRUE)

# ---- Proteome fasta for locating sites within proteins --------------------

# add_peptide_positions_from_cleavage() re-digests the protein sequences to
# place each peptide, and hence each site, within its protein.
#
# The sequences are taken from the search database rather than fetched with
# make_fasta(), which is what data-raw/psm_tmt_phospho.R does. Those are the
# sequences MaxQuant actually matched against, so every peptide places exactly;
# a current UniProt release would differ from the 2021 database this search
# used wherever a sequence has since been revised, and those peptides would
# silently fail to place. The construct entries are dropped with everything
# else that is not a retained master protein.
#
# The copy on disk is one revision later than the one `mqpar.xml` names, so a
# handful of accessions MaxQuant reported are absent from it. None of them
# survive PSM filtering, which is why the check below is made against the
# accessions that are actually placed rather than against every accession
# mentioned.
search_db <- source_path("2025/2025_30_*/raw/*/unip_human_*.fasta")

proteome <- Biostrings::readAAStringSet(search_db)
names(proteome) <- gsub("(sp|tr)\\|(\\S*)\\|.*", "\\2", names(proteome))

# MaxQuant prefixes its own contaminant entries with `CON__` and decoys with
# `REV__`, and neither is in the search database proper. Both are removed
# before sites are placed, so neither needs a sequence here.
fasta_accessions <- unique(psm_tmt_phospho_mq$Leading.razor.protein)
fasta_accessions <- sort(fasta_accessions[
  fasta_accessions != "" & !grepl(";", fasta_accessions) &
    !grepl("^(CON|REV)__", fasta_accessions)])

# Every master protein that survives PSM filtering must have a sequence, or
# its sites cannot be placed and the vignette loses them silently.
placeable <- unique(psm_tmt_phospho_mq$Leading.razor.protein[
  psm_tmt_phospho_mq$Potential.contaminant == "" &
    psm_tmt_phospho_mq$Reverse == "" &
    !grepl(";", psm_tmt_phospho_mq$Leading.proteins)])

missing <- setdiff(placeable, names(proteome))
if (length(missing) > 0) {
  stop(sprintf("%i placeable master proteins are absent from the search database: %s",
               length(missing), paste(head(missing), collapse = ", ")))
}

proteome <- proteome[intersect(fasta_accessions, names(proteome))]

fasta_out <- "inst/extdata/tmt_phospho_mq_proteome.fasta.gz"
Biostrings::writeXStringSet(proteome, fasta_out, compress = TRUE)

cat(sprintf("tmt_phospho_mq_proteome.fasta.gz: %d sequences\n",
            length(proteome)))
