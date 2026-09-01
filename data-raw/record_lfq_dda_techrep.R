# Records inst/extdata/lfq_dda_techrep_peptides.txt.gz, a MaxQuant peptide-level
# export with technical replicates, for the technical replicates vignette.
#
# No other example dataset in the package has repeated measurements of the same
# biological sample, so none of them can show what technical replicates do or
# what goes wrong when they are treated as biological ones.
#
# Source: project 2026_09, a bait pulldown against a tag-only control, three
# biological replicates, each acquired twice. The design is not recorded in a
# spreadsheet here; it is encoded in the raw file names, which also carry the
# analyst's name and the bait, so it is parsed structurally and only the
# structure is written out. Condition labels are generic, matching the
# convention used for the other enrichment dataset in the package. See
# data-raw/source_paths.R for why the input paths are globbed rather than
# written out.

library(dplyr)

set.seed(42)

source("data-raw/source_paths.R")

# The project holds two acquisitions and only one of them was run with
# technical replicates. Rather than naming it, each is parsed and the one whose
# design has every biological replicate acquired twice is selected. The
# informative part of a raw file name is its suffix: a biological replicate
# number, a letter for the condition, and an optional trailing "r" marking the
# repeat acquisition. Everything before that is discarded unparsed.
parse_design <- function(mqpar_inf) {
  mqpar <- xml2::read_xml(mqpar_inf)
  raw_files <- basename(gsub("\\\\", "/", xml2::xml_text(
    xml2::xml_find_all(mqpar, ".//filePaths/string"))))
  experiments <- xml2::xml_text(
    xml2::xml_find_all(mqpar, ".//experiments/string"))

  suffix <- sub(".*_(\\d+[A-Z]r?)(_\\d+)?\\.raw$", "\\1", raw_files)
  if (any(suffix == raw_files)) return(NULL)

  data.frame(
    experiments,
    Replicate = as.integer(sub("^(\\d+).*", "\\1", suffix)),
    condition_letter = sub("^\\d+([A-Z]).*", "\\1", suffix),
    Tech_rep = ifelse(grepl("r$", suffix), 2L, 1L))
}

has_technical_replicates <- function(design) {
  !is.null(design) &&
    all(table(design$Replicate, design$condition_letter) == 2)
}

mqpar_files <- Sys.glob(path.expand(file.path(
  projects_root, "2026/2026_09_*/raw/*/mqpar.xml")))
designs <- lapply(mqpar_files, parse_design)
selected <- which(vapply(designs, has_technical_replicates, logical(1)))

stopifnot(length(selected) == 1)

design <- designs[[selected]]
pep_inf <- file.path(dirname(mqpar_files[selected]), "peptides.txt")
stopifnot(file.exists(pep_inf))

# Two conditions, the control sorting first by letter. Assert the shape rather
# than trusting it silently.
source_conditions <- sort(unique(design$condition_letter))
stopifnot(length(source_conditions) == 2, nrow(design) == 12)

condition_map <- setNames(c("Control", "IP"), source_conditions)

design <- design %>%
  mutate(Condition = unname(condition_map[condition_letter]),
         Sample = paste(Condition, Replicate, Tech_rep, sep = "_")) %>%
  arrange(Condition, Replicate, Tech_rep)

intensity_cols <- paste0("Intensity.", design$experiments)

# ---- Columns to retain ---------------------------------------------------

infdf <- read.delim(pep_inf)

keep_cols <- c(
  "Sequence", "Length", "Missed.cleavages", "Mass", "Proteins",
  "Leading.razor.protein", "Unique..Groups.",
  "Unique..Proteins.", "PEP", "Score", intensity_cols,
  "Reverse", "Potential.contaminant", "MS.MS.Count")

stopifnot(all(keep_cols %in% colnames(infdf)))
infdf <- infdf[, keep_cols]

# The intensity columns are named for the samples, since the mapping from
# MaxQuant's experiment letters lives only in the raw file names
colnames(infdf)[match(intensity_cols, colnames(infdf))] <- design$Sample

# ---- Subsample to a teaching-sized dataset -------------------------------

unique_razor <- unique(
  infdf$Leading.razor.protein[infdf$Potential.contaminant != "+" &
                                infdf$Reverse != "+" &
                                infdf$Unique..Proteins. == "yes"])
keep_proteins <- sample(unique_razor, min(400, length(unique_razor)))

contaminant_proteins <- unique(
  infdf$Leading.razor.protein[infdf$Potential.contaminant == "+"])
keep_contaminants <- sample(contaminant_proteins,
                            min(20, length(contaminant_proteins)))

non_unique_ix <- which(infdf$Unique..Proteins. != "yes")
keep_non_unique_ix <- sample(non_unique_ix, min(200, length(non_unique_ix)))

peptides <- infdf %>%
  filter(Leading.razor.protein %in% c(keep_proteins, keep_contaminants) |
           row_number() %in% keep_non_unique_ix) %>%
  distinct()

cat(sprintf("lfq_dda_techrep_peptides.txt.gz: %d peptides, %d razor proteins\n",
            nrow(peptides), length(unique(peptides$Leading.razor.protein))))
print(design[, c("Sample", "Condition", "Replicate", "Tech_rep")])

out <- gzfile("inst/extdata/lfq_dda_techrep_peptides.txt.gz", "w")
write.table(peptides, out, sep = "\t", quote = FALSE, row.names = FALSE)
close(out)
