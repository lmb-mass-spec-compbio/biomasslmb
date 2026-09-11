# Records inst/extdata/lfq_dda_mq_peptides.txt.gz, the MaxQuant peptide-level
# output for the same six runs as inst/extdata/lfq_dda_pd_PeptideGroups.txt,
# used by the search engine comparison vignette.
#
# Source: 2026_05 project, the same acquisition as the Proteome Discoverer
# export recorded by data-raw/record_lfq_dda_peptide_groups.R. Both files
# therefore describe identical raw data, and any difference between them is a
# difference between the two search and quantification pipelines rather than
# between experiments.
#
# The two subsets are matched on protein: this file is restricted to the master
# proteins already present in the PD subset (plus contaminants and
# non-unique-razor-protein peptides, so those filtering steps have something
# real to remove). That makes the pair suitable for comparing quantification
# and testing, and unsuitable for comparing identification depth, since neither
# engine was given the chance to contribute proteins the other's subset lacks.
#
# The intensity columns keep the MaxQuant experiment names (Intensity.A ...),
# because relabelling them from a design table is the first thing the vignette
# does. In mqpar.xml the experiments are listed in acquisition order against
# raw files numbered 1-6, which are samples 1-6 of raw/Sample_Info.xlsx, i.e.
# the three wild type replicates followed by the three mutant replicates.
# Condition labels are generic, since the experiment is unpublished. See
# data-raw/source_paths.R for why the input paths are globbed rather than
# written out.

library(dplyr)

set.seed(42)

source("data-raw/source_paths.R")

pep_inf <- source_path("2026/2026_05_*/raw/peptides.txt")
mqpar_inf <- source_path("2026/2026_05_*/raw/mqpar.xml")
design_inf <- source_path("2026/2026_05_*/raw/Sample_Info.xlsx")

# ---- Check the experiments line up with the design -----------------------

mqpar <- xml2::read_xml(mqpar_inf)

raw_files <- basename(gsub("\\\\", "/",
                           xml2::xml_text(xml2::xml_find_all(mqpar, ".//filePaths/string"))))
experiments <- xml2::xml_text(xml2::xml_find_all(mqpar, ".//experiments/string"))

# The raw files are numbered for the sample they hold
sample_numbers <- as.numeric(sub(".*_(\\d+)\\.raw$", "\\1", raw_files))

exp_design <- readxl::read_excel(design_inf)

# The two conditions sort with the mutant before the wild type, a property of
# this design sheet, so assert the shape rather than trusting it silently.
source_conditions <- sort(unique(exp_design$condition))
stopifnot(length(source_conditions) == 2,
          identical(sort(sample_numbers), sort(exp_design$Sample)))

condition_map <- setNames(c("Mutant", "WT"), source_conditions)

sample_order <- exp_design %>%
  mutate(condition = unname(condition_map[condition])) %>%
  arrange(desc(condition), replicate) %>%
  pull(Sample)

# The experiment letters, in the order the vignette's design table expects
ordered_experiments <- experiments[match(sample_order, sample_numbers)]

# ---- Columns to retain ---------------------------------------------------

infdf <- read.delim(pep_inf)

intensity_cols <- paste0("Intensity.", ordered_experiments)

keep_cols <- c(
  "Sequence", "Length", "Missed.cleavages", "Mass", "Proteins",
  "Leading.razor.protein", "Start.position", "End.position", "Gene.names",
  "Protein.names", "Unique..Groups.", "Unique..Proteins.", "Charges", "PEP",
  "Score", intensity_cols, "Reverse", "Potential.contaminant", "MS.MS.Count")

stopifnot(all(keep_cols %in% colnames(infdf)))

infdf <- infdf[, keep_cols]

# ---- Subsample to match the PD export ------------------------------------

pd_proteins <- unique(read.delim(
  "inst/extdata/lfq_dda_pd_PeptideGroups.txt")$Master.Protein.Accessions)

contaminant_proteins <- unique(
  infdf$Leading.razor.protein[infdf$Potential.contaminant == "+"])
keep_contaminants <- sample(contaminant_proteins,
                            min(20, length(contaminant_proteins)))

non_unique_ix <- which(infdf$Unique..Proteins. != "yes")
keep_non_unique_ix <- sample(non_unique_ix, min(300, length(non_unique_ix)))

mq_peptides <- infdf %>%
  filter(Leading.razor.protein %in% c(pd_proteins, keep_contaminants) |
           row_number() %in% keep_non_unique_ix) %>%
  distinct()

cat(sprintf("lfq_dda_mq_peptides.txt.gz: %d peptides, %d razor proteins (%d of the %d PD master proteins)\n",
            nrow(mq_peptides),
            length(unique(mq_peptides$Leading.razor.protein)),
            length(intersect(mq_peptides$Leading.razor.protein, pd_proteins)),
            length(pd_proteins)))

out <- gzfile("inst/extdata/lfq_dda_mq_peptides.txt.gz", "w")
write.table(mq_peptides, out, sep = "\t", quote = FALSE, row.names = FALSE)
close(out)
