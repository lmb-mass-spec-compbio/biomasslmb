# Records inst/extdata/lfq_dda_ibaq_PeptideGroups.txt.gz and
# inst/extdata/lfq_dda_ibaq_proteome.fasta.gz for the absolute quantification
# vignette.
#
# Source: 2026_05 project, the same Proteome Discoverer export as
# inst/extdata/lfq_dda_pd_PeptideGroups.txt. A second subset is needed rather
# than reusing that one because iBAQ has to be checked against a set of
# proteins whose relative molar abundance is known in advance, and a random
# few hundred proteins of a proteome contain almost none of any complex. Here
# every cytosolic ribosomal protein is retained on top of the random sample, so
# the vignette's check has something to check.
#
# The FASTA holds the sequences of the retained proteins. iBAQ divides summed
# intensity by the number of peptides a protein could yield, which means
# digesting the sequences in silico, and shipping them means the vignette does
# not query UniProt while building. Its headers are also where the vignette
# identifies the ribosomal proteins, since the PD peptide export carries no
# protein description column.
#
# Condition labels are generic, since the experiment is unpublished. See
# data-raw/source_paths.R for why the input paths are globbed rather than
# written out.

library(dplyr)

devtools::load_all()

set.seed(42)

source("data-raw/source_paths.R")

pep_inf <- source_path("2026/2026_05_*/raw/PD/*_PeptideGroups.txt")

infdf <- read.delim(pep_inf)

keep_cols <- c(
  "Confidence", "Sequence", "Annotated.Sequence", "Modifications",
  "Contaminant", "Number.of.Protein.Groups", "Number.of.Proteins",
  "Number.of.PSMs", "Master.Protein.Accessions", "Protein.Accessions",
  "Number.of.Missed.Cleavages",
  grep("^Abundance\\.F", colnames(infdf), value = TRUE),
  "PSM.Search.Engine.Rank.by.Search.Engine.CHIMERYS")

stopifnot(all(keep_cols %in% colnames(infdf)))
infdf <- infdf[, keep_cols]

# ---- Identify the ribosomal proteins -------------------------------------

# The peptide export has no description column, so the protein names come from
# UniProt. This is the only place the full accession list is needed.
all_accessions <- unique(unlist(strsplit(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "False"], "; ")))
all_accessions <- all_accessions[nzchar(all_accessions)]

cat(sprintf("querying UniProt for %d accessions\n", length(all_accessions)))

all_fasta <- tempfile(fileext = ".fasta")
make_fasta(all_accessions, file = all_fasta)
all_proteome <- Biostrings::readAAStringSet(all_fasta)

# UniProt names the cytosolic ribosome's subunits by family: the `e`
# (eukaryote-specific) and `u` (universal) prefixes are the cytosolic ribosome,
# while `m` is the mitochondrial ribosome, a separate complex at a quite
# different cellular abundance, so it is excluded. Older entries still carry
# the previous "40S/60S ribosomal protein" naming, so both forms are matched.
is_cytosolic_ribosomal <- grepl(
  "(Small|Large) ribosomal subunit protein [eu][LS][0-9]|^(40S|60S) ribosomal protein",
  sub("^\\S+\\s+", "", names(all_proteome))) &
  !grepl("mitochondrial", names(all_proteome))

accession_of <- gsub("(sp|tr)\\|(\\S*)\\|.*", "\\2", names(all_proteome))
ribosomal_proteins <- accession_of[is_cytosolic_ribosomal]

cat(sprintf("cytosolic ribosomal proteins in the search results: %d\n",
            length(ribosomal_proteins)))

# ---- Subsample to a teaching-sized dataset -------------------------------

unique_master_proteins <- unique(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "False" &
                                    infdf$Number.of.Protein.Groups == 1])
keep_proteins <- sample(unique_master_proteins, 400)
keep_proteins <- union(keep_proteins, ribosomal_proteins)

contaminant_proteins <- unique(
  infdf$Master.Protein.Accessions[infdf$Contaminant == "True"])
keep_contaminants <- sample(contaminant_proteins,
                            min(20, length(contaminant_proteins)))

non_unique_ix <- which(infdf$Number.of.Protein.Groups > 1)
keep_non_unique_ix <- sample(non_unique_ix, min(200, length(non_unique_ix)))

peptide_groups <- infdf %>%
  filter(Master.Protein.Accessions %in% c(keep_proteins, keep_contaminants) |
           row_number() %in% keep_non_unique_ix) %>%
  distinct()

cat(sprintf("lfq_dda_ibaq_PeptideGroups.txt.gz: %d peptides, %d master proteins\n",
            nrow(peptide_groups),
            length(unique(peptide_groups$Master.Protein.Accessions))))

out <- gzfile("inst/extdata/lfq_dda_ibaq_PeptideGroups.txt.gz", "w")
write.table(peptide_groups, out, sep = "\t", quote = FALSE, row.names = FALSE)
close(out)

# ---- Proteome fasta for the theoretical digest ----------------------------

retained <- unique(unlist(strsplit(
  peptide_groups$Master.Protein.Accessions, "; ")))
retained <- retained[nzchar(retained)]

proteome <- all_proteome[accession_of %in% retained]
Biostrings::writeXStringSet(proteome, "inst/extdata/lfq_dda_ibaq_proteome.fasta.gz",
                            compress = TRUE)
file.remove(all_fasta)

cat(sprintf("lfq_dda_ibaq_proteome.fasta.gz: %d sequences of %d accessions\n",
            length(proteome), length(retained)))
