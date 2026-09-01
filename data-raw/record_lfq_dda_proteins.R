# Records inst/extdata/lfq_dda_pd_Proteins.txt.gz, the Proteome Discoverer
# protein-level export matching inst/extdata/lfq_dda_pd_PeptideGroups.txt.
#
# Only needed for filter_by_protein_fdr(), which reads the protein-level FDR
# confidence out of a separate file rather than from the peptide export. Just
# the columns that function uses, plus the description, are kept.
#
# Source: 2026_05 project, the same search as the peptide export recorded by
# data-raw/record_lfq_dda_peptide_groups.R. See data-raw/source_paths.R for why
# the input paths are globbed rather than written out.

library(dplyr)

source("data-raw/source_paths.R")

prot_inf <- source_path("2026/2026_05_*/raw/PD/*_Proteins.txt")

infdf <- read.delim(prot_inf)

keep_cols <- c("Accession", "Description", "Contaminant",
               "Protein.FDR.Confidence.Combined",
               "Protein.Group.FDR.Confidence.Combined",
               "Number.of.Peptides", "Number.of.Unique.Peptides")

stopifnot(all(keep_cols %in% colnames(infdf)))

# The peptide export was subset to a sample of proteins; keep the protein rows
# for those, so the two files describe the same set.
peptide_groups <- read.delim("inst/extdata/lfq_dda_pd_PeptideGroups.txt")
peptide_proteins <- unique(unlist(strsplit(
  peptide_groups$Master.Protein.Accessions, "; ")))

proteins <- infdf[infdf$Accession %in% peptide_proteins, keep_cols]

cat(sprintf("lfq_dda_pd_Proteins.txt.gz: %d proteins\n", nrow(proteins)))
print(table(proteins$Protein.FDR.Confidence.Combined))

out <- gzfile("inst/extdata/lfq_dda_pd_Proteins.txt.gz", "w")
write.table(proteins, out, sep = "\t", quote = FALSE, row.names = FALSE)
close(out)
