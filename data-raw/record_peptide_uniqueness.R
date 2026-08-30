# Records inst/extdata/peptide_uniqueness.rds and
# inst/extdata/C_vibrioides_proteome.fasta.gz for the peptide-to-protein
# vignette.
#
# The vignette's argument is a comparison between a large proteome, where a
# tryptic peptide often maps to several proteins, and a small one, where it
# rarely does. Running that comparison live would mean shipping a human
# reference proteome (7.4 MB compressed) and digesting 20,000 proteins while
# the vignette builds, so the human side is precomputed here and cached. The
# Caulobacter proteome is small enough to ship and digest live, so the vignette
# can show the code actually running and the reader can rerun it.
#
# Both inputs are public UniProt reference proteomes:
#   UP000005640 Homo sapiens
#   UP000001364 Caulobacter vibrioides
#
# Cached, per proteome:
#   peptides: for every tryptic peptide of observable length, the number of
#     proteins in the proteome it could have come from
#   proteins: for every protein, how many of its peptides are unique to it and
#     what fraction of its peptides those are
#   n_proteins: the size of the proteome

library(dplyr)
library(tidyr)

source("data-raw/source_paths.R")

workshop_root <- "~/git_repos/training/Proteomics_workshop_data_exploration"

human_fasta <- source_path("raw/UP000005640_*.fasta.gz", root = workshop_root)
cvib_fasta <- source_path("raw/UP000001364_*.fasta.gz", root = workshop_root)

# Peptides shorter than 6 or longer than 30 residues are rarely identified, so
# counting them would overstate how much of the proteome is actually at stake.
MIN_LENGTH <- 6
MAX_LENGTH <- 30

digest_proteome <- function(fasta_filepath) {
  proteins <- Biostrings::readAAStringSet(fasta_filepath)
  names(proteins) <- gsub('(sp|tr)\\|(\\S*)\\|.*', '\\2', names(proteins))

  peptides <- cleaver::cleave(proteins, enzym = "trypsin", unique = FALSE)

  pep2prot <- data.frame(
    protein = rep(names(peptides), lengths(peptides)),
    peptide = unlist(peptides, use.names = FALSE)) %>%
    filter(nchar(peptide) >= MIN_LENGTH, nchar(peptide) <= MAX_LENGTH) %>%
    distinct()

  per_peptide <- pep2prot %>%
    count(peptide, name = "n_proteins")

  per_protein <- pep2prot %>%
    left_join(per_peptide, by = "peptide") %>%
    group_by(protein) %>%
    summarise(n_peptides = n(),
              n_unique = sum(n_proteins == 1),
              perc_unique = 100 * mean(n_proteins == 1),
              .groups = "drop")

  list(
    # only the distribution is needed, not the peptide sequences
    peptides = per_peptide %>% count(n_proteins, name = "n_peptides"),
    proteins = per_protein,
    n_proteins = length(proteins))
}

peptide_uniqueness <- list(
  `Homo sapiens` = digest_proteome(human_fasta),
  `Caulobacter vibrioides` = digest_proteome(cvib_fasta))

for (species in names(peptide_uniqueness)) {
  res <- peptide_uniqueness[[species]]
  cat(sprintf("%s: %d proteins, %d distinct peptides, %.1f%% peptide-unique, median %d unique peptides per protein\n",
              species, res$n_proteins, sum(res$peptides$n_peptides),
              100 * sum(res$peptides$n_peptides[res$peptides$n_proteins == 1]) /
                sum(res$peptides$n_peptides),
              median(res$proteins$n_unique)))
}

saveRDS(peptide_uniqueness, "inst/extdata/peptide_uniqueness.rds")

file.copy(cvib_fasta, "inst/extdata/C_vibrioides_proteome.fasta.gz",
          overwrite = TRUE)

cat(sprintf("peptide_uniqueness.rds: %.0f KB\n",
            file.size("inst/extdata/peptide_uniqueness.rds") / 1024))
