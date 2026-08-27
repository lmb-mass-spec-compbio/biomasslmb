# Records inst/extdata/protein_annotation_cache.rds, the cached UniProt query
# results displayed by vignettes/protein_annotation.Rmd. Requires real network
# access to rest.uniprot.org.
#
# The vignette shows each query call but does not run it, reading the result
# from this cache instead, so that building the vignette does not depend on
# UniProt being reachable. Re-run this script to refresh the cache against a
# newer UniProt release; `release` records which release it was built from.
#
# Objects cached:
#   uniprot2details: get_uniprot_details() for every accession in the
#     lfq_dda_pd_PeptideGroups.txt master protein accessions. The 11 Cont_
#     prefixed contaminant accessions do not map and are absent from the
#     result -- this is deliberate, the vignette uses it to make the point.
#   uniparc_example: get_uniprot_details() for a stable accession and a
#     retired one, to show the Annotation.Source column taking both values.
#   go_res / go_res_all: get_go_terms() without and with ancestor expansion.
#   tm_topology: get_protein_tm_topology() for two 6-pass mitochondrial
#     carriers and one soluble protein.

devtools::load_all()

pep_inf <- system.file(
  "extdata", "lfq_dda_pd_PeptideGroups.txt",
  package = "biomasslmb"
)

protein_groups <- unique(read.delim(pep_inf)$Master.Protein.Accessions)

accessions <- protein_groups %>%
  strsplit('; ') %>%
  unlist() %>%
  unique()

protein_annotation_cache <- list(
  uniprot2details = get_uniprot_details(accessions),
  uniparc_example = get_uniprot_details(c('O76024', 'A0A5G2QPJ4')),
  go_res = get_go_terms(accessions),
  go_res_all = get_go_terms(accessions, expand_terms = TRUE),
  tm_topology = get_protein_tm_topology(c('P48962', 'P51881', 'P62827')),
  release = check_uniprot_release()
)

saveRDS(
  protein_annotation_cache,
  'inst/extdata/protein_annotation_cache.rds',
  compress = 'xz'
)
