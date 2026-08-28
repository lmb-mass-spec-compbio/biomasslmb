# Records inst/extdata/go_enrichment_cache.rds, the GO term annotations used by
# vignettes/functional_enrichment.Rmd. Requires real network access to
# rest.uniprot.org.
#
# The vignette shows the get_go_terms() call but does not run it, reading the
# result from this cache instead, so that building the vignette does not depend
# on UniProt being reachable. Ancestor expansion is the slow part -- it takes
# several minutes for this many proteins -- which is the other reason for
# caching it. Re-run this script to refresh against a newer UniProt release;
# `release` records which release it was built from.
#
# Objects cached:
#   go_res_all: get_go_terms(expand_terms = TRUE) for every accession in the
#     dia_qf protein assay, i.e. the background set the vignette tests against.

devtools::load_all()

accessions <- rownames(dia_qf[['protein']])

go_enrichment_cache <- list(
  go_res_all = get_go_terms(accessions, expand_terms = TRUE),
  release = check_uniprot_release()
)

saveRDS(
  go_enrichment_cache,
  'inst/extdata/go_enrichment_cache.rds',
  compress = 'xz'
)
