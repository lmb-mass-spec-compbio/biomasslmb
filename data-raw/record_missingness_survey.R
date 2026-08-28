# Records inst/extdata/missingness_survey.rds, the cross-dataset summary of
# missingness against abundance displayed by vignettes/handling_missing_values.Rmd.
#
# The source is a survey of 56 real experiments processed by the facility,
# assembled in the Proteomics_workshop_data_exploration repository. Each
# experiment contributes one feature-level abundance and one missingness value
# per feature, plus three dataset-level scores. Only the summary is recorded
# here: the underlying experiments are unpublished, and the vignette needs the
# shape of the relationship rather than the data behind it.
#
# What is deliberately not carried over:
#   infile, assay_name, groups -- the file paths name collaborators and the
#     group columns name baits. Nothing in the vignette needs them.
# `organism` is normalised: the source has both 'M. musculus' and 'M. Musculus'.
#
# Features are subsampled to `n_per_exp` per experiment so that the shipped file
# stays small and no single large experiment dominates a density plot. The
# source object is ~49 MB.
#
# Columns recorded:
#   features: exp, mean_log_intensity (min-max scaled within experiment),
#     perc_missing (fraction of samples with no value for that feature)
#   experiments: exp, ms_type, sample_type, organism,
#     tjur_r2_intensity_only (how well abundance alone predicts missingness),
#     condition_miss_index (coverage-penalised), condition_miss_weighted_mean

set.seed(42)

n_per_exp <- 1200

workshop <- '~/git_repos/training/Proteomics_workshop_data_exploration/notebooks'

features_all <- readRDS(file.path(workshop, 'mean_vs_perc_miss_combined.rds'))
attrs <- readRDS(file.path(workshop, 'attr_summary.rds'))

features <- features_all %>%
  dplyr::select(exp, mean_log_intensity, perc_missing) %>%
  dplyr::filter(is.finite(mean_log_intensity), is.finite(perc_missing)) %>%
  dplyr::group_by(exp) %>%
  dplyr::slice_sample(n = n_per_exp) %>%
  dplyr::ungroup() %>%
  as.data.frame()

experiments <- attrs %>%
  dplyr::transmute(
    exp,
    ms_type,
    sample_type,
    organism = sub('M. Musculus', 'M. musculus', organism, fixed = TRUE),
    tjur_r2_intensity_only,
    condition_miss_index = mnar_i_cp,
    condition_miss_weighted_mean = mnar_i_no_cp
  )

missingness_survey <- list(
  features = features,
  experiments = experiments,
  n_per_exp = n_per_exp
)

saveRDS(
  missingness_survey,
  'inst/extdata/missingness_survey.rds',
  compress = 'xz'
)
