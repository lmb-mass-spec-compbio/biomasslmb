# Records inst/extdata/spectronaut_report.tsv.gz, the Spectronaut precursor
# report read by vignettes/LFQ_DIA_Precursor_QC_Summarisation.Rmd.
#
# Source: the total-lysate arm of project 2026_17. Two of its conditions are
# kept, three biological replicates each, since the vignette uses this file
# only to show how a Spectronaut export is read into a QFeatures object rather
# than to work through a second analysis.
#
# The export is a 'BGS Factory Report (Normal)' already reduced to the 16
# columns needed downstream. It is 91 MB compressed at full size, so protein
# groups are subsampled and float precision is reduced here.
#
# Runs are renamed to Control_1..3 / Treated_1..3, per the usual rule for
# shipped example data. See data-raw/source_paths.R for why the input path is
# globbed rather than written out.

source("data-raw/source_paths.R")

set.seed(42)

n_protein_groups <- 250

infile <- source_path("2026/2026_17_*/raw/lysate/*_filtered.tsv.gz")

# Sample numbers taken from the project's sample designation sheet.
run_map <- c(
  '7'  = 'Control_1',
  '17' = 'Control_2',
  '27' = 'Control_3',
  '8'  = 'Treated_1',
  '18' = 'Treated_2',
  '28' = 'Treated_3'
)

precursor_df <- read.delim(infile)

sample_number <- sub('.*_sample', '', precursor_df$R.FileName)

precursor_df <- precursor_df[sample_number %in% names(run_map), ]

precursor_df$R.FileName <- unname(
  run_map[sub('.*_sample', '', precursor_df$R.FileName)])

keep_groups <- sample(unique(precursor_df$PG.ProteinGroups), n_protein_groups)

precursor_df <- precursor_df[precursor_df$PG.ProteinGroups %in% keep_groups, ]

precursor_df <- precursor_df[order(precursor_df$R.FileName,
                                   precursor_df$PG.ProteinGroups), ]

# Spectronaut writes 15 significant figures. Rounding costs nothing downstream
# and roughly halves the compressed file.
for (col in c('PG.Qvalue', 'PG.PEP', 'EG.Qvalue', 'EG.PEP')) {
  precursor_df[[col]] <- signif(precursor_df[[col]], 4)
}
for (col in c('EG.TotalQuantity..Settings.', 'FG.Quantity')) {
  precursor_df[[col]] <- signif(precursor_df[[col]], 6)
}

outfile <- 'inst/extdata/spectronaut_report.tsv'

write.table(precursor_df, file = outfile,
            sep = '\t', row.names = FALSE, quote = FALSE)

R.utils::gzip(outfile, overwrite = TRUE)
