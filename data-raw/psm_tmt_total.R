# Load libraries
library(camprotR)
library(dplyr)
# psm_tmt_total is part of the camprotR package which we will reuse here for now

# Spectrum.File holds the original acquisition filenames, which name a
# collaborator. Each maps 1:1 to a File.ID, so the files are renamed to
# run_NN.raw in acquisition order, matching the convention used by the other
# shipped example datasets.
file_ids <- unique(psm_tmt_total$File.ID)
file_ids <- file_ids[order(as.numeric(sub('^F\\d+\\.', '', file_ids)))]
run_names <- setNames(sprintf('run_%02i.raw', seq_along(file_ids)), file_ids)

# Output .rda file
psm_tmt_total <- psm_tmt_total %>%
  mutate(Quan.Info=ifelse(is.na(Quan.Info), '', Quan.Info)) %>%
  mutate(Spectrum.File=unname(run_names[File.ID]))
 usethis::use_data(psm_tmt_total, overwrite = TRUE)
