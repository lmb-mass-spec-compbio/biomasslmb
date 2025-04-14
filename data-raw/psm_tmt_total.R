# Load libraries
library(camprotR)
library(dplyr)
# psm_tmt_total is part of the camprotR package which we will reuse here for now

# Output .rda file
psm_tmt_total <- psm_tmt_total %>%
  mutate(Quan.Info=ifelse(is.na(Quan.Info), '', Quan.Info))
 usethis::use_data(psm_tmt_total, overwrite = TRUE)

