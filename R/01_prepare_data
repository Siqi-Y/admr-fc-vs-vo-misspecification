# ============================================================
# 01_prepare_data.R
# Prepare aggregate data
# ============================================================

library(admr)
library(dplyr)
library(tidyr)

data("examplomycin")

# Convert individual-level long data to subject-by-time wide format
examplomycin_wide <- examplomycin %>%
  filter(EVID != 101) %>%
  select(ID, TIME, DV) %>%
  pivot_wider(names_from = TIME, values_from = DV) %>%
  select(-ID)

# Full covariance aggregate data
examplomycin_aggregated <- examplomycin_wide %>%
  admr::meancov()

# Variance-only aggregate data
examplomycin_agg_var <- examplomycin_aggregated
examplomycin_agg_var$V <- diag(diag(examplomycin_aggregated$V))

# Create result folders
dir.create("results", showWarnings = FALSE)
dir.create("results/FC-results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/VO-results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/dosing-plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
