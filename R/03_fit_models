# ============================================================
# 03_fit_models.R
# Fit FC and VO models under IIV misspecification
# ============================================================

library(purrr)

removed_iiv_list <- c("cl", "v1", "v2", "q", "ka")

# Full covariance options
fc_opts <- create_base_opts(no_cov = FALSE)

# Variance-only options
vo_opts <- create_base_opts(no_cov = TRUE)

# ============================================================
# Fit full covariance models
# ============================================================

fc_results <- map(
  removed_iiv_list,
  ~ fit_removed_iiv_model(
    base_opts = fc_opts,
    obs_data = examplomycin_aggregated,
    removed_iiv = .x
  )
)

names(fc_results) <- paste0("fit_", removed_iiv_list)

saveRDS(
  fc_results,
  file = "results/FC-results/fc_simulation.rds"
)

# ============================================================
# Fit variance-only models
# ============================================================

vo_results <- map(
  removed_iiv_list,
  ~ fit_removed_iiv_model(
    base_opts = vo_opts,
    obs_data = examplomycin_agg_var,
    removed_iiv = .x
  )
)

names(vo_results) <- paste0("fit_var_", removed_iiv_list)

saveRDS(
  vo_results,
  file = "results/VO-results/vo_simulation.rds"
)
