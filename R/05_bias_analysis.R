# ============================================================
# 05_bias_analysis.R
# Extract parameter bias and BSV bias
# ============================================================

library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

fc_results <- readRDS("results/FC-results/fc_simulation.rds")
vo_results <- readRDS("results/VO-results/vo_simulation.rds")

true_theta <- c(
  cl = 5,
  v1 = 10,
  v2 = 30,
  q = 10,
  ka = 1
)

true_bsv <- 30.4

# ============================================================
# Extract results from admr output
# ============================================================

extract_admr_results <- function(fit_list, method_label) {
  
  imap_dfr(fit_list, function(fit, fit_name) {
    
    removed <- toupper(sub("fit_(var_)?", "", fit_name))
    
    fit$param_df %>%
      filter(Parameter != "Residual Error") %>%
      mutate(
        Est_bt = as.numeric(str_extract(
          `Back-transformed(95%CI)`,
          "^[0-9.]+"
        )),
        
        CI_low = as.numeric(str_extract(
          `Back-transformed(95%CI)`,
          "(?<=\\()[0-9.]+"
        )),
        
        CI_high = as.numeric(str_extract(
          `Back-transformed(95%CI)`,
          "(?<=, )[0-9.]+"
        )),
        
        True_Val = true_theta[Parameter],
        Rel_Bias = (Est_bt - True_Val) / True_Val * 100,
        
        True_BSV = true_bsv,
        BSV_Estimated = as.numeric(`BSV(CV%)`),
        BSV_Bias = (BSV_Estimated - True_BSV) / True_BSV * 100,
        
        Method = method_label,
        Removed_IIV = removed
      )
  })
}

comparison_df <- bind_rows(
  extract_admr_results(fc_results, "Full Covariance"),
  extract_admr_results(vo_results, "Variance Only")
)

saveRDS(
  comparison_df,
  file = "results/tables/comparison_df.rds"
)

write.csv(
  comparison_df,
  file = "results/tables/comparison_df.csv",
  row.names = FALSE
)

# ============================================================
# Parameter bias heatmap
# ============================================================

p_param_heatmap <- ggplot(
  comparison_df,
  aes(x = Parameter, y = Removed_IIV, fill = Rel_Bias)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Rel_Bias, 1)), size = 3) +
  facet_wrap(~ Method) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0
  ) +
  labs(
    title = "Relative Bias of Parameter Estimates",
    x = "Parameter",
    y = "Removed IIV",
    fill = "Bias (%)"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "results/tables/parameter_bias_heatmap.png",
  p_param_heatmap,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# BSV bias heatmap
# ============================================================

p_bsv_heatmap <- ggplot(
  comparison_df,
  aes(x = Parameter, y = Removed_IIV, fill = BSV_Bias)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = ifelse(is.na(BSV_Bias), "", round(BSV_Bias, 1))),
    size = 3
  ) +
  facet_wrap(~ Method) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    na.value = "white"
  ) +
  labs(
    title = "Relative Bias of BSV",
    x = "Parameter",
    y = "Removed IIV",
    fill = "Bias (%)"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "results/tables/bsv_bias_heatmap.png",
  p_bsv_heatmap,
  width = 10,
  height = 6,
  dpi = 300
)
