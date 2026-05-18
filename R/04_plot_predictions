# ============================================================
# 04_plot_predictions.R
# Generate dosing prediction plots
# ============================================================

library(admr)
library(rxode2)
library(dplyr)
library(ggplot2)

fc_results <- readRDS("results/FC-results/fc_simulation.rds")
vo_results <- readRDS("results/VO-results/vo_simulation.rds")

time_points <- seq(0, 18, by = 0.1)

ev <- eventTable(amount.units = "mg", time.units = "hours")
ev$add.dosing(dose = 100, nbr.doses = 3, dosing.interval = 6)
ev$add.sampling(time_points)

plot_cols <- c(
  "True parameters" = "grey25",
  "Full Covariance" = "#4183c8",
  "Variance Only" = "#ff7f00"
)

# ============================================================
# Build true model
# ============================================================

build_true_model <- function(params) {
  
  rx_model <- function() {
    ini({
      cl <- params$beta["cl"]
      v1 <- params$beta["v1"]
      v2 <- params$beta["v2"]
      q  <- params$beta["q"]
      ka <- params$beta["ka"]
      
      eta_cl ~ params$Omega[1, 1]
      eta_v1 ~ params$Omega[2, 2]
      eta_v2 ~ params$Omega[3, 3]
      eta_q  ~ params$Omega[4, 4]
      eta_ka ~ params$Omega[5, 5]
    })
    
    model({
      cl <- cl * exp(eta_cl)
      v1 <- v1 * exp(eta_v1)
      v2 <- v2 * exp(eta_v2)
      q  <- q  * exp(eta_q)
      ka <- ka * exp(eta_ka)
      
      cp = linCmt(cl, v1, v2, q, ka)
    })
  }
  
  rxode2(rx_model())$simulationModel
}

# ============================================================
# Build misspecified model
# ============================================================

build_misspecified_model <- function(params, removed_iiv) {
  
  all_params <- c("cl", "v1", "v2", "q", "ka")
  random_params <- setdiff(all_params, removed_iiv)
  
  model_fun <- function() {
    ini({
      cl <- params$beta["cl"]
      v1 <- params$beta["v1"]
      v2 <- params$beta["v2"]
      q  <- params$beta["q"]
      ka <- params$beta["ka"]
      
      if ("cl" %in% random_params) eta_cl ~ params$Omega[which(random_params == "cl"), which(random_params == "cl")]
      if ("v1" %in% random_params) eta_v1 ~ params$Omega[which(random_params == "v1"), which(random_params == "v1")]
      if ("v2" %in% random_params) eta_v2 ~ params$Omega[which(random_params == "v2"), which(random_params == "v2")]
      if ("q"  %in% random_params) eta_q  ~ params$Omega[which(random_params == "q"),  which(random_params == "q")]
      if ("ka" %in% random_params) eta_ka ~ params$Omega[which(random_params == "ka"), which(random_params == "ka")]
    })
    
    model({
      if ("cl" %in% random_params) cl <- cl * exp(eta_cl)
      if ("v1" %in% random_params) v1 <- v1 * exp(eta_v1)
      if ("v2" %in% random_params) v2 <- v2 * exp(eta_v2)
      if ("q"  %in% random_params) q  <- q  * exp(eta_q)
      if ("ka" %in% random_params) ka <- ka * exp(eta_ka)
      
      cp = linCmt(cl, v1, v2, q, ka)
    })
  }
  
  rxode2(model_fun())$simulationModel
}

# ============================================================
# Simulate prediction intervals
# ============================================================

simulate_ci <- function(model, label) {
  
  sim <- rxSolve(
    model,
    events = ev,
    cores = 0,
    nSub = 10000
  )
  
  as.data.frame(confint(sim, "cp", level = 0.95)) %>%
    mutate(Model = label)
}

# ============================================================
# Plot function
# ============================================================

plot_ci <- function(ci_data, title) {
  
  ggplot(ci_data, aes(x = time)) +
    
    geom_ribbon(
      aes(ymin = p2.5, ymax = p97.5, fill = Model),
      alpha = 0.18
    ) +
    
    geom_line(
      aes(y = p50, colour = Model, linetype = Model),
      linewidth = 1
    ) +
    
    scale_colour_manual(values = plot_cols) +
    scale_fill_manual(values = plot_cols) +
    
    scale_linetype_manual(
      values = c(
        "True parameters" = "dashed",
        "Full Covariance" = "solid",
        "Variance Only" = "solid"
      )
    ) +
    
    labs(
      title = title,
      x = "Time (h)",
      y = "Central concentration (mg/L)"
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}

# ============================================================
# Generate plots
# ============================================================

for (par in removed_iiv_list) {
  
  message("Generating plots for removed IIV on: ", par)
  
  params_fc <- fc_results[[paste0("fit_", par)]]$transformed_params
  params_vo <- vo_results[[paste0("fit_var_", par)]]$transformed_params
  
  ci_true <- simulate_ci(
    build_true_model(params_true),
    "True parameters"
  )
  
  ci_fc <- simulate_ci(
    build_misspecified_model(params_fc, par),
    "Full Covariance"
  )
  
  ci_vo <- simulate_ci(
    build_misspecified_model(params_vo, par),
    "Variance Only"
  )
  
  p_compare <- plot_ci(
    bind_rows(ci_true, ci_fc, ci_vo),
    paste("Remove IIV on", toupper(par))
  )
  
  ggsave(
    filename = paste0("results/dosing-plots/", par, "_compare.png"),
    plot = p_compare,
    width = 10,
    height = 6,
    dpi = 300
  )
}
