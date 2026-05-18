# ============================================================
# 04_plot_predictions.R
# Generate dosing prediction plots
# ============================================================

library(admr)
library(rxode2)
library(dplyr)
library(ggplot2)
library(units)

fc_results <- readRDS("results/FC-results/fc_simulation.rds")
vo_results <- readRDS("results/VO-results/vo_simulation.rds")

time_points <- seq(0, 18, by = 0.1)

ev <- eventTable(amount.units = "mg", time.units = "hours")
ev$add.dosing(dose = 100, nbr.doses = 3, dosing.interval = 6)
ev$add.sampling(time_points)

plot_cols <- c(
  "True parameters" = "grey50",
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
  
  rx_model <- switch(
    removed_iiv,
    
    "cl" = function() {
      ini({
        cl <- params$beta["cl"]
        v1 <- params$beta["v1"]
        v2 <- params$beta["v2"]
        q  <- params$beta["q"]
        ka <- params$beta["ka"]
        
        eta_v1 ~ params$Omega[1, 1]
        eta_v2 ~ params$Omega[2, 2]
        eta_q  ~ params$Omega[3, 3]
        eta_ka ~ params$Omega[4, 4]
      })
      model({
        v1 <- v1 * exp(eta_v1)
        v2 <- v2 * exp(eta_v2)
        q  <- q  * exp(eta_q)
        ka <- ka * exp(eta_ka)
        
        cp = linCmt(cl, v1, v2, q, ka)
      })
    },
    
    "v1" = function() {
      ini({
        cl <- params$beta["cl"]
        v1 <- params$beta["v1"]
        v2 <- params$beta["v2"]
        q  <- params$beta["q"]
        ka <- params$beta["ka"]
        
        eta_cl ~ params$Omega[1, 1]
        eta_v2 ~ params$Omega[2, 2]
        eta_q  ~ params$Omega[3, 3]
        eta_ka ~ params$Omega[4, 4]
      })
      model({
        cl <- cl * exp(eta_cl)
        v2 <- v2 * exp(eta_v2)
        q  <- q  * exp(eta_q)
        ka <- ka * exp(eta_ka)
        
        cp = linCmt(cl, v1, v2, q, ka)
      })
    },
    
    "v2" = function() {
      ini({
        cl <- params$beta["cl"]
        v1 <- params$beta["v1"]
        v2 <- params$beta["v2"]
        q  <- params$beta["q"]
        ka <- params$beta["ka"]
        
        eta_cl ~ params$Omega[1, 1]
        eta_v1 ~ params$Omega[2, 2]
        eta_q  ~ params$Omega[3, 3]
        eta_ka ~ params$Omega[4, 4]
      })
      model({
        cl <- cl * exp(eta_cl)
        v1 <- v1 * exp(eta_v1)
        q  <- q  * exp(eta_q)
        ka <- ka * exp(eta_ka)
        
        cp = linCmt(cl, v1, v2, q, ka)
      })
    },
    
    "q" = function() {
      ini({
        cl <- params$beta["cl"]
        v1 <- params$beta["v1"]
        v2 <- params$beta["v2"]
        q  <- params$beta["q"]
        ka <- params$beta["ka"]
        
        eta_cl ~ params$Omega[1, 1]
        eta_v1 ~ params$Omega[2, 2]
        eta_v2 ~ params$Omega[3, 3]
        eta_ka ~ params$Omega[4, 4]
      })
      model({
        cl <- cl * exp(eta_cl)
        v1 <- v1 * exp(eta_v1)
        v2 <- v2 * exp(eta_v2)
        ka <- ka * exp(eta_ka)
        
        cp = linCmt(cl, v1, v2, q, ka)
      })
    },
    
    "ka" = function() {
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
      })
      model({
        cl <- cl * exp(eta_cl)
        v1 <- v1 * exp(eta_v1)
        v2 <- v2 * exp(eta_v2)
        q  <- q  * exp(eta_q)
        
        cp = linCmt(cl, v1, v2, q, ka)
      })
    }
  )
  
  rxode2(rx_model())$simulationModel
}

# ============================================================
# VPC-style simulation with CI around percentiles
# ============================================================

simulate_vpc_ci <- function(
    model,
    label,
    n_rep = 200,
    n_sub = 500,
    ci_level = 0.95
) {
  
  alpha <- (1 - ci_level) / 2
  
  one_rep <- function(i) {
    
    sim <- rxSolve(
      model,
      events = ev,
      cores = 0,
      nSub = n_sub
    )
    
    as.data.frame(sim) %>%
      mutate(
        time = as.numeric(time),
        cp = as.numeric(cp)
      ) %>%
      group_by(time) %>%
      summarise(
        q5  = quantile(cp, 0.05, na.rm = TRUE),
        q50 = quantile(cp, 0.50, na.rm = TRUE),
        q95 = quantile(cp, 0.95, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(rep = i)
  }
  
  vpc_raw <- purrr::map_dfr(1:n_rep, one_rep)
  
  vpc_ci <- vpc_raw %>%
    pivot_longer(
      cols = c(q5, q50, q95),
      names_to = "Percentile",
      values_to = "Value"
    ) %>%
    group_by(time, Percentile) %>%
    summarise(
      Line = median(Value, na.rm = TRUE),
      Lower = quantile(Value, alpha, na.rm = TRUE),
      Upper = quantile(Value, 1 - alpha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Model = label,
      Percentile = factor(
        Percentile,
        levels = c("q5", "q50", "q95"),
        labels = c("5th percentile", "50th percentile", "95th percentile")
      )
    )
  
  return(vpc_ci)
}

# ============================================================
# Plot VPC-style percentiles with CI bands
# ============================================================

plot_vpc_ci <- function(vpc_data, title) {
  
  ggplot(
    vpc_data,
    aes(
      x = time,
      y = Line,
      colour = Model,
      fill = Model,
      linetype = Model,
      group = interaction(Model, Percentile)
    )
  ) +
    
    # CI ribbons
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper),
      alpha = 0.25,
      colour = NA
    ) +
    
    geom_line(
      linewidth = 0.9
    ) +
    
    scale_colour_manual(values = plot_cols) +
    scale_fill_manual(values = plot_cols) +
    
    # Model line types
    scale_linetype_manual(values = c(
      "True parameters" = "dashed",
      "Full Covariance" = "solid",
      "Variance Only" = "solid"
    )) +
    
    labs(
      title = title,
      x = "Time (h)",
      y = "Central concentration (mg/L)",
      colour = "Model",
      fill = "Model",
      linetype = "Model"
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Generate plots
# ============================================================


for (par in removed_iiv_list) {
  
  message("Generating VPC plots for removed IIV on: ", par)
  
  params_fc <- fc_results[[paste0("fit_", par)]]$transformed_params
  params_vo <- vo_results[[paste0("fit_var_", par)]]$transformed_params
  
  vpc_true <- simulate_vpc_ci(
    model = build_true_model(params_true),
    label = "True parameters",
    n_rep = 200,
    n_sub = 500
  )
  
  vpc_fc <- simulate_vpc_ci(
    model = build_misspecified_model(params_fc, par),
    label = "Full Covariance",
    n_rep = 200,
    n_sub = 500
  )
  
  vpc_vo <- simulate_vpc_ci(
    model = build_misspecified_model(params_vo, par),
    label = "Variance Only",
    n_rep = 200,
    n_sub = 500
  )
  
  # Combined comparison plot
  p_compare <- plot_vpc_ci(
    bind_rows(vpc_true, vpc_fc, vpc_vo),
    paste("Remove IIV on", toupper(par))
  )
  
  ggsave(
    filename = paste0("results/dosing-plots/", par, "_vpc_compare.png"),
    plot = p_compare,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  # FC separate plot
  p_fc <- plot_vpc_ci(
    bind_rows(vpc_true, vpc_fc),
    paste("Remove IIV on", toupper(par), "| Full Covariance")
  )
  
  ggsave(
    filename = paste0("results/dosing-plots/", par, "_vpc_fc.png"),
    plot = p_fc,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # VO separate plot
  p_vo <- plot_vpc_ci(
    bind_rows(vpc_true, vpc_vo),
    paste("Remove IIV on", toupper(par), "| Variance Only")
  )
  
  ggsave(
    filename = paste0("results/dosing-plots/", par, "_vpc_vo.png"),
    plot = p_vo,
    width = 10,
    height = 6,
    dpi = 300
  )
}
