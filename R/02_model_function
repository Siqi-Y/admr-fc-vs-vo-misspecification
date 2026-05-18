# ============================================================
# 02_model_functions.R
# Model functions and helper functions
# ============================================================

library(admr)
library(rxode2)
library(mnorm)
library(dplyr)
library(purrr)

rxode2::rxSetSilentErr(1)

# Observation times used in aggregate data
obs_times <- c(0.1, 0.25, 0.5, 1, 2, 3, 5, 8, 12)

# True parameter values
params_true <- list(
  beta = c(cl = 5, v1 = 10, v2 = 30, q = 10, ka = 1),
  Omega = diag(rep(0.09, 5)),
  Sigma_prop = 0.04
)

# Initial parameter values for fitting
initial_params <- list(
  beta = c(cl = 4, v1 = 12, v2 = 25, q = 12, ka = 1.2),
  Omega = omegas(0.09, 0, 4),
  Sigma_prop = 0.04
)

# ============================================================
# Base rxode2 model
# ============================================================

rx_model_base <- function() {
  model({
    cp = linCmt(
      cl,
      v1,
      v2,
      q,
      ka
    )
  })
}

rx_model <- rxode2(rx_model_base)
rx_model <- rx_model$simulationModel

# ============================================================
# Prediction function for admr
# ============================================================

predder <- function(time, theta_i, dose = 100) {
  
  n_individuals <- nrow(theta_i)
  
  if (is.null(n_individuals)) {
    n_individuals <- 1
  }
  
  ev <- eventTable(amount.units = "mg", time.units = "hours")
  ev$add.dosing(dose = dose, nbr.doses = 1, start.time = 0)
  ev$add.sampling(time)
  
  out <- rxSolve(
    rx_model,
    params = theta_i,
    events = ev,
    cores = 0
  )
  
  cp_matrix <- matrix(
    out$cp,
    nrow = n_individuals,
    ncol = length(time),
    byrow = TRUE
  )
  
  return(cp_matrix)
}

# ============================================================
# Create admr options
# ============================================================

create_base_opts <- function(no_cov = FALSE) {
  
  list(
    time = obs_times,
    p = initial_params,
    nsim = 10000,
    n = 500,
    fo_appr = FALSE,
    omega_expansion = 1,
    f = predder,
    no_cov = no_cov
  )
}

# ============================================================
# Fit model with one IIV term removed
# ============================================================

fit_removed_iiv_model <- function(base_opts, obs_data, removed_iiv) {
  
  all_params <- c("cl", "v1", "v2", "q", "ka")
  random_params <- setdiff(all_params, removed_iiv)
  
  opts <- do.call(
    genopts,
    modifyList(
      base_opts,
      list(
        g = function(beta, bi = rep(0, length(random_params)), ai) {
          
          theta <- beta[all_params]
          names(theta) <- all_params
          
          for (j in seq_along(random_params)) {
            theta[random_params[j]] <- beta[random_params[j]] * exp(bi[j])
          }
          
          return(theta)
        },
        
        single_betas = all_params == removed_iiv,
        
        p_thetai = function(p, origbeta, bi) {
          
          dmnorm(
            bi,
            mean = log(p$beta[random_params] / origbeta[random_params]),
            sigma = p$Omega,
            log = TRUE
          )$den
        }
      )
    )
  )
  
  fit <- admr::fitMC(
    opts = opts,
    obs = obs_data,
    chains = 1,
    maxiter = 200,
    use_grad = TRUE
  )
  
  return(fit)
}
