#script for running simulation study for variable selection Project

#read in Simulation_Functions script
source("./Project 4/Code/Simulation_Functions.R")

#packages needed for sim
library(furrr)
library(future)
##simulation parameters
#set seed for reproducibility
set.seed(6624)

#data generating parameters
n_vals <- c(250, 500)
rho_vals <- c(0, 0.35, 0.7)
true_betas <- c(seq(0.5/3, 2.5/3, 0.5/3), rep(0, 15))
p <- 20
p1 <- 5
nreps <- 10000


true_vars <- paste0("X", 1:5)
noise_vars <- paste0("X", 6:20)

# Create simulation grid
sim_grid <- expand.grid(
  n = n_vals,
  rho = rho_vals,
  rep_id = 1:nreps
)

# For testing first:
# sim_grid <- sim_grid[1:20, ]

# Parallel setup
plan(multisession, workers = parallel::detectCores() - 1)

# Run simulation in parallel
final_results_df <- future_pmap_dfr(
  sim_grid,
  function(n, rho, rep_id) {
    run_one_rep(
      rep_id = rep_id,
      n = n,
      rho = rho,
      true_betas = true_betas,
      p = p,
      p1 = p1
    )
  },
  .options = furrr_options(seed = TRUE)
)

# Add helper columns
final_results_df <- final_results_df %>%
  mutate(
    truly_active = variable %in% true_vars,
    noise_variable = variable %in% noise_vars,
    ci_covers = if_else(
      selected,
      ci_low <= true_beta & ci_hi >= true_beta,
      NA
    ),
    reject_null = if_else(
      selected,
      p_value < 0.05,
      FALSE
    )
  )

# Save compact output
saveRDS(
  final_results_df,
  "./Project 4/Data/p4_simulation_results.rds"
)


# Selection frequency by variable
selection_summary <- final_results_df %>%
  group_by(n, rho, method, variable, true_beta, truly_active) %>%
  summarise(
    selection_rate = mean(selected),
    .groups = "drop"
  )

# True positive rates for X1-X5
true_positive_summary <- final_results_df %>%
  filter(truly_active) %>%
  group_by(n, rho, method, variable, true_beta) %>%
  summarise(
    true_positive_rate = mean(selected),
    .groups = "drop"
  )

# False positive selection rates for X6-X20
false_positive_summary <- final_results_df %>%
  filter(noise_variable) %>%
  group_by(n, rho, method, variable) %>%
  summarise(
    false_positive_selection_rate = mean(selected),
    .groups = "drop"
  )

# Bias and coverage among selected variables only
bias_coverage_summary <- final_results_df %>%
  filter(selected) %>%
  group_by(n, rho, method, variable, true_beta, truly_active) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - true_beta, na.rm = TRUE),
    coverage = mean(ci_covers, na.rm = TRUE),
    .groups = "drop"
  )

# Type I error for noise variables
type1_summary <- final_results_df %>%
  filter(noise_variable) %>%
  group_by(n, rho, method, variable) %>%
  summarise(
    type1_error = mean(selected & p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# Power and type II error for true variables
power_type2_summary <- final_results_df %>%
  filter(truly_active) %>%
  group_by(n, rho, method, variable, true_beta) %>%
  summarise(
    power = mean(selected & p_value < 0.05, na.rm = TRUE),
    type2_error = 1 - power,
    .groups = "drop"
  )

# Overall method-level summaries
method_summary <- final_results_df %>%
  group_by(n, rho, method) %>%
  summarise(
    avg_true_positive_rate = mean(selected[truly_active]),
    avg_false_positive_rate = mean(selected[noise_variable]),
    avg_type1_error = mean((selected & p_value < 0.05)[noise_variable], na.rm = TRUE),
    avg_power = mean((selected & p_value < 0.05)[truly_active], na.rm = TRUE),
    avg_type2_error = 1 - avg_power,
    .groups = "drop"
  )
