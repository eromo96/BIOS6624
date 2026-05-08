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

## For testing first:
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

#add labels to method column
final_results_df <- final_results_df %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "backward_pvalue", "backward_aic", "backward_bic",
        "lasso_min", "lasso_1se", "enet_min", "enet_1se"
      ),
      labels = c("Backward P-Value", "Backward AIC", "Backward BIC",
                 "LASSO Min", "LASSO 1SE", "Elastic Net Min", "Elastic Net 1SE")
    )
  )


#replication level selection summary
rep_level_selection <- final_results_df %>%
  group_by(n, rho, method, rep_id) %>%
  summarise(
    #selection metrics
    TPR = mean(selected[truly_active]),
    FPR = mean(selected[noise_variable]),
    .groups = "drop"
  )

#aggregate summary to method level across reps
method_level_selection <- rep_level_selection %>%
  group_by(n, rho, method) %>%
  summarise(
    TPR = mean(TPR),
    FPR = mean(FPR),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = paste0("n=", n, ", rho=", rho),
    scenario = factor(
      scenario,
      levels = c(
        "n=250, rho=0", "n=250, rho=0.35", "n=250, rho=0.7",
        "n=500, rho=0", "n=500, rho=0.35", "n=500, rho=0.7"
      )
    )
  )

#heat map for proportion of correct selections and incorrect selections
heatmap_selection <- method_level_selection %>%
  select(method, scenario, TPR, FPR) %>%
  pivot_longer(
    cols = c(TPR, FPR),
    names_to = "metric",
    values_to = "value"
  )



#heatmap for tpr/fpr
ggplot(heatmap_selection, aes(x = scenario, y = method, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  facet_wrap(~ metric) +
  labs(
    title = "True Positive and False Positive Selection Rates",
    x = "Scenario",
    y = "Selection Method",
    fill = "Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#variable level type 1/type 2 error
variable_error_summary <- final_results_df %>%
  mutate(
    covariate_group = case_when(
      variable == "X1" ~ "X1",
      variable == "X2" ~ "X2",
      variable == "X3" ~ "X3",
      variable == "X4" ~ "X4",
      variable == "X5" ~ "X5",
      variable %in% noise_vars ~ "X6-X20"
    ),
    error = case_when(
      truly_active ~ !reject_null,
      noise_variable ~ reject_null
    ),
    error_type = case_when(
      truly_active ~ "Type 2 Error",
      noise_variable ~ "Type 1 Error"
    )
  ) %>%
  group_by(n, rho, method, covariate_group, error_type) %>%
  summarise(
    error_rate = round(mean(error, na.rm = TRUE), 4),
    .groups = "drop"
  )

#pivot to wide format
variable_error_tab_wide <- variable_error_summary %>%
  select(n, rho, method, covariate_group, error_rate) %>%
  pivot_wider(
    names_from = c(n, covariate_group),
    values_from = error_rate
  ) %>%
  arrange(rho, method)
#prep the table to be ready for knitting
variable_error_tab_display <- variable_error_tab_wide %>%
  rename_with(~ gsub("250_", "n250_", .x)) %>%
  rename_with(~ gsub("500_", "n500_", .x)) %>%
  mutate(
    rho = paste0("$\\rho = ", rho, "$")
  )
variable_error_tab_display <- variable_error_tab_display %>%
  select(
    rho, method,
    n250_X1, n250_X2, n250_X3, n250_X4, n250_X5, `n250_X6-X20`,
    n500_X1, n500_X2, n500_X3, n500_X4, n500_X5, `n500_X6-X20`
  )

#bias and coverage summary 
bias_coverage_summary <- final_results_df %>%
  filter(truly_active, selected) %>%
  group_by(rho, n, method, variable, true_beta) %>%
  summarise(
    bias = round(mean(estimate - true_beta, na.rm = TRUE), 4),
    coverage = round(mean(ci_covers, na.rm = TRUE), 4),
    .groups = "drop"
  )
#create table 
bias_coverage_mega_table <- bias_coverage_summary %>%
  select(rho, n, method, variable, bias, coverage) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(bias, coverage)
  ) %>%
  arrange(rho, n, method)
#prep table to be ready for knitting
bias_coverage_tab_display <- bias_coverage_mega_table %>%
  mutate(
    rho = paste0("$\\rho = ", rho, "$"),
    n = paste0("$n = ", n, "$")
  ) %>%
  select(
    rho, n, method,
    bias_X1, coverage_X1,
    bias_X2, coverage_X2,
    bias_X3, coverage_X3,
    bias_X4, coverage_X4,
    bias_X5, coverage_X5
  )
