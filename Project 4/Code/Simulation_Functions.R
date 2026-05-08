#create all the functions that will be used to run the simulation for project 4

#load necessary packages
library(glmnet)
library(MASS)
library(purrr)
library(tibble)
library(broom)
library(tidyverse)

#test out gen_data function from hdmr
# testdata<- hdrm::gen_data(100, 10, 5)

# Global variable names
x_vars <- paste0("X", 1:20)

#data generating function
gen_dataset <- function(n, rho, p = 20, p1 = 5, true_betas) {
  #use gen_data function from hdmr package
  dat <- hdrm::gen_data(
    n = n,
    p = p,
    p1 = p1,
    beta = true_betas,
    family = "gaussian",
    corr = "exchangeable",
    rho = rho
  )
  
  X <- as.data.frame(dat$X)
  names(X) <- paste0("X", 1:p)
  
  sim_dat <- data.frame(
    y = dat$y,
    X
  )
  return(sim_dat)
}

#function that refits selected variables from model selection method with lm()
refit_extract <- function(sim_dat, selected_vars, method, n, rho, rep_id, true_betas) {
  
  all_vars <- paste0("X", seq_along(true_betas))
  
  # Base output: one row for every variable
  results <- tibble(
    rep_id = rep_id,
    n = n,
    rho = rho,
    method = method,
    variable = all_vars,
    selected = variable %in% selected_vars,
    estimate = NA_real_,
    p_value = NA_real_,
    ci_low = NA_real_,
    ci_hi = NA_real_,
    true_beta = true_betas
  )
  
  # If no variables selected, return all rows with selected = FALSE
  if (length(selected_vars) == 0) {
    return(results)
  }
  
  # Refit final no-intercept model
  form <- as.formula(
    paste("y ~ 0 +", paste(selected_vars, collapse = " + "))
  )
  
  final_fit <- lm(form, data = sim_dat)
  
  coefs_df <- broom::tidy(final_fit, conf.int = TRUE) %>%
    rename(
      variable = term,
      ci_low = conf.low,
      ci_hi = conf.high,
      p_value = p.value
    ) %>%
    select(variable, estimate, p_value, ci_low, ci_hi)
  
  # Join estimates back onto all 20 variables
  results <- results %>%
    select(-estimate, -p_value, -ci_low, -ci_hi) %>%
    left_join(coefs_df, by = "variable")
  
  return(results)
}

#function for backward selection using F-test p-values
select_backward_p <- function(sim_dat, p_thresh = 0.1) {
  #start with full no-intercept model
  current_vars <- x_vars
  while (length(current_vars) > 0) {
    
    form <- as.formula(
      paste("y ~ 0 +", paste(current_vars, collapse = " + "))
    )
    
    fit <- lm(form, data = sim_dat)
    pvals <- summary(fit)$coefficients[, "Pr(>|t|)"]
    
    max_p <- max(pvals, na.rm = TRUE)
    
    if (max_p <= p_thresh) {
      break
    }
    
    var_remove <- names(which.max(pvals))
    current_vars <- setdiff(current_vars, var_remove)
  }
  return(current_vars)
}

#function for aic backward selection
select_backward_aic <- function(sim_dat) {
  
  full_fit <- lm(
    as.formula(paste("y ~ 0 +", paste(x_vars, collapse = " + "))),
    data = sim_dat
  )
  
  null_fit <- lm(y ~ 0, data = sim_dat)
  
  step_fit <- step(
    full_fit,
    scope = list(
      lower = formula(null_fit),
      upper = formula(full_fit)
    ),
    direction = "backward",
    k = 2,
    trace = 0
  )
  
  selected_vars <- names(coef(step_fit))
  return(selected_vars)
}

#function for bic backward selection
select_backward_bic <- function(sim_dat) {
  
  n <- nrow(sim_dat)
  
  full_fit <- lm(
    as.formula(paste("y ~ 0 +", paste(x_vars, collapse = " + "))),
    data = sim_dat
  )
  
  null_fit <- lm(y ~ 0, data = sim_dat)
  
  step_fit <- step(
    full_fit,
    scope = list(
      lower = formula(null_fit),
      upper = formula(full_fit)
    ),
    direction = "backward",
    k = log(n),
    trace = 0
  )
  
  selected_vars <- names(coef(step_fit))
  return(selected_vars)
}

# Fit lasso once and extract lambda.min and lambda.1se
select_lasso_both <- function(sim_dat) {
  
  X <- as.matrix(sim_dat[, x_vars])
  y <- sim_dat$y
  
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = 1,
    intercept = FALSE,
    standardize = TRUE
  )
  
  beta_min <- as.matrix(coef(cv_fit, s = "lambda.min"))[-1, , drop = FALSE]
  beta_1se <- as.matrix(coef(cv_fit, s = "lambda.1se"))[-1, , drop = FALSE]
  
  list(
    lasso_min = rownames(beta_min)[beta_min[, 1] != 0],
    lasso_1se = rownames(beta_1se)[beta_1se[, 1] != 0]
  )
}

# Fit elastic net once and extract lambda.min and lambda.1se
select_enet_both <- function(sim_dat) {
  
  X <- as.matrix(sim_dat[, x_vars])
  y <- sim_dat$y
  
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = 0.5,
    intercept = FALSE,
    standardize = TRUE
  )
  
  beta_min <- as.matrix(coef(cv_fit, s = "lambda.min"))[-1, , drop = FALSE]
  beta_1se <- as.matrix(coef(cv_fit, s = "lambda.1se"))[-1, , drop = FALSE]
  
  list(
    enet_min = rownames(beta_min)[beta_min[, 1] != 0],
    enet_1se = rownames(beta_1se)[beta_1se[, 1] != 0]
  )
}

#function to run all model selection methods on one simulated dataset
run_one_rep <- function(rep_id, n, rho, true_betas, p = 20, p1 = 5) {
  
  sim_dat <- gen_dataset(
    n = n,
    rho = rho,
    p = p,
    p1 = p1,
    true_betas = true_betas
  )
  
  all_results <- list()
  
  # Backward p-value
  selected <- select_backward_p(sim_dat, p_thresh = 0.15)
  all_results[["backward_pvalue"]] <- refit_extract(
    sim_dat, selected, "backward_pvalue", n, rho, rep_id, true_betas
  )
  
  # Backward AIC
  selected <- select_backward_aic(sim_dat)
  all_results[["backward_aic"]] <- refit_extract(
    sim_dat, selected, "backward_aic", n, rho, rep_id, true_betas
  )
  
  # Backward BIC
  selected <- select_backward_bic(sim_dat)
  all_results[["backward_bic"]] <- refit_extract(
    sim_dat, selected, "backward_bic", n, rho, rep_id, true_betas
  )
  
  # Lasso
  lasso_selected <- select_lasso_both(sim_dat)
  
  all_results[["lasso_min"]] <- refit_extract(
    sim_dat, lasso_selected$lasso_min, "lasso_min", n, rho, rep_id, true_betas
  )
  
  all_results[["lasso_1se"]] <- refit_extract(
    sim_dat, lasso_selected$lasso_1se, "lasso_1se", n, rho, rep_id, true_betas
  )
  
  # Elastic net
  enet_selected <- select_enet_both(sim_dat)
  
  all_results[["enet_min"]] <- refit_extract(
    sim_dat, enet_selected$enet_min, "enet_min", n, rho, rep_id, true_betas
  )
  
  all_results[["enet_1se"]] <- refit_extract(
    sim_dat, enet_selected$enet_1se, "enet_1se", n, rho, rep_id, true_betas
  )
  
  bind_rows(all_results)
}