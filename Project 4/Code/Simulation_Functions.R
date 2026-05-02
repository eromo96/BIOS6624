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

#data generating function
gen_dataset <- function(n, rho) {
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
  names(X) <- paste0("X", 1:20)
  
  sim_dat <- data.frame(
    y = dat$y,
    X
  )
  return(sim_dat)
}

#function that refits selected variables from model selection method with lm()
refit_extract <- function(sim_dat, selected_vars, method_name, n, rho, rep_id) {
  if(length(selected_vars) == 0) {
    print("No variables selected")
  }
  
  #build formula with no intercept
  form <- as.formula(paste("y ~ 0 +", paste(selected_vars, collapse = " + ")))
  
  #refit final model
  final_fit <- lm(formula = form, data = sim_dat)
  
  #extract coefficient table
  coefs_df <- broom::tidy(final_fit, conf.int = TRUE)
  
  #create compact results df
  results <- data.frame(
    rep_id = rep_id,
    n = n,
    rho = rho,
    method_name = method_name,
    variable = coefs_df$term,
    estimate = coefs_df$estimate,
    p_value = coefs_df$p.value,
    ci_low = coefs_df$conf.low,
    ci_hi = coefs_df$conf.high,
    selected = TRUE
  )
  return(results)
}

#function for backward selection using F-test p-values
select_backward_p <- function(sim_dat, p_thresh = 0.15) {
  #start with full no-intercept model
  current_vars <- paste0("X", 1:20)
  continue <- TRUE
  while (continue) {
    #fit current model
    form <- as.formula(paste("y ~ 0 +", paste(current_vars, collapse = " + ")))
    fit <- lm(form, data = sim_dat)
    pvals <- summary(fit)$coefficients[, "Pr(>|t|)"]
    max_p <- max(pvals)
    #if max p is less than threshold then all vars are sig so stop selection
    if(max_p <= p_thresh) {
      continue <- FALSE
    }
    else {
      var_remove <- names(which.max(pvals))
      current_vars <- setdiff(current_vars, var_remove)
    }
    if(length(current_vars) == 0) {
      continue <- FALSE
    }
  }
  return(current_vars)
}

#function for aic backward selection
select_backward_aic <- function(sim_dat) {
  all_vars <- paste0("X", 1:20)
  full_fit <- lm(formula = as.formula(paste("y ~ 0 +", paste(all_vars, collapse = " + "))),
                 data = sim_dat)
  null_fit <- lm(y ~ 0, data = sim_dat)
  #backward selection via aic
  step_fit <- step(
    full_fit,
    scope = list(lower = formula(null_fit), upper = formula(full_fit)),
    direction = "backward",
    trace = 0
  )
  #extract coefficients
  selected_vars <- names(coef(step_fit))
  return(selected_vars)
}

#function for bic backward selection
select_backward_bic <- function(sim_dat) {
  n <- nrow(sim_dat)
  all_vars <- paste0("X", 1:20)
  full_fit <- lm(formula = as.formula(paste("y ~ 0 +", paste(all_vars, collapse = " + "))),
                 data = sim_dat)
  null_fit <- lm(y ~ 0, data = sim_dat)
  #backward selection via aic
  step_fit <- step(
    full_fit,
    scope = list(lower = formula(null_fit), upper = formula(full_fit)),
    direction = "backward",
    trace = 0,
    k = log(n)
  )
  #extract coefficients
  selected_vars <- names(coef(step_fit))
  return(selected_vars)
}

#function for lasso using cv.glmnet with alpha = 1
select_lasso <- function(sim_dat, lambda_choice) {
  X <- as.matrix(sim_dat[, paste0("X", 1:20)])
  y <- sim_dat$y
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = 1,
    intercept = FALSE,
    standardize = TRUE
  )
  if(lambda_choice == "lambda.min") {
    lambda_use <- cv_fit$lambda.min
  }
  if(lambda_choice == "lambda.1se") {
    lambda_use <- cv_fit$lambda.1se
  }
  beta_hat <- coef(cv_fit, s = lambda_use)
  
  #remove intercept row
  beta_hat <- beta_hat[-1, ]
  selected_vars <- names(beta_hat[beta_hat != 0])
  return(selected_vars)
}

#function for elastic net
select_elastic_net <- function(sim_dat, lambda_choice) {
  X <- as.matrix(sim_dat[, paste0("X", 1:20)])
  y <- sim_dat$y
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = 0.5,
    intercept = FALSE,
    standardize = TRUE
  )
  if(lambda_choice == "lambda.min") {
    lambda_use <- cv_fit$lambda.min
  }
  if(lambda_choice == "lambda.1se") {
    lambda_use <- cv_fit$lambda.1se
  }
  beta_hat <- coef(cv_fit, s = lambda_use)
  
  #remove intercept row
  beta_hat <- beta_hat[-1, ]
  selected_vars <- names(beta_hat[beta_hat != 0])
  return(selected_vars)
}

#function to run all model selection methods on one simulated dataset
run_one_rep <- function(rep_id, n, rho) {
  
}