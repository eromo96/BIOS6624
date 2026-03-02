####CODE FOR DATA ANALYSIS 
library(cmdstanr)
library(parallel)
library(brms)
library(posterior)
library(bayestestR)
library(loo)
library(broom)
library(kableExtra)
#read in analysis dataset created in Exploration.R script
hiv_andata <- read_csv("./Project1/Data/hiv_andata.csv")
hiv_andata <- hiv_andata %>%
  mutate(
    hard_drugs_cat_0 = factor(hard_drugs_cat_0),
    smoke_cat_0      = factor(smoke_cat_0),
    edu_cat_0        = factor(edu_cat_0),
    race_cat_0       = factor(race_cat_0),
    adh_cat_2        = factor(adh_cat_2)
  )
#set default priors for models but will need to specify different prior for d_LEU3N and r_logVLOAD
default_priors <- c(
  prior(normal(0, 100), class = "Intercept"),
  prior(normal(0, 100), class = "b"),
  prior(normal(0, 100), class = "sigma", lb = 0)
)

#create function that sets the same non-informative priors across models, saves model to "Models" folder, and runs model
run_brm <- function(formula, data, file_name, family = gaussian(), priors = default_priors,
                    iter = 20000, warmup = 5000, chains = 4, cores = min(4, parallel::detectCores()),
                    seed = 6624, adapt_delta = 0.95, max_treedepth = 12, models_dir = "./Project1/Models",
                    force_refit = FALSE){
  file_path <- file.path(models_dir, paste0("brms_", file_name, ".rds"))
  if (file.exists(file_path) && !force_refit) {
    message("Reading cached model: ", file_path)
    return(readRDS(file_path))
  }
  message("Fitting model: ", deparse(formula), "\nSaving to: ", file_path)
  
  fit <- brm(formula = formula, data = data, family = family, prior = priors, backend = "cmdstanr",
             iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
             control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
             save_pars = save_pars(all = TRUE))
  
  saveRDS(fit, file_path)
  fit
}

#non adherence model formulas
formulas_no_adh <- list(
  d_AGG_MENT = d_AGG_MENT ~ hard_drugs_cat_0 + AGG_MENT_0 + age_0 + BMI_0 + smoke_cat_0 + edu_cat_0 + race_cat_0,
  d_AGG_PHYS = d_AGG_PHYS ~ hard_drugs_cat_0 + AGG_PHYS_0 + age_0 + BMI_0 + smoke_cat_0 + edu_cat_0 + race_cat_0,
  d_LEU3N = d_LEU3N ~ hard_drugs_cat_0 + LEU3N_0 + age_0 + BMI_0 + smoke_cat_0 + edu_cat_0 + race_cat_0,
  r_logVLOAD = r_logVLOAD ~ hard_drugs_cat_0 + logVLOAD_0 + age_0 + BMI_0 + smoke_cat_0 + edu_cat_0 + race_cat_0
)
#adherence model formulas
formulas_with_adh <- list(
  d_AGG_MENT = update(formulas_no_adh$d_AGG_MENT, . ~ . + adh_cat_2),
  d_AGG_PHYS = update(formulas_no_adh$d_AGG_PHYS, . ~ . + adh_cat_2),
  d_LEU3N = update(formulas_no_adh$d_LEU3N, . ~ . + adh_cat_2),
  r_logVLOAD = update(formulas_no_adh$r_logVLOAD, . ~ . + adh_cat_2)
)

#initialize list that will store brms fits
brms_fits <- list()

#loop through each formula and run our personalized brms function
for(mods in names(formulas_no_adh)) {
  brms_fits[[paste0(mods, "_no_adh")]] <- 
    run_brm(formula = formulas_no_adh[[mods]], data = hiv_andata, file_name = paste0(mods, "_no_adh"),
            force_refit = TRUE)
  brms_fits[[paste0(mods, "_with_adh")]] <- 
    run_brm(formula = formulas_with_adh[[mods]], data = hiv_andata, file_name = paste0(mods, "_with_adh"),
            force_refit = TRUE)
}


#function that runs linear regression and saves model output to models folder
run_lm <- function(formula, data, file_name, models_dir = "./Project1/Models",
                   force_refit = FALSE) {
  file_path <- file.path(models_dir, paste0("lm_", file_name, ".rds"))
  if (file.exists(file_path) && !force_refit) {
    message("Reading cached lm: ", file_path)
    return(readRDS(file_path))
  }
  message("Fitting lm: ", deparse(formula), "\nSaving to: ", file_path)
  fit <- lm(formula = formula, data = data)
  saveRDS(fit, file_path)
  fit
}

#initialize list for lm fits
lm_fits <- list()

#loop through each formula and run our personalized lm function
for(mods in names(formulas_no_adh)) {
  lm_fits[[paste0(mods, "_no_adh")]] <- 
    run_lm(formula = formulas_no_adh[[mods]], data = hiv_andata, file_name = paste0(mods, "_no_adh"),
           force_refit = TRUE)
  lm_fits[[paste0(mods, "_with_adh")]] <- 
    run_lm(formula = formulas_with_adh[[mods]], data = hiv_andata, file_name = paste0(mods, "_with_adh"),
           force_refit = TRUE)
}


# clinically meaningful thresholds
clin_thresh <- c(
  d_AGG_MENT = 2,
  d_AGG_PHYS = 2,
  d_LEU3N = 50,
  r_logVLOAD = 0.5
)

#create nice outcome labels for tables
outcome_labs <- function(x) {
  recode(x,
         d_AGG_MENT = "Δ AGG_MENT",
         d_AGG_PHYS = "Δ AGG_PHYS",
         d_LEU3N    = "Δ CD4 (LEU3N)",
         r_logVLOAD = "log10(VLOAD2/VLOAD0)")
}

#loo helper function
looic_cached_brm <- function(fit,
                             file_name,
                             models_dir = "./Project1/Models",
                             reloo = FALSE,
                             force_refit = FALSE) {
  
  dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
  file_path <- file.path(models_dir, paste0("brms_", file_name, ".rds"))
  
  # If we already have a saved brmsfit with loo criterion added, reuse it
  if (file.exists(file_path) && !force_refit) {
    fit2 <- readRDS(file_path)
    if (!is.null(fit2$criteria$loo)) {
      loo <- fit2$criteria$loo
      looic <- loo$estimates["looic", "Estimate"]
      looic_se <- loo$estimates["looic", "SE"]
      return(list(looic = looic, looic_se = looic_se, fit = fit2))
    }
    # else: fall through and add criterion below
    fit <- fit2
  }
  
  # Add criterion using brms (uses log_lik; does NOT require rstan)
  fit2 <- brms::add_criterion(fit, "loo", reloo = reloo)
  
  # Save updated fit so you never recompute LOO again
  saveRDS(fit2, file_path)
  
  loo <- fit2$criteria$loo
  looic <- loo$estimates["looic", "Estimate"]
  looic_se <- loo$estimates["looic", "SE"]
  
  list(looic = looic, looic_se = looic_se, fit = fit2)
}

delta_looic_vals <- function(looic_full, looic_red) {
  looic_red - looic_full
}

#extract frequentist model output
lm_term_stats <- function(fit, term) {
  
  sm <- summary(fit)
  coefs <- sm$coefficients
  
  if (!term %in% rownames(coefs)) {
    stop("Term not found in summary(fit): ", term,
         "\nAvailable: ", paste(rownames(coefs), collapse = ", "))
  }
  
  est <- coefs[term, "Estimate"]
  se  <- coefs[term, "Std. Error"]
  p   <- coefs[term, "Pr(>|t|)"]
  
  ci <- confint(fit, parm = term)
  
  tibble::tibble(
    est_f = est,
    se_f  = se,
    l95_f = ci[1],
    u95_f = ci[2],
    p_f   = p
  )
}

#find the unique coefficient term in lm() by pattern
find_lm_term <- function(fit, pattern) {
  trms <- names(coef(fit))
  hit <- grep(pattern, trms, value = TRUE)
  if (length(hit) != 1) stop ("Pattern matched ", length(hit), " terms: ", paste(hit, collapse = ", "))
  hit
}

#extract bayesian coefficient stats
brms_coef_stats <- function(fit, coef_regex, clin_thresh,
                            include_mcse = FALSE) {
  coef_names <- rownames(brms::fixef(fit))
  hit <- grep(coef_regex, coef_names, value = TRUE)
  
  if (length(hit) != 1) {
    stop("Coefficient regex matched ", length(hit), " terms. Regex=", coef_regex,
         "\nMatches: ", paste(hit, collapse = ", "))
  }
  
  coef_name <- hit
  draws <- posterior::as_draws_df(fit)
  colname <- paste0("b_", coef_name)
  if (!colname %in% names(draws)) stop("Posterior draw column not found: ", colname)
  
  b <- draws[[colname]]
  
  est <- mean(b)
  sd  <- sd(b)
  
  h <- bayestestR::hdi(b, ci = 0.95, method = "HPD")
  
  # Keep ONE clinically meaningful posterior probability:
  p_abs <- mean(abs(b) > clin_thresh)
  
  mcse <- if (include_mcse) posterior::mcse_mean(posterior::as_draws_vector(b)) else NA_real_
  
  tibble::tibble(
    est_b = est,
    sd_b  = sd,
    mcse_b = mcse,
    l95_b = h$CI_low,
    u95_b = h$CI_high,
    pr_abs_gt_clin = p_abs
  )
}

#function to create table for main analysis (effect of hard drug use on treatment effect)
build_main_table <- function(outcomes, formulas_no_adh, lm_fits, brms_fits,
                             run_brm, hiv_andata,
                             models_dir = "./Project1/Models",
                             include_mcse = FALSE,
                             reloo = FALSE) {
  
  # reduced models without hard drugs
  formulas_no_drugs <- lapply(formulas_no_adh, function(f) update(f, . ~ . - hard_drugs_cat_0))
  brms_red <- list()
  
  for (y in outcomes) {
    brms_red[[y]] <- run_brm(
      formula = formulas_no_drugs[[y]],
      data = hiv_andata,
      file_name = paste0(y, "_no_adh_no_drugs"),
      models_dir = models_dir
    )
  }
  
  tab <- purrr::map_dfr(outcomes, function(y) {
    
    # Frequentist (no-adh)
    lm_fit <- lm_fits[[paste0(y, "_no_adh")]]
    hd_term <- find_lm_term(lm_fit, "^hard_drugs_cat_0")
    lm_row  <- lm_term_stats(lm_fit, term = hd_term)
    
    # Bayesian (no-adh)
    brm_fit <- brms_fits[[paste0(y, "_no_adh")]]
    b_row <- brms_coef_stats(
      brm_fit,
      coef_regex = "^hard_drugs_cat_0",
      clin_thresh = clin_thresh[[y]],
      include_mcse = include_mcse
    )
    
    # LOOIC via brms::add_criterion (cached in saved brmsfit)
    loo_full <- looic_cached_brm(
      brm_fit,
      file_name = paste0(y, "_no_adh"),     
      models_dir = models_dir,
      reloo = reloo
    )
    
    loo_red <- looic_cached_brm(
      brms_red[[y]],
      file_name = paste0(y, "_no_adh_no_drugs"),
      models_dir = models_dir,
      reloo = reloo
    )
    
    tibble::tibble(
      outcome = y,
      delta_looic_red_minus_full = delta_looic_vals(loo_full$looic, loo_red$looic)
    ) %>%
      dplyr::bind_cols(lm_row, b_row)
  })
  
  tab_fmt <- tab %>%
    dplyr::mutate(
      Outcome = outcome_labs(outcome),
      
      `Freq: Est (SE)` = sprintf("%.2f (%.2f)", est_f, se_f),
      `Freq: 95% CI`   = sprintf("[%.2f, %.2f]", l95_f, u95_f),
      `Freq: p`        = format.pval(p_f, digits = 3, eps = 0.001),
      
      `Bayes: Mean (SD)` = sprintf("%.2f (%.2f)", est_b, sd_b),
      `Bayes: 95% HPDI`  = sprintf("[%.2f, %.2f]", l95_b, u95_b),
      `Pr(|β| > clin)`   = sprintf("%.3f", pr_abs_gt_clin),
      
      `ΔLOOIC (red-full)`= sprintf("%.1f", delta_looic_red_minus_full)
    ) %>%
    dplyr::select(
      Outcome,
      `Freq: Est (SE)`, `Freq: 95% CI`, `Freq: p`,
      `Bayes: Mean (SD)`, `Bayes: 95% HPDI`, `Pr(|β| > clin)`,
      `ΔLOOIC (red-full)`
    )
  
  tab_fmt
}

#function for secondary analysis table ("mediation")
build_secondary_table <- function(outcomes, lm_fits,
                                  models_dir = "./Project1/Models",
                                  brms_fits = brms_fits,
                                  include_mcse = FALSE,
                                  reloo = FALSE) {
  
  tab <- purrr::map_dfr(outcomes, function(y) {
    lm_with <- lm_fits[[paste0(y, "_with_adh")]]
    hd_term_with <- find_lm_term(lm_with, "^hard_drugs_cat_0")
    lm_row_with  <- lm_term_stats(lm_with, term = hd_term_with)
    
    # Also grab hard drug estimate from NO-ADH model to compute change in coef
    lm_no <- lm_fits[[paste0(y, "_no_adh")]]
    hd_term_no <- find_lm_term(lm_no, "^hard_drugs_cat_0")
    lm_row_no  <- lm_term_stats(lm_no, term = hd_term_no)
    
    delta_beta_f <- lm_row_with$est_f - lm_row_no$est_f
    brm_with <- brms_fits[[paste0(y, "_with_adh")]]
    b_row_with <- brms_coef_stats(
      brm_with,
      coef_regex  = "^hard_drugs_cat_0",
      clin_thresh = clin_thresh[[y]],
      include_mcse = include_mcse
    )
    
    # Also compute Bayesian change in posterior mean (with_adh - no_adh)
    brm_no <- brms_fits[[paste0(y, "_no_adh")]]
    b_row_no <- brms_coef_stats(
      brm_no,
      coef_regex  = "^hard_drugs_cat_0",
      clin_thresh = clin_thresh[[y]],
      include_mcse = include_mcse
    )
    
    delta_beta_b <- b_row_with$est_b - b_row_no$est_b
    loo_full <- looic_cached_brm(
      brm_with,
      file_name = paste0(y, "_with_adh"),
      models_dir = models_dir,
      reloo = reloo
    )
    
    loo_red <- looic_cached_brm(
      brm_no,
      file_name = paste0(y, "_no_adh"),
      models_dir = models_dir,
      reloo = reloo
    )
    
    tibble::tibble(
      outcome = y,
      delta_beta_f = delta_beta_f,
      delta_beta_b = delta_beta_b,
      delta_looic_red_minus_full = delta_looic_vals(loo_full$looic, loo_red$looic)
    ) %>%
      dplyr::bind_cols(lm_row_with, b_row_with)
  })
  
  tab_fmt <- tab %>%
    dplyr::mutate(
      Outcome = outcome_labs(outcome),
      
      `Freq: Est (SE)` = sprintf("%.2f (%.2f)", est_f, se_f),
      `Freq: 95% CI`   = sprintf("[%.2f, %.2f]", l95_f, u95_f),
      `Freq: p`        = format.pval(p_f, digits = 3, eps = 0.001),
      
      # change in hard drug coefficient after adding adherence
      `Δβ Freq (with-no)` = sprintf("%.2f", delta_beta_f),
      
      `Bayes: Mean (SD)` = sprintf("%.2f (%.2f)", est_b, sd_b),
      `Bayes: 95% HPDI`  = sprintf("[%.2f, %.2f]", l95_b, u95_b),
      `Pr(|β| > clin)`   = sprintf("%.3f", pr_abs_gt_clin),
      
      `Δβ Bayes (with-no)` = sprintf("%.2f", delta_beta_b),
      
      `ΔLOOIC (no-with)` = sprintf("%.1f", delta_looic_red_minus_full)
    ) %>%
    dplyr::select(
      Outcome,
      `Freq: Est (SE)`, `Freq: 95% CI`, `Freq: p`, `Δβ Freq (with-no)`,
      `Bayes: Mean (SD)`, `Bayes: 95% HPDI`, `Pr(|β| > clin)`, `Δβ Bayes (with-no)`,
      `ΔLOOIC (no-with)`
    )
  
  tab_fmt
}

#create both tables
outcomes <- names(formulas_no_adh)

main_tbl_fmt <- build_main_table(
  outcomes = outcomes,
  formulas_no_adh = formulas_no_adh,
  lm_fits = lm_fits,
  brms_fits = brms_fits,
  run_brm = run_brm,
  hiv_andata = hiv_andata,
  models_dir = "./Project1/Models",
  include_mcse = FALSE,
  reloo = FALSE
)

secondary_tbl_fmt <- build_secondary_table(
  outcomes = outcomes,
  lm_fits = lm_fits,
  brms_fits = brms_fits,
  models_dir = "./Project1/Models",
  include_mcse = FALSE,
  reloo = FALSE
)

#use kableExtra to print nicely formatted tables and will work in rmarkdown to knit to pdf
kbl(main_tbl_fmt, booktabs = TRUE,
    caption = "Main analysis: hard drug use coefficient by outcome (Frequentist vs Bayesian).") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

kbl(secondary_tbl_fmt, booktabs = TRUE,
    caption = "Secondary analysis (naive mediation): adherence coefficient by outcome (Frequentist vs Bayesian).") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))



#### model diagnostics

#make diagnostic checking more efficient using loop
for (nm in names(brms_fits)) {
  
  fit <- brms_fits[[nm]]
  
  diag_tbl <- posterior::summarise_draws(
    posterior::as_draws_array(fit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  
  cat("\n", nm, "\n")
  cat("Max Rhat:", max(diag_tbl$rhat, na.rm = TRUE), "\n")
  cat("Min Bulk ESS:", min(diag_tbl$ess_bulk, na.rm = TRUE), "\n")
  cat("Min Tail ESS:", min(diag_tbl$ess_tail, na.rm = TRUE), "\n")
  
  # ACF plot
  p1 <- brms::mcmc_plot(
    fit,
    type = "acf",
    variable = c("b_hard_drugs_cat_0Yes", "sigma")
  )
  print(p1)
  
  # Trace plot
  p2 <- brms::mcmc_plot(
    fit,
    type = "trace",
    variable = c("b_hard_drugs_cat_0Yes", "sigma")
  )
  print(p2)
}

for (nm in names(brms_fits)) {
  
  cat("\n====================\n")
  cat("MODEL:", nm, "\n")
  cat("====================\n")
  
  fit <- brms_fits[[nm]]
  
  print(summary(fit))
  
}

