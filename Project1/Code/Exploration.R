####CODE FOR DATA EXPLORATION INCLUDING CREATING ANALYSIS DATASET AND EXPLORATION
library(knitr)
library(kableExtra)
library(stringr)
#read in code form cleaning script so that we dont have to redo data cleaning here
source("./Project1/Code/Cleaning.R")

#create analysis dataset
hivp1_data <- hiv_data %>%
  #keep only baseline and year 2 measurements
  filter(years %in% c(0,2)) %>%
  #keep only the columns we need
  #covariates will become *_0 and *_2 after pivoting
  pivot_wider(
    id_cols = newid,
    names_from = years,
    values_from = c(AGG_MENT, AGG_PHYS, LEU3N, VLOAD, income, BMI, SMOKE, smoke_cat,
                    ADH, adh_cat, RACE, race_cat, EDUCBAS, edu_cat, age, hard_drugs,
                    hard_drugs_cat),
    names_sep = "_"
  ) %>%
  #compute year2 - baseline diffs for outcomes
  mutate(
    d_AGG_MENT = AGG_MENT_2 - AGG_MENT_0,
    d_AGG_PHYS = AGG_PHYS_2 - AGG_PHYS_0,
    d_LEU3N = LEU3N_2 - LEU3N_0,
    d_VLOAD = VLOAD_2 - VLOAD_0
  ) %>%
  select(newid, ends_with("_0"), ends_with("_2"), starts_with("d_"))


#run histograms of outcomes to see if normally distributed
hist(hivp1_data$d_AGG_MENT) #looks approx normal
hist(hivp1_data$d_AGG_PHYS) #looks approx normal
hist(hivp1_data$d_LEU3N) #looks approx normal
hist(hivp1_data$d_VLOAD) #very skewed
# try log10 transformation of ratio of VLOAD_2 and VLOAD_0
hivp1_data <- hivp1_data %>%
  mutate(
    r_logVLOAD = log10(VLOAD_2/VLOAD_0),
    logVLOAD_0 = log10(VLOAD_0)
  )
hist(hivp1_data$r_logVLOAD) #looks a lot more normal now

# Boxplots of changes by baseline hard drug use
par(mfrow=c(2,2))
boxplot(d_AGG_MENT ~ hard_drugs_cat_0, data=hivp1_data, main="Δ AGG_MENT (Year2 - Baseline)", xlab="", ylab="Change")
boxplot(d_AGG_PHYS ~ hard_drugs_cat_0, data=hivp1_data, main="Δ AGG_PHYS (Year2 - Baseline)", xlab="", ylab="Change")
boxplot(d_LEU3N    ~ hard_drugs_cat_0, data=hivp1_data, main="Δ LEU3N (Year2 - Baseline)", xlab="", ylab="Change")
boxplot(r_logVLOAD ~ hard_drugs_cat_0, data=hivp1_data, main="log(VLOAD2/VLOAD0)", xlab="", ylab="Log ratio")
par(mfrow=c(1,1))

#missingness table code
vars_check <- c("income_0","BMI_0","age_0","smoke_cat_0","adh_cat_2","race_cat_0","edu_cat_0",
                "d_AGG_MENT","d_AGG_PHYS","d_LEU3N","r_logVLOAD")

miss_tbl <- hivp1_data %>%
  summarise(across(all_of(vars_check), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing_n") %>%
  mutate(
    N = nrow(hivp1_data),
    Missing_pct = 100 * Missing_n / N
  ) %>%
  arrange(desc(Missing_pct))

kbl(miss_tbl, format="latex", booktabs=TRUE, digits=1,
    caption="Missingness counts and percents for key analysis variables.",
    escape=FALSE) %>%
  kable_styling(latex_options=c("hold_position"))

#exclude observations that have missing values for outcomes for the descriptive
#summaries and analysis
hivp1_data <- hivp1_data %>%
  drop_na(d_AGG_MENT, d_AGG_PHYS, d_LEU3N, r_logVLOAD)

### ---- code for creating table 1 --------
# --- helper functions ---
fmt_n_pct <- function(n, denom) {
  if (is.na(n) | is.na(denom) | denom == 0) return(NA_character_)
  sprintf("%d (%.1f\\%%)", n, 100 * n / denom)
}

fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_missing <- function(x) {
  n_miss <- sum(is.na(x))
  p_miss <- mean(is.na(x)) * 100
  sprintf("%d (%.1f\\%%)", n_miss, p_miss)
}

# --- main Table 1 function ---
make_table1 <- function(df,
                        group_var = hard_drugs_cat,
                        cont_vars = c("BMI_0","age_0", "LEU3N_0", "logVLOAD_0", "AGG_MENT_0", "AGG_PHYS_0"),
                        cat_vars  = c("smoke_cat_0","adh_cat_2","race_cat_0","edu_cat_0"),
                        var_labels = list(
                          BMI_0    = "BMI",
                          age_0    = "Age",
                          smoke_cat_0 = "Smoking status",
                          adh_cat_2   = "Adherence category (year 2)",
                          race_cat_0  = "Race/ethnicity",
                          edu_cat_0   = "Education",
                          LEU3N_0 = "CD4+ Count",
                          logVLOAD_0 = "Log Viral Load",
                          AGG_MENT_0 = "Mental Quality of Life",
                          AGG_PHYS_0 = "Physical Quality of Life"
                        )) {
  
  # ensure group is a factor with nice labels
  df <- df %>%
    mutate(.group = as.factor({{ group_var }}))
  
  # group counts (non-missing group)
  n_overall <- nrow(df)
  n_bygrp <- df %>%
    filter(!is.na(.group)) %>%
    count(.group, name = "N")
  
  # ---------- continuous summaries ----------
  cont_tbl <- lapply(cont_vars, function(v) {
    lab <- var_labels[[v]] %||% v
    
    overall_mean_sd   <- fmt_mean_sd(df[[v]])
    overall_miss      <- fmt_missing(df[[v]])
    
    bygrp <- df %>%
      group_by(.group) %>%
      summarise(
        mean_sd = fmt_mean_sd(.data[[v]]),
        miss    = fmt_missing(.data[[v]]),
        .groups = "drop"
      ) %>%
      right_join(n_bygrp, by = ".group") %>%
      arrange(.group)
    
    # build rows: Mean (SD), Missing
    out <- tibble(
      Variable = c(lab, ""),
      Level = c("Mean (SD)", "Missing"),
      Overall = c(overall_mean_sd, overall_miss)
    )
    
    # add each group as columns (No/Yes)
    for (g in levels(df$.group)) {
      val_mean_sd <- bygrp$mean_sd[bygrp$.group == g]
      val_miss    <- bygrp$miss[bygrp$.group == g]
      out[[as.character(g)]] <- c(val_mean_sd, val_miss)
    }
    
    out
  }) %>% bind_rows()
  
  # ---------- categorical summaries ----------
  cat_tbl <- lapply(cat_vars, function(v) {
    lab <- var_labels[[v]] %||% v
    
    # denominators for percents: non-missing within each column
    denom_overall <- sum(!is.na(df[[v]]))
    denom_bygrp <- df %>%
      group_by(.group) %>%
      summarise(denom = sum(!is.na(.data[[v]])), .groups = "drop")
    
    # all levels (excluding NA)
    levs <- df[[v]] %>% as.factor() %>% levels()
    
    # counts overall
    overall_counts <- df %>%
      filter(!is.na(.data[[v]])) %>%
      count(.data[[v]], name = "n") %>%
      rename(Level = 1)
    
    # counts by group
    grp_counts <- df %>%
      filter(!is.na(.group)) %>%
      count(.group, .data[[v]], name = "n") %>%
      rename(Level = 2)
    
    # base rows for each level
    out <- tibble(
      Variable = c(lab, rep("", length(levs) - 1), ""),
      Level = c(levs, "Missing")
    )
    
    # Overall column
    out$Overall <- c(
      sapply(levs, function(L) {
        nL <- overall_counts$n[overall_counts$Level == L]
        nL <- ifelse(length(nL) == 0, 0, nL)
        fmt_n_pct(nL, denom_overall)
      }),
      fmt_missing(df[[v]])
    )
    
    # group columns
    for (g in levels(df$.group)) {
      denom_g <- denom_bygrp$denom[denom_bygrp$.group == g]
      denom_g <- ifelse(length(denom_g) == 0, NA_integer_, denom_g)
      
      out[[as.character(g)]] <- c(
        sapply(levs, function(L) {
          nL <- grp_counts$n[grp_counts$.group == g & grp_counts$Level == L]
          nL <- ifelse(length(nL) == 0, 0, nL)
          fmt_n_pct(nL, denom_g)
        }),
        fmt_missing(df[[v]][df$.group == g])
      )
    }
    
    out
  }) %>% bind_rows()
  
  # ---------- sample size row ----------
  Ns <- tibble(
    Variable = "Sample size",
    Level = "N",
    Overall = as.character(n_overall)
  )
  for (g in levels(df$.group)) {
    Ns[[as.character(g)]] <- as.character(sum(df$.group == g, na.rm = TRUE))
  }
  
  # combine
  bind_rows(Ns, cont_tbl, cat_tbl)
}

tab1 <- make_table1(
  hivp1_data,
  group_var = hard_drugs_cat_0,   # baseline hard drug use category column after pivot
  cont_vars = c("BMI_0","age_0", "LEU3N_0", "logVLOAD_0", "AGG_MENT_0", "AGG_PHYS_0"),
  cat_vars  = c("smoke_cat_0","adh_cat_2","race_cat_0","edu_cat_0")
)

# Nice column names
colnames(tab1) <- str_replace_all(colnames(tab1), "^No$", "No hard drugs")
colnames(tab1) <- str_replace_all(colnames(tab1), "^Yes$", "Yes hard drugs")

#code that will be used to knit table1 to pdf from rmarkdown in chunk
table1 <- kbl(
  tab1,
  format = "latex",
  booktabs = TRUE,
  longtable = TRUE,
  align = "llccc",
  caption = "Baseline characteristics stratified by baseline hard drug use."
  ,
  escape = FALSE,
  linesep = "") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), font_size = 8) %>%
  column_spec(1, width = "4.0cm") %>%
  column_spec(2, width = "3.0cm") %>%
  column_spec(3:ncol(tab1), width = "2.2cm")