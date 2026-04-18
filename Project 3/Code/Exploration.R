#data exploration and descriptives
library(kableExtra)
library(survival)
library(survminer)
#read in code form cleaning script so that we dont have to redo data cleaning here
source("./Project 3/Code/Cleaning.R")

#=============
#Basic descriptive data exploration
#=============
#some basic data checks/descriptives
fram_base %>% count(RANDID) %>% count(n)
#sample size by sex
fram_base %>% count(sex)
#number of 10 year stroke events by sex
fram_base %>% count(sex, stroke_10)
##missingness
#count and percent of missingness
fram_base %>%
  summarise(
    across(
      c(AGE, DIABETES, SYSBP, PREVCHD, BPMEDS, CURSMOKE, TOTCHOL, BMI),
      ~ paste0(
        sum(is.na(.)), " (",
        round(100 * mean(is.na(.)), 1), "%)"
      )
    )
  )
#pct missing of total dataset
fram_base %>%
  mutate(complete_case = if_all(
    c(AGE, DIABETES, SYSBP, PREVCHD, BPMEDS, CURSMOKE, TOTCHOL, BMI),
    ~ !is.na(.)
  )) %>%
  summarise(
    complete = sum(complete_case),
    missing = sum(!complete_case),
    pct_missing = round(100 * mean(!complete_case), 1)
  )
#missinness by sex
fram_base %>%
  group_by(sex) %>%
  summarise(
    across(
      c(AGE, DIABETES, SYSBP, PREVCHD, BPMEDS, CURSMOKE, TOTCHOL, BMI),
      ~ round(100 * mean(is.na(.)), 1)
    )
  )
#label stroke variable
fram_base <- fram_base %>%
  mutate(
    stroke_10_fac = factor(stroke_10, levels = c(0, 1), labels = c("No", "Yes"))
  )

##create grid (2x5) of kaplan meier curves stratified by gender and 5 risk factors of interest
#create separate dataset of plotting variables
km_data <- fram_base %>%
  mutate(
    bmi_group = factor(
      if_else(BMI >= 30, "Obese (>=30)", "Not Obese (<30)"),
      levels = c("Not Obese (<30)", "Obese (>=30)")
    ),
    chol_group = factor(
      if_else(TOTCHOL >= 240, "High (>=240)", "Not high (<240)"),
      levels = c("Not high (<240)", "High (>=240)")
    )
  ) %>%
  select(
    RANDID, sex, stroke_time_yrs, stroke_10, prevchd, bpmeds,
    cursmoke, chol_group, bmi_group
  )
#reshape to long
km_long <- km_data %>%
  pivot_longer(
    cols = c(prevchd, bpmeds, cursmoke, chol_group, bmi_group),
    names_to = "risk_factor",
    values_to = "group"
  ) %>%
  mutate(
    risk_factor = recode(
      risk_factor,
      prevchd = "Coronary heart disease",
      bpmeds = "Blood pressure meds",
      cursmoke = "Current smoker",
      chol_group = "Total cholesterol",
      bmi_group = "BMI"
    ),
    risk_factor = factor(
      risk_factor,
      levels = c(
        "Coronary heart disease",
        "Blood pressure meds",
        "Current smoker",
        "Total cholesterol",
        "BMI"
      )
    )
  )
#function for survfit output
tidy_survfit <- function(data) {
  fit <- survfit(Surv(stroke_time_yrs, stroke_10) ~ group, data = data)
  
  s <- summary(fit)
  
  out <- tibble(
    time = s$time,
    surv = s$surv,
    strata = s$strata
  )
  
  # extract group name from strings like "group=Yes"
  out <- out %>%
    mutate(group = sub("^group=", "", strata))
  
  out
}
#fit km curves for each sex x risk factor
km_plot_data <- km_long %>%
  group_by(sex, risk_factor) %>%
  nest() %>%
  mutate(km = map(data, tidy_survfit)) %>%
  select(-data) %>%
  unnest(km)
#plot 2x5 faceted km figure
ggplot(km_plot_data, aes(x = time, y = surv, linetype = group)) +
  geom_step(linewidth = 0.6) +
  facet_grid(sex ~ risk_factor) +
  labs(
    x = "Follow-up time (years)",
    y = "Survival probability",
    linetype = NULL,
    title = "Kaplan-Meier curves for stroke-free survival by sex and baseline risk factors"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

#=============
#Table 1
#=============

# Sample sizes for headers
# sample sizes for column headers
n_male <- sum(fram_base$sex == "Male")
n_female <- sum(fram_base$sex == "Female")
n_overall <- nrow(fram_base)

male_header <- paste0("Male\nN=", n_male)
female_header <- paste0("Female\nN=", n_female)
overall_header <- paste0("Overall\nN=", n_overall)

# formatting helpers
fmt_mean_sd <- function(x) {
  if (all(is.na(x))) return("")
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_n_pct <- function(n, denom) {
  if (is.na(denom) || denom == 0) return("")
  sprintf("%d (%.1f\\%%)", n, 100 * n / denom)
}

# continuous block
make_cont_block <- function(data, var, label) {
  x_all <- data[[var]]
  x_m <- data %>% filter(sex == "Male") %>% pull(all_of(var))
  x_f <- data %>% filter(sex == "Female") %>% pull(all_of(var))
  
  out <- tibble(
    Variable = label,
    Male = fmt_mean_sd(x_m),
    Female = fmt_mean_sd(x_f),
    Overall = fmt_mean_sd(x_all)
  )
  
  miss_all <- sum(is.na(x_all))
  miss_m <- sum(is.na(x_m))
  miss_f <- sum(is.na(x_f))
  
  if (miss_all > 0) {
    out <- bind_rows(
      out,
      tibble(
        Variable = "\\hspace{1em}Missing, n (\\%)",
        Male = fmt_n_pct(miss_m, length(x_m)),
        Female = fmt_n_pct(miss_f, length(x_f)),
        Overall = fmt_n_pct(miss_all, length(x_all))
      )
    )
  }
  
  out
}

# categorical block
make_cat_block <- function(data, var, label) {
  x_all <- data[[var]]
  x_m <- data %>% filter(sex == "Male") %>% pull(all_of(var))
  x_f <- data %>% filter(sex == "Female") %>% pull(all_of(var))
  
  levs <- levels(droplevels(factor(x_all)))
  
  rows <- lapply(levs, function(lv) {
    denom_all <- sum(!is.na(x_all))
    denom_m <- sum(!is.na(x_m))
    denom_f <- sum(!is.na(x_f))
    
    n_all <- sum(x_all == lv, na.rm = TRUE)
    n_m <- sum(x_m == lv, na.rm = TRUE)
    n_f <- sum(x_f == lv, na.rm = TRUE)
    
    tibble(
      Variable = paste0("\\hspace{1em}", lv),
      Male = fmt_n_pct(n_m, denom_m),
      Female = fmt_n_pct(n_f, denom_f),
      Overall = fmt_n_pct(n_all, denom_all)
    )
  })
  
  out <- bind_rows(
    tibble(
      Variable = label,
      Male = "",
      Female = "",
      Overall = ""
    ),
    bind_rows(rows)
  )
  
  miss_all <- sum(is.na(x_all))
  miss_m <- sum(is.na(x_m))
  miss_f <- sum(is.na(x_f))
  
  if (miss_all > 0) {
    out <- bind_rows(
      out,
      tibble(
        Variable = "\\hspace{1em}Missing, n (\\%)",
        Male = fmt_n_pct(miss_m, length(x_m)),
        Female = fmt_n_pct(miss_f, length(x_f)),
        Overall = fmt_n_pct(miss_all, length(x_all))
      )
    )
  }
  
  out
}

# build table
table1 <- bind_rows(
  make_cat_block(fram_base, "stroke_10_fac", "Stroke during follow-up"),
  make_cont_block(fram_base, "stroke_time_yrs", "Follow-up time (years)"),
  make_cont_block(fram_base, "AGE", "Age"),
  make_cont_block(fram_base, "SYSBP", "Systolic blood pressure"),
  make_cat_block(fram_base, "bpmeds", "Blood pressure meds"),
  make_cont_block(fram_base, "TOTCHOL", "Total cholesterol"),
  make_cont_block(fram_base, "BMI", "BMI"),
  make_cat_block(fram_base, "diabetes", "Diabetes"),
  make_cat_block(fram_base, "prevchd", "Coronary heart disease"),
  make_cat_block(fram_base, "cursmoke", "Current smoker")
)

colnames(table1) <- c(
  "Variable",
  male_header,
  female_header,
  overall_header
)
