#data exploration and descriptives
library(kableExtra)

#read in code form cleaning script so that we dont have to redo data cleaning here
source("./Project 3/Code/Cleaning.R")

#some basic data checks/descriptives
fram_base %>% count(RANDID) %>% count(n)
#sample size by sex
fram_base %>% count(sex)
#number of 10 year stroke events by sex
fram_base %>% count(sex, stroke_10)
#missingness in analysis variables
fram_base %>%
  summarise(
    across(c(AGE, DIABETES, SYSBP, PREVCHD, BPMEDS, CURSMOKE, TOTCHOL, BMI),
           ~ sum(is.na(.)))
  )
#label stroke variable
fram_base <- fram_base %>%
  mutate(
    stroke_10 = factor(stroke_10, levels = c(0, 1), labels = c("No", "Yes"))
  )
####Table 1
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
  make_cat_block(fram_base, "stroke_10", "Stroke during follow-up"),
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
