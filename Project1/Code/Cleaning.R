#### CODE FOR DATA CLEANING OF PROJECT1
library(tidyverse)
#read in data
hiv_data <- read_csv("./Project1/Data/hiv_6624_final.csv") %>%
  select(-`...1`)

#see unique race and education levels in data to see which labels are needed
unique(hiv_data$RACE)
unique(hiv_data$EDUCBAS)

#collapse race variable to White vs Other
hiv_data <- hiv_data %>%
  mutate(
    race_cat = factor(case_when(
      RACE %in% c(1, 2) ~ 1,
      RACE %in% c(3, 4, 7, 8) ~ 2
    ),
  levels = c(1, 2), labels = c("White", "Other")))
#collapse education variable to College Degree+ vs Non College Degree
hiv_data <- hiv_data %>%
  mutate(
    edu_cat = factor(case_when(
      EDUCBAS %in% c(1, 2, 3, 4) ~ 1,
      EDUCBAS %in% c(5, 6, 7) ~ 2
    ),
    levels = c(1, 2), labels = c("No College Degree", "College Degree")
    )
  )

#recode abnormal BMI values so that they are missing/NAs
hiv_data <- hiv_data %>%
  mutate(
    BMI = if_else(BMI <= 0 | BMI > 250, NA_real_, BMI)
  )
#create labelled versions of categorical variables
hiv_data <- hiv_data %>%
  mutate(
    smoke_cat = factor(case_when(
      SMOKE %in% c(1,2) ~ 1,
      SMOKE == 3 ~ 2
    ),
      levels = c(1, 2),
      labels = c("Not current smoker",
                 "Current smoker")
    ),
    adh_cat = factor(case_when(
      ADH %in% c(1, 2) ~ 1,
      ADH %in% c(3, 4) ~ 2
    ),
      levels = c(1, 2),
      labels = c(">= 95% Adherence",
                 "< 95% Adherence")
    ),
    hard_drugs_cat = factor(
      hard_drugs,
      levels = c(0, 1),
      labels = c("No",
                 "Yes")
    )
  )


