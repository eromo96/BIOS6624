#### CODE FOR DATA CLEANING OF PROJECT1
library(tidyverse)
#read in data
hiv_data <- read_csv("./Project1/Data/hiv_6624_final.csv") %>%
  select(-`...1`)

#see unique race and education levels in data to see which labels are needed
unique(hiv_data$RACE)

#create labelled versions of categorical variables
hiv_data <- hiv_data %>%
  mutate(
    race_cat = factor(
      RACE,
      levels = c(1, 2, 3, 4, 7, 8),
      labels = c("White, non-Hispanic",
                 "White, Hispanic",
                 "Black, non-Hispanic",
                 "Black, Hispanic",
                 "Other",
                 "Other Hispanic")
    ),
    edu_cat = factor(
      EDUCBAS,
      levels = c(1, 2, 3, 4, 5, 6, 7),
      labels = c("8th grade or less",
                 "9, 10, or 11th grade",
                 "12th grade",
                 "At least one year college but no degree",
                 "Four years college/ got degree",
                 "Some graduate work",
                 "Post-graduate degree")
    ),
    smoke_cat = factor(
      SMOKE,
      levels = c(1, 2, 3),
      labels = c("Never smoked",
                 "Former smoker",
                 "Current smoker")
    ),
    adh_cat = factor(
      ADH,
      levels = c(1, 2, 3, 4),
      labels = c("100%",
                 "95-99%",
                 "75-94%",
                 "<75%")
    ),
    hard_drugs_cat = factor(
      hard_drugs,
      levels = c(0, 1),
      labels = c("No",
                 "Yes")
    )
  )


