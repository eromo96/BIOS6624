#data cleaning for project0 data
library(tidyverse)

#read in data
p0_data <- read_csv("./Project0/Data/Project0_Clean_v2.csv")

#remove columns that wont be used for data analysis
p0_data <- p0_data %>%
  select(-`Booklet: Sample interval`, -`Booklet: Sample interval Decimal Time (mins)`, -`MEMs: Sample interval`,
         -`MEMs: Sample interval Decimal Time (mins)`, -`Cortisol (ug/dl)`, -`DHEA (pg/dl)`)

#convert date columns to appropriate date format
p0_data <- p0_data %>%
  mutate(
    `Collection Date` = mdy(`Collection Date`)
  ) %>%
  rename(
    date = `Collection Date`,
    sample = `Collection Sample`,
    wake_tm = `Sleep Diary reported wake time`,
    booklet_tm = `Booket: Clock Time`,
    mems_tm = `MEMs: Clock Time`,
    cortisol = `Cortisol (nmol/L)`,
    dhea = `DHEA (nmol/L)`,
    day = DAYNUMB
  )

#explore missingness of data
miss_p0 <- p0_data %>%
  summarise(across(everything(),
                   ~mean(is.na(.)))) %>%
  pivot_longer(everything(),
               names_to = "variable",
               values_to = "prop_missing") %>%
  arrange(desc(prop_missing))
print(miss_p0)  

#manually calculate sample interval times for MEMs and Booklet 
p0_data <- p0_data %>%
  group_by(SubjectID, date) %>%
  arrange(sample, .by_group = TRUE) %>%
  fill(wake_tm, .direction = "down") %>%
  ungroup() %>%
  mutate(
    mems_interval = as.numeric(mems_tm - wake_tm)/3600,
    booklet_interval = as.numeric(booklet_tm - wake_tm)/3600
  )

#create truth columns of +/- 7.5 and 15 mins to check adherence
p0_data <- p0_data %>%
  mutate(
    truth = case_when(
      sample == 2 ~ .5,
      sample == 4 ~ 10,
      .default = NA
    ),
    deviation_mems = abs(mems_interval - truth),
    deviation_book = abs(booklet_interval - truth)
  )

#check adherence
p0_data <- p0_data %>%
  mutate(
    adherence_mems7_5 = if_else(deviation_mems <= 7.5/60, 1, 0),
    adherence_book7_5 = if_else(deviation_book <= 7.5/60, 1, 0),
    adherence_mems15 = if_else(deviation_mems <= 15/60, 1, 0),
    adherence_book15 = if_else(deviation_book <= 15/60, 1, 0)
  )

#calc adherence proportions
adherence <- p0_data %>%
  filter(!is.na(truth)) %>%
  summarise(
    within_7_5_mems = mean(adherence_mems7_5, na.rm = TRUE),
    within_15_mems = mean(adherence_mems15, na.rm = TRUE),
    within_7_5_book = mean(adherence_book7_5, na.rm = TRUE),
    within_15_book = mean(adherence_book15, na.rm = TRUE)
  )


