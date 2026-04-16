#Code for cleaning framingham data
library(tidyverse)

#read in data
fram_data <- read_csv("./Project 3/Data/frmgham2.csv")

#select columns that are of interest
fram_data <- fram_data %>%
  select(
    RANDID, PERIOD, TIME, SEX,
    AGE, DIABETES, SYSBP,
    PREVCHD, BPMEDS, CURSMOKE, TOTCHOL, BMI,
    PREVSTRK, STROKE, TIMESTRK
  )
#recode variables for readibility
fram_data <- fram_data %>%
  mutate(
    sex = factor(SEX, levels = c(1,2), labels = c("Male", "Female")),
    diabetes = factor(DIABETES, levels = c(0,1), labels = c("No", "Yes")),
    bpmeds = factor(BPMEDS, levels = c(0,1), labels = c("No", "Yes")),
    cursmoke = factor(CURSMOKE, levels = c(0, 1), labels = c("No", "Yes")),
    prevchd = factor(PREVCHD, levels = c(0,1), labels = c("No", "Yes")),
    prevstrk = factor(PREVSTRK, levels = c(0,1), labels = c("No", "Yes"))
  )

#create main analysis dataset
fram_base <- fram_data %>%
  #first exam period only and exclude those with stroke at baseline
  filter(PERIOD == 1, PREVSTRK == 0) %>%
  mutate(
    ten_yrdys = 10*365.25,
    #event within first 10 years
    stroke_10 = if_else(STROKE == 1 & TIMESTRK <= ten_yrdys, STROKE, 0),
    #truncate follow up time at 10 years
    stroke_time_days = pmin(TIMESTRK, ten_yrdys),
    stroke_time_yrs = stroke_time_days/365.25
  )
#create separate baseline datasets by gender for data exploration
fram_base_male <- fram_base %>%
  filter(sex == "Male")
fram_base_female <- fram_base %>%
  filter(sex == "Female")
#complete case analysis dataset for aim 1
fram_model <- fram_base %>%
  select(
    RANDID, sex, AGE, diabetes, SYSBP,
    prevchd, bpmeds, cursmoke, TOTCHOL, BMI,
    stroke_10, stroke_time_days, stroke_time_yrs
  ) %>%
  drop_na(AGE, diabetes, SYSBP, prevchd, bpmeds, cursmoke, TOTCHOL, BMI)

#create separate complete case analysis dataset for male and female
fram_mod_male <- fram_model %>%
  filter(sex == "Male")
fram_mod_female <- fram_model %>%
  filter(sex == "Female")
#create dataset for descriptive time varying covariate part
fram_long_desc <- fram_data %>%
  select(RANDID, PERIOD, TIME, sex, AGE, diabetes, SYSBP) %>%
  mutate(
    time_yrs = TIME/365.25,
    period = factor(PERIOD, levels = c(1, 2, 3),
                    labels = c("Exam 1", "Exam 2", "Exam 3"))
  )

