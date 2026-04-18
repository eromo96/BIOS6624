library(survival)
library(survminer)

#read in code for cleaning script so that we dont have to redo data cleaning here
source("./Project 3/Code/Cleaning.R")
#==============
#Cox Proportional Hazards Modeleling
#==============
#fit full models for each gender separately
fullcox_fem <- coxph(Surv(stroke_time_yrs, stroke_10) ~ AGE + diabetes +
                       SYSBP + prevchd + bpmeds + cursmoke + TOTCHOL + BMI,
                     data = fram_mod_female)
fullcox_male <- coxph(Surv(stroke_time_yrs, stroke_10) ~ AGE + diabetes +
                        SYSBP + prevchd + bpmeds + cursmoke + TOTCHOL + BMI,
                      data = fram_mod_male)
#run backward selection via AIC for each gender separately
back1_fem <- step(fullcox_fem, direction = "backward", 
                  scope = list(lower = ~ AGE + diabetes + SYSBP))
back1_male <- step(fullcox_male, direction = "backward", 
                   scope = list(lower = ~ AGE + diabetes + SYSBP))
#females final model did not select any risk factors of interest
#males final model only selected cursmoker
#BUT females model with cursmoker is only lower by .46 AIC so 
#final models for both genders will be the same to make comparisons easier
#in future inferential work
final_mod_fem <- coxph(Surv(stroke_time_yrs, stroke_10) ~ AGE + diabetes +
                         SYSBP + cursmoke, data = fram_mod_female)
final_mod_male <- coxph(Surv(stroke_time_yrs, stroke_10) ~ AGE + diabetes +
                          SYSBP + cursmoke, data = fram_mod_male)

#==============
#Proportion Hazards Check
#==============
#check proportional hazards assumption for each final model for females
ph_res_fem <- cox.zph(final_mod_fem)
ggcoxzph(ph_res_fem, se = F)
#check proportional hazards assumption for each final model for males
ph_res_male <- cox.zph(final_mod_male)
ggcoxzph(ph_res_male, se = F)

#==============
#10 year Probability Risk Profiles
#==============
#Risk profile 1 values: baseline average
mean_sysbp_f <- mean(fram_mod_female$SYSBP, na.rm = TRUE)
mean_sysbp_m <- mean(fram_mod_male$SYSBP, na.rm = TRUE)

#Female risk profiles table
#create empty results
res_fem <- matrix(NA, nrow = 5, ncol = 3)
rownames(res_fem) <- c(
  "Baseline",
  "High BP (>=160)",
  "Diabetes",
  "High BP + Diabetes",
  "Current Smoker"
)

colnames(res_fem) <- c("Age 40", "Age 50", "Age 60")

ages <- c(40, 50, 60)

for (i in 1:3) {
  age_val <- ages[i]
  
  profiles <- list(
    c(age_val, 0, mean_sysbp_f, 0), #baseline
    c(age_val, 0, 160, 0), #high bp
    c(age_val, 1, mean_sysbp_f, 0), #diabetes
    c(age_val, 1, 160, 0), #hpb + diabetes
    c(age_val, 0, mean_sysbp_f, 1) #smoker
  )
  for (j in 1:5) {
    newdata <- data.frame(
      AGE = profiles[[j]][1],
      diabetes = factor(profiles[[j]][2], levels = c(0,1), labels = c("No","Yes")),
      SYSBP = profiles[[j]][3],
      cursmoke = factor(profiles[[j]][4], levels = c(0,1), labels = c("No","Yes"))
    )
    fit <- survfit(final_mod_fem, newdata = newdata)
    surv_10 <- summary(fit, times = 10)$surv
    res_fem[j, i] <- 1 - surv_10
  }
}

res_fem_df <- as.data.frame(res_fem)
res_fem_df <- round(res_fem_df, 3)


#Male risk profiles table
#create empty results
res_male <- matrix(NA, nrow = 5, ncol = 3)
rownames(res_male) <- c(
  "Baseline",
  "High BP (>=160)",
  "Diabetes",
  "High BP + Diabetes",
  "Current Smoker"
)

colnames(res_male) <- c("Age 40", "Age 50", "Age 60")

for (i in 1:3) {
  age_val <- ages[i]
  
  profiles <- list(
    c(age_val, 0, mean_sysbp_m, 0), #baseline
    c(age_val, 0, 160, 0), #high bp
    c(age_val, 1, mean_sysbp_m, 0), #diabetes
    c(age_val, 1, 160, 0), #hpb + diabetes
    c(age_val, 0, mean_sysbp_m, 1) #smoker
  )
  for (j in 1:5) {
    newdata <- data.frame(
      AGE = profiles[[j]][1],
      diabetes = factor(profiles[[j]][2], levels = c(0,1), labels = c("No","Yes")),
      SYSBP = profiles[[j]][3],
      cursmoke = factor(profiles[[j]][4], levels = c(0,1), labels = c("No","Yes"))
    )
    fit <- survfit(final_mod_male, newdata = newdata)
    surv_10 <- summary(fit, times = 10)$surv
    res_male[j, i] <- 1 - surv_10
  }
}

res_male_df <- as.data.frame(res_male)
res_male_df <- round(res_male_df, 3)



#==============
#Change in Risk factors over 3 periods
#==============

# diabetes by period: n (%)
diab_tab <- fram_long_desc %>%
  group_by(period) %>%
  summarise(
    diab_n = sum(diabetes == "Yes", na.rm = TRUE),
    diab_denom = sum(!is.na(diabetes)),
    diabetes_stat = sprintf("%d (%.1f\\%%)", diab_n, 100 * diab_n / diab_denom)
  ) %>%
  select(period, diabetes_stat)

# systolic blood pressure by period: mean (SD)
sysbp_tab <- fram_long_desc %>%
  group_by(period) %>%
  summarise(
    sysbp_stat = sprintf("%.1f (%.1f)",
                         mean(SYSBP, na.rm = TRUE),
                         sd(SYSBP, na.rm = TRUE))
  ) %>%
  select(period, sysbp_stat)

# combine into one table
change_tab <- data.frame(
  `Risk Factor` = c("Diabetes, n (\\%)", "Systolic blood pressure, mean (SD)"),
  `Exam 1` = c(
    diab_tab$diabetes_stat[diab_tab$period == "Exam 1"],
    sysbp_tab$sysbp_stat[sysbp_tab$period == "Exam 1"]
  ),
  `Exam 2` = c(
    diab_tab$diabetes_stat[diab_tab$period == "Exam 2"],
    sysbp_tab$sysbp_stat[sysbp_tab$period == "Exam 2"]
  ),
  `Exam 3` = c(
    diab_tab$diabetes_stat[diab_tab$period == "Exam 3"],
    sysbp_tab$sysbp_stat[sysbp_tab$period == "Exam 3"]
  ),
  check.names = FALSE
)
