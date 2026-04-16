library(survival)
library(survminer)

#read in code form cleaning script so that we dont have to redo data cleaning here
source("./Project 3/Code/Cleaning.R")

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
#in futrue inferential work
final_mod_fem <- coxph(Surv(stroke_time_yrs, stroke_10) ~ AGE + diabetes +
                         SYSBP + cursmoke, data = fram_mod_female)
final_mod_male <- coxph(Surv(stroke_time_yrs, stroke_10) ~ AGE + diabetes +
                          SYSBP + cursmoke, data = fram_mod_male)
