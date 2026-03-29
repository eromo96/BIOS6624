#Code to review/investigate preliminary data

library(tidyverse)
library(car)
#read in preliminary data
prelim_data <- read_csv("./Project 2/Data/PrelimData.csv")

#run basic summary of prelim data to get basic descriptive statistics
summary(prelim_data)

#run mlr models between predictors and each outcome separately to get R^2
cort_mod <- lm(CORT_CNG3 ~ IL_6 + MCP_1, data = prelim_data)
summary(cort_mod)
vif(cort_mod)
cvlt_mod <- lm(CVLT_CNG3 ~ IL_6 + MCP_1, data = prelim_data)
summary(cvlt_mod)
vif(cvlt_mod)
# Quick correlation matrix
cor_mat <- cor(prelim_data, use = "complete.obs")
print(round(cor_mat, 3))
