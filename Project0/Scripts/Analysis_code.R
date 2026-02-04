#data cleaning for project0 data
library(kableExtra)
library(tidyverse)
library(lme4)
library(lmerTest)
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



##########--------Table 1------#############

#create table 1 of summary statistics for difference between mems and booklet
#time, cortisol, and dhea stratified by sample

# ---- Prep analysis dataset ----
table1_data <- p0_data %>%
  mutate(
    Sample = factor(sample,
                    levels = 1:4,
                    labels = c("Sample 1\n(Wake)",
                               "Sample 2\n(+30 min)",
                               "Sample 3\n(Lunch)",
                               "Sample 4\n(+10 hr)")
    ),
    time_diff = mems_interval - booklet_interval
  )

# ---- Helper formatters ----
fmt_mean_sd <- function(x, digits = 2){
  if(all(is.na(x))) return(NA_character_)
  sprintf(paste0("%.",digits,"f (%.",digits,"f)"),
          mean(x, na.rm = TRUE),
          sd(x, na.rm = TRUE))
}

fmt_median_iqr <- function(x, digits = 2){
  if(all(is.na(x))) return(NA_character_)
  q <- quantile(x, c(.25,.5,.75), na.rm = TRUE)
  sprintf(paste0("%.",digits,"f [%.",digits,"f, %.",digits,"f]"),
          q[2], q[1], q[3])
}

fmt_missing <- function(x){
  n_miss <- sum(is.na(x))
  p_miss <- mean(is.na(x)) * 100
  sprintf("%d (%.1f%%)", n_miss, p_miss)
}

# ---- By-sample summaries (LONG) ----
table1_long <- table1_data %>%
  group_by(Sample) %>%
  summarise(
    N = as.character(n()),
    
    diff_mean_sd = fmt_mean_sd(time_diff),
    diff_med_iqr = fmt_median_iqr(time_diff),
    diff_miss    = fmt_missing(time_diff),
    
    cort_mean_sd = fmt_mean_sd(cortisol),
    cort_med_iqr = fmt_median_iqr(cortisol),
    cort_miss    = fmt_missing(cortisol),
    
    dhea_mean_sd = fmt_mean_sd(dhea),
    dhea_med_iqr = fmt_median_iqr(dhea),
    dhea_miss    = fmt_missing(dhea),
    
    .groups = "drop"
  ) %>%
  pivot_longer(-Sample, names_to = "Statistic", values_to = "Value")

# ---- Overall summaries (LONG) ----
overall_long <- table1_data %>%
  summarise(
    N = as.character(n()),
    
    diff_mean_sd = fmt_mean_sd(time_diff),
    diff_med_iqr = fmt_median_iqr(time_diff),
    diff_miss    = fmt_missing(time_diff),
    
    cort_mean_sd = fmt_mean_sd(cortisol),
    cort_med_iqr = fmt_median_iqr(cortisol),
    cort_miss    = fmt_missing(cortisol),
    
    dhea_mean_sd = fmt_mean_sd(dhea),
    dhea_med_iqr = fmt_median_iqr(dhea),
    dhea_miss    = fmt_missing(dhea)
  ) %>%
  pivot_longer(everything(), names_to = "Statistic", values_to = "Value") %>%
  mutate(Sample = "Overall")

# ---- Combine + order + relabel stats + flip to wide ----
stat_order <- c(
  "N",
  "diff_mean_sd", "diff_med_iqr", "diff_miss",
  "cort_mean_sd", "cort_med_iqr", "cort_miss",
  "dhea_mean_sd", "dhea_med_iqr", "dhea_miss"
)

table1_final <- bind_rows(table1_long, overall_long) %>%
  mutate(
    Statistic = factor(Statistic, levels = stat_order),
    Statistic = recode(as.character(Statistic),
                       N = "N",
                       
                       diff_mean_sd = "Time diff (MEMS − Booklet), mean (SD) [hours]",
                       diff_med_iqr = "Time diff (MEMS − Booklet), median [IQR] [hours]",
                       diff_miss    = "Time diff missing, n (%)",
                       
                       cort_mean_sd = "Cortisol (nmol/L), mean (SD)",
                       cort_med_iqr = "Cortisol (nmol/L), median [IQR]",
                       cort_miss    = "Cortisol missing, n (%)",
                       
                       dhea_mean_sd = "DHEA (nmol/L), mean (SD)",
                       dhea_med_iqr = "DHEA (nmol/L), median [IQR]",
                       dhea_miss    = "DHEA missing, n (%)"
    )
  ) %>%
  arrange(match(Statistic, levels(factor(Statistic)))) %>%
  pivot_wider(names_from = Sample, values_from = Value)

# ---- Render to PDF using kableExtra ----
kable(
  table1_final,
  format = "latex",
  booktabs = TRUE,
  align = "lccccc",
  caption = "Table 1. Descriptive statistics by sample collection time"
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 10
  )





##########--------QUESTION 1------#############

#run linear mixed effects model of booklet time vs mems time with random intercept
q1_mod1 <- lmer(mems_interval ~ booklet_interval + (1|SubjectID), data = p0_data)

#create table of lmm for question 1
#get fixed effects from model
fixef_tab <- summary(q1_mod1)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  select(
    Term,
    Estimate,
    `Std. Error`,
    `Pr(>|t|)`
  )
#95% ci
ci_tab <- confint(q1_mod1) %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  rename(
    LL = `2.5 %`,
    UL = `97.5 %`
  ) %>%
  filter(Term %in% fixef_tab$Term)
#putting table together
model_table <- fixef_tab %>%
  left_join(ci_tab, by = "Term") %>%
  mutate(
    CI_95 = sprintf("[%.3f, %.3f]", LL, UL),
    Estimate = round(Estimate, 3),
    `Std. Error` = round(`Std. Error`, 3),
    `Pr(>|t|)` = signif(`Pr(>|t|)`, 3)
  ) %>%
  select(
    Term,
    Estimate,
    CI_95,
    `Pr(>|t|)`
  ) %>%
  rename(
    `95% CI` = CI_95,
    `p-value` = `Pr(>|t|)`
  )
#extract hypothesis test of whether beta1 is equal to 1
contest_res <- contest(q1_mod1, L = c(0, 1), rhs = 1)

contest_table <- tibble(
  Term = "Slope test: β1 = 1",
  Estimate = NA,
  `95% CI` = NA,
  `p-value` = signif(contest_res$`Pr(>F)`[1], 3)
)
#final table
final_table <- bind_rows(
  model_table,
  contest_table
)
#change row names so that they are appropriate for final report
rownames(final_table) <- c("Intercept", "Booklet Interval", "Slope Test = 1")
#use kbl to knit to pdf nicely (kbl is the best)
final_table %>%
  kbl(
    caption = "Linear mixed model results for MEMS interval vs booklet interval",
    align = "lccc",
    booktabs = TRUE,
    linesep = ""
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 10
  )

#scatterplot of agreement between mems and booklet times
ggplot(p0_data, aes(x= booklet_interval, y= mems_interval)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_abline(slope = fixef(q1_mod1)["booklet_interval"], intercept = fixef(q1_mod1)["(Intercept)"],
              linewidth = 1, col = "blue") +
  labs(title = "Agreement Between MEMs and Booklet Sampling Times",
       x = "Booklet Time from Wake", y= "MEMs Time from Wake")


##########--------QUESTION 2------#############

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

#transform adherence table to long so that it's ready to print
adherence_table <- adherence %>%
  pivot_longer(
    cols = everything(),
    names_to = "name",
    values_to = "prop"
  ) %>%
  mutate(
    Method = if_else(grepl("_mems$", name), "MEMS", "Booklet"),
    Window = case_when(
      grepl("^within_7_5_", name) ~ "±7.5 minutes",
      grepl("^within_15_",  name) ~ "±15 minutes"
    )
  ) %>%
  select(Window, Method, prop) %>%
  pivot_wider(names_from = Method, values_from = prop)
#convert to percents
adherence_table <- adherence_table %>%
  mutate(
    MEMS = scales::percent(MEMS, accuracy = 0.1),
    Booklet = scales::percent(Booklet, accuracy = 0.1)
  )
#use kbl to print adherence table so that it is formatted nicely
adherence_table %>%
  kbl(
    caption = "Adherence of MEMS and Booklet Sampling Times Within Specified Time Windows",
    align = c("l", "c", "c"),
    booktabs = TRUE,
    linesep = ""
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 10
  )

##########--------QUESTION 3------#############

#check to see subjects with DHEA values equal to 5.205
p0_data %>% filter(dhea == 5.205)

#subject 3037 has several DHEA measures at detection limit so remove from Q3 analysis
q3_data <- p0_data %>%
  filter(SubjectID != 3037)

#cortisol data filtering out samples with cortisol values over 80 due to lab error
cortisol_data <- q3_data %>%
  filter(cortisol <=80)

#dhea data filtering out samples with dhea == 5.205
dhea_data <- q3_data %>%
  filter(dhea < 5.205)

#check to see distribution of cortisol
hist(cortisol_data$cortisol) #found cortisol to be heavily right skewed and not normal
#log transform cortisol and check distribution
hist(log(cortisol_data$cortisol)) #appears more normally distributed

#boxplot to further visualize outliers
ggplot(cortisol_data, aes(x = factor(sample), y = cortisol)) +
  geom_boxplot() +
  labs(x="Sample (1–4)", y="Cortisol (nmol/L)", title="Cortisol by sample")

#create piecewise at 30 mins for mems time
cortisol_data <- cortisol_data %>%
  mutate(
    booklet_step = pmax(0, booklet_interval-0.5)
  )

#plot spaghetti plots of cortisol trajectories for each subject
ggplot(cortisol_data, aes(x = sample, y = cortisol, group = SubjectID)) + 
  geom_line(alpha = 0.6, color = "gray40") +
  theme_minimal()

#distribution of cortisol is not normal whereas log(x) is approx normal so use log transformed
q3_cortmod1 <- lmer(log(cortisol) ~ booklet_interval + booklet_step + (1|SubjectID),
                    data = cortisol_data)
summary(q3_cortmod1)

#confidence interval
cort_conf <- as.data.frame(confint(q3_cortmod1)[3:5,]) 

#create table for model results
cortmod_tab <- as.data.frame(coef(summary(q3_cortmod1))) %>%
  select(Estimate, "Pr(>|t|)")%>%
  mutate("95% CI" = paste0("[", round(exp(cort_conf$`2.5 %`), 2), ", ", 
                        round(exp(cort_conf$`97.5 %`), 2), "]"),
         "exp(Estimate)" = exp(Estimate),
         "P-value" = round(`Pr(>|t|)`, 4)) %>%
  select("exp(Estimate)", "95% CI", "P-value")

#cortisol prediction grid
cort_grid <- data.frame(
  booklet_interval = seq(0, 12, length = 200)
)

cort_grid$booklet_step <- pmax(0, cort_grid$booklet_interval - 0.5)

cort_grid$pred <- predict(q3_cortmod1,
                          newdata = cort_grid,
                          re.form = NA)
#plot of predicted over observed
ggplot(cortisol_data,
       aes(booklet_interval, log(cortisol))) +
  geom_point(alpha = 0.4) +
  geom_line(data = cort_grid,
            aes(booklet_interval, pred),
            linewidth = 1.2) +
  labs(
    x = "Time Since Wake (hours)",
    y = "log(Cortisol)",
    title = "Predicted Cortisol Diurnal Curve"
  )


#check to see distribution of dhea
hist(dhea_data$dhea) #not normal and heavily right skewed
#log transform dhea and check distribution
hist(log(dhea_data$dhea)) #still right skewed but a bit more normal

#boxplot to visualize outliers
ggplot(dhea_data, aes(x = factor(sample), y = dhea)) +
  geom_boxplot() +
  labs(x="Sample (1–4)", y="DHEA (nmol/L)", title="DHEA by sample")

#create piecewise at 30 mins for mems time
dhea_data <- dhea_data %>%
  mutate(
    booklet_step = pmax(0, booklet_interval-0.5)
  )

#spaghetti plot of subject's dhea trajectories over time
ggplot(dhea_data, aes(x = sample, y = dhea, group = SubjectID)) + 
  geom_line(alpha = 0.6, color = "gray40") +
  theme_minimal()

#distribution of dhea is not normal whereas log(x+1) is a bit better
q3_dheamod1 <- lmer(log(dhea) ~ booklet_interval + booklet_step + (1|SubjectID),
                    data = dhea_data)
summary(q3_dheamod1)

#confidence interval
dhea_conf <- as.data.frame(confint(q3_dheamod1)[3:5,]) 

#create table for model results
dheamod_tab <- as.data.frame(coef(summary(q3_dheamod1))) %>%
  select(Estimate, "Pr(>|t|)")%>%
  mutate("95% CI" = paste0("[", round(exp(dhea_conf$`2.5 %`), 2), ", ", 
                        round(exp(dhea_conf$`97.5 %`), 2), "]"),
         "exp(Estimate)" = exp(Estimate),
         "P-value" = round(`Pr(>|t|)`, 4)) %>%
  select("exp(Estimate)", "95% CI", "P-value")

#dhea predicted grid
dhea_grid <- data.frame(
  booklet_interval = seq(0, 12, length = 200)
)

dhea_grid$booklet_step <- pmax(0, dhea_grid$booklet_interval - 0.5)

dhea_grid$pred <- predict(q3_dheamod1,
                          newdata = dhea_grid,
                          re.form = NA)
#plot of predicted over observed
ggplot(dhea_data,
       aes(booklet_interval, log(dhea))) +
  geom_point(alpha = 0.4) +
  geom_line(data = dhea_grid,
            aes(booklet_interval, pred),
            linewidth = 1.2) +
  labs(
    x = "Time Since Wake (hours)",
    y = "log(DHEA)",
    title = "Predicted DHEA Diurnal Curve"
  )


#combine cortisol and dhea model outputs into single table
tab_combined <- cbind(
  `Cortisol Model` = cortmod_tab,
  `DHEA Model`     = dheamod_tab
)
#change rownames so they are appropriate for report
rownames(tab_combined) <- c("Intercept", "Booklet Time Pre-30m", "Booklet Time Post-30m")
# Flatten the weird colnames created by cbind() into clean names
colnames(tab_combined) <- c(
  "exp(Estimate)", "95\\% CI", "P-value",
  "exp(Estimate)", "95\\% CI", "P-value"
)

kbl(
  tab_combined,
  booktabs = TRUE,
  align = "lcccccc",
  caption = "Hormone Fluctuation Over Time via Linear Mixed Effects",
  escape = FALSE
) |>
  add_header_above(c(" " = 1, "Cortisol Model" = 3, "DHEA Model" = 3)) |>
  kable_styling(latex_options = c("hold_position", "striped")) |>
  column_spec(1, width = "3.0cm") |>
  column_spec(2:7, width = "2.0cm") |>
  footnote(
    general = "exp(Estimate) represents exponentiated fixed-effect estimates",
    threeparttable = TRUE
  )


##code for creating table of back transformed and interpretable relative change in hormones pre and post 30 mins
############### CORTISOL ################

b_cort <- fixef(q3_cortmod1)
V_cort <- vcov(q3_cortmod1)

b1_c <- b_cort["booklet_interval"]
b2_c <- b_cort["booklet_step"]

# Slopes
slope_c_0_30 <- b1_c
slope_c_30p  <- b1_c + b2_c

# SEs
se_c_0_30 <- sqrt(V_cort["booklet_interval","booklet_interval"])
se_c_30p  <- sqrt(
  V_cort["booklet_interval","booklet_interval"] +
    V_cort["booklet_step","booklet_step"] +
    2*V_cort["booklet_interval","booklet_step"]
)

# df + p-values from model summary
coef_c <- summary(q3_cortmod1)$coefficients
df_c <- coef_c["booklet_interval","df"]
p_c_0_30 <- coef_c["booklet_interval","Pr(>|t|)"]

# Wald p for combined slope
t_c_30p <- slope_c_30p / se_c_30p
p_c_30p <- 2 * pt(abs(t_c_30p), df=df_c, lower.tail=FALSE)

crit_c <- qt(0.975, df_c)

############### DHEA ################

b_dhea <- fixef(q3_dheamod1)
V_dhea <- vcov(q3_dheamod1)

b1_d <- b_dhea["booklet_interval"]
b2_d <- b_dhea["booklet_step"]

slope_d_0_30 <- b1_d
slope_d_30p  <- b1_d + b2_d

se_d_0_30 <- sqrt(V_dhea["booklet_interval","booklet_interval"])
se_d_30p  <- sqrt(
  V_dhea["booklet_interval","booklet_interval"] +
    V_dhea["booklet_step","booklet_step"] +
    2*V_dhea["booklet_interval","booklet_step"]
)

coef_d <- summary(q3_dheamod1)$coefficients
df_d <- coef_d["booklet_interval","df"]
p_d_0_30 <- coef_d["booklet_interval","Pr(>|t|)"]

t_d_30p <- slope_d_30p / se_d_30p
p_d_30p <- 2 * pt(abs(t_d_30p), df=df_d, lower.tail=FALSE)

crit_d <- qt(0.975, df_d)

############### BUILD TABLE ################

final_q3_table <- tibble(
  `Time Window` = c("0–30 min", "30+ min"),
  `Cortisol % Change` = c(
    round((exp(slope_c_0_30 * 0.5) - 1) * 100, 1),   
    round((exp(slope_c_30p) - 1) * 100, 1)           
  ),
  
  `Cortisol 95% CI` = c(
    paste0("[",
           round((exp((slope_c_0_30 - crit_c*se_c_0_30) * 0.5) - 1) * 100, 1),
           ", ",
           round((exp((slope_c_0_30 + crit_c*se_c_0_30) * 0.5) - 1) * 100, 1),
           "]"
    ),
    paste0("[",
           round((exp(slope_c_30p - crit_c*se_c_30p) - 1) * 100, 1),
           ", ",
           round((exp(slope_c_30p + crit_c*se_c_30p) - 1) * 100, 1),
           "]"
    )
  ),
  
  `Cortisol P-value` = ifelse(
    c(p_c_0_30, p_c_30p) < 0.0001,
    "<0.0001",
    formatC(c(p_c_0_30, p_c_30p), format="f", digits=4)
  ),
  
  `DHEA % Change` = c(
    round((exp(slope_d_0_30 * 0.5) - 1) * 100, 1),
    round((exp(slope_d_30p) - 1) * 100, 1)
  ),
  
  `DHEA 95% CI` = c(
    paste0("[",
           round((exp((slope_d_0_30 - crit_d*se_d_0_30) * 0.5) - 1) * 100, 1),
           ", ",
           round((exp((slope_d_0_30 + crit_d*se_d_0_30) * 0.5) - 1) * 100, 1),
           "]"
    ),
    paste0("[",
           round((exp(slope_d_30p - crit_d*se_d_30p) - 1) * 100, 1),
           ", ",
           round((exp(slope_d_30p + crit_d*se_d_30p) - 1) * 100, 1),
           "]"
    )
  ),
  
  `DHEA P-value` = ifelse(
    c(p_d_0_30, p_d_30p) < 0.0001,
    "<0.0001",
    formatC(c(p_d_0_30, p_d_30p), format="f", digits=4)
  )
)


############### PRINT TABLE ################

kbl(
  final_q3_table,
  format="latex",
  booktabs=TRUE,
  align="lccc|ccc",
  caption="Percent change in hormone concentrations"
) %>%
  add_header_above(c(" " = 1, "Cortisol Model" = 3, "DHEA Model" = 3)) %>%
  kable_styling(
    latex_options=c("hold_position","scale_down"),   
    font_size=9                                      
  ) %>%
  column_spec(1, width="2.8cm") %>%
  column_spec(2:7, width="1.8cm")
