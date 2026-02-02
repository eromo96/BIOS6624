#data cleaning for project0 data
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


##########--------QUESTION 1------#############

#create scatterplot of mems vs booklet interval times
ggplot(p0_data, aes(x= booklet_interval, y= mems_interval)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Booklet Time from Wake", y= "MEMs Time from Wake")

#run linear mixed effects model of booklet time vs mems time with random intercept
q1_mod1 <- lmer(mems_interval ~ booklet_interval + (1|SubjectID), data = p0_data)

#output of model; p-value for intercept tests whether its equal to 0 or not(of interest)
summary(q1_mod1)
confint(q1_mod1)

#test whether beta1 is equal to 1 or not (of interest)
contest(q1_mod1, L = c(0, 1), rhs = 1)

#now fit a model of just diff between MEMs and Booklet
q1_mod2 <- lmer(mems_interval - booklet_interval ~ 1 + (1|SubjectID), data = p0_data)

summary(q1_mod2)
confint(q1_mod2)



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
hist(log1p(cortisol_data$cortisol)) #appears more normally distributed

#boxplot to further visualize outliers
ggplot(cortisol_data, aes(x = factor(sample), y = cortisol)) +
  geom_boxplot() +
  labs(x="Sample (1–4)", y="Cortisol (nmol/L)", title="Cortisol by sample")

#create piecewise at 30 mins for mems time
cortisol_data <- cortisol_data %>%
  mutate(
    mems_step = pmax(0, mems_interval-0.5)
  )

#distribution of cortisol is not normal whereas log(x+1) is approx normal so use log transformed
q3_cortmod1 <- lmer(log1p(cortisol) ~ mems_interval + mems_step + (1|SubjectID),
                    data = cortisol_data)
summary(q3_cortmod1)

#cortisol prediction grid
cort_grid <- data.frame(
  mems_interval = seq(0, 12, length = 200)
)

cort_grid$mems_step <- pmax(0, cort_grid$mems_interval - 0.5)

cort_grid$pred <- predict(q3_cortmod1,
                          newdata = cort_grid,
                          re.form = NA)
#plot of predicted over observed
ggplot(cortisol_data,
       aes(mems_interval, log1p(cortisol))) +
  geom_point(alpha = 0.4) +
  geom_line(data = cort_grid,
            aes(mems_interval, pred),
            linewidth = 1.2) +
  labs(
    x = "Time Since Wake (hours)",
    y = "log(1 + Cortisol)",
    title = "Predicted Cortisol Diurnal Curve"
  )


#check to see distribution of dhea
hist(dhea_data$dhea) #not normal and heavily right skewed
#log transform dhea and check distribution
hist(log1p(dhea_data$dhea)) #still right skewed but a bit more normal

#boxplot to visualize outliers
ggplot(dhea_data, aes(x = factor(sample), y = dhea)) +
  geom_boxplot() +
  labs(x="Sample (1–4)", y="DHEA (nmol/L)", title="DHEA by sample")

#create piecewise at 30 mins for mems time
dhea_data <- dhea_data %>%
  mutate(
    mems_step = pmax(0, mems_interval-0.5)
  )

#distribution of dhea is not normal whereas log(x+1) is a bit better
q3_dheamod1 <- lmer(log1p(dhea) ~ mems_interval + mems_step + (1|SubjectID),
                    data = dhea_data)
summary(q3_dheamod1)

#dhea predicted grid
dhea_grid <- data.frame(
  mems_interval = seq(0, 12, length = 200)
)

dhea_grid$mems_step <- pmax(0, dhea_grid$mems_interval - 0.5)

dhea_grid$pred <- predict(q3_dheamod1,
                          newdata = dhea_grid,
                          re.form = NA)
#plot of predicted over observed
ggplot(dhea_data,
       aes(mems_interval, log1p(dhea))) +
  geom_point(alpha = 0.4) +
  geom_line(data = dhea_grid,
            aes(mems_interval, pred),
            linewidth = 1.2) +
  labs(
    x = "Time Since Wake (hours)",
    y = "log(1 + DHEA)",
    title = "Predicted DHEA Diurnal Curve"
  )
