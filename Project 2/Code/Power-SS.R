library(powertools)
# Power/Sample Size code

####Aim 1 power calculation
N_total <- 175
##Branch A: joint model
#range of Rsquared values based on preliminary data

aim1_joint_grid <- tibble(
  Rsq = c(0.05, 0.1, 0.2, 0.35, 0.5, 0.6)
) %>%
  mutate(
    power_p4 = map_dbl(
      Rsq,
      ~ mlrF.overall(
        N = N_total,
        p = 4,
        Rsq = .x,
        alpha = 0.025
      )
    )
  )

aim1_joint_grid

##Branch B: separate models if predictors highly collinear
aim1_corr_grid <- tibble(
  rhoA = c(0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9)
) %>%
  mutate(
    power = map_dbl(
      rhoA,
      ~ corr.1samp(
        N = N_total,
        rho0 = 0,
        rhoA = .x,
        alpha = 0.025,
        sides = 2
      )
    )
  )

aim1_corr_grid


####Aim 2 power calculation
aim2_power_final <- tribble(
  ~n1, ~n2, ~rho1, ~rho2,
  150,  25,  0.10,  0.20,
  150,  25,  0.10,  0.30,
  150,  25,  0.10,  0.40,
  150,  25,  0.10,  0.50,
  150,  25,  0.10,  0.60,
  125,  50,  0.10,  0.20,
  125,  50,  0.10,  0.30,
  125,  50,  0.10,  0.40,
  125,  50,  0.10,  0.50,
  125,  50,  0.10,  0.60,
  100,  75,  0.10,  0.20,
  100,  75,  0.10,  0.30,
  100,  75,  0.10,  0.40,
  100,  75,  0.10,  0.50,
  100,  75,  0.10,  0.60,
  75, 100,  0.10,  0.20,
  75, 100,  0.10,  0.30,
  75, 100,  0.10,  0.40,
  75, 100,  0.10,  0.50,
  75, 100,  0.10,  0.60,
  50, 125,  0.10,  0.20,
  50, 125,  0.10,  0.30,
  50, 125,  0.10,  0.40,
  50, 125,  0.10,  0.50,
  50, 125,  0.10,  0.60,
  25, 150,  0.10,  0.20,
  25, 150,  0.10,  0.30,
  25, 150,  0.10,  0.40,
  25, 150,  0.10,  0.50,
  25, 150,  0.10,  0.60
) %>%
  mutate(
    n.ratio = n2 / n1,
    diff_rho = rho2 - rho1,
    power = pmap_dbl(
      list(n1, n.ratio, rho1, rho2),
      function(n1, n.ratio, rho1, rho2) {
        corr.2samp(
          n1 = n1,
          n.ratio = n.ratio,
          rho1 = rho1,
          rho2 = rho2,
          alpha = 0.025,
          sides = 2
        )
      }
    ),
    `Sample Size Split` = paste0(n1, "/", n2),
    `Correlation Difference` = round(diff_rho, 2),
    Power = round(power, 3)
  ) %>%
  select(`Correlation Difference`, `Sample Size Split`, Power) %>%
  pivot_wider(
    names_from = `Sample Size Split`,
    values_from = Power
  ) %>%
  arrange(`Correlation Difference`)


#plots for all three power calculation
ggplot(aim1_joint_grid, aes(x = Rsq, y = power_p4)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Aim 1 joint-model power for omnibus test of 4 inflammatory predictors",
    x = expression(R^2),
    y = "Power"
  ) +
  theme_minimal()

ggplot(aim1_corr_grid, aes(x = rhoA, y = power)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Aim 1 separate-model power for single predictor-outcome association",
    x = "Assumed correlation",
    y = "Power"
  ) +
  theme_minimal()


