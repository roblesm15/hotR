setwd("hotR")

devtools::load_all()
library(MASS)
library(splines)
library(tidyverse)
pr_counts <- puerto_rico_counts_tmax

pr_counts <- pr_counts |>
  mutate(city = "Puerto Rico")

fit_long_term <- glm.nb(deaths ~ ns(date, df = 2 * length(unique(year(date)))) + offset(log(population)),
          data = pr_counts)

# residuals

residuals_long_term <- pr_counts$deaths- predict(fit_long_term, type = "response")

plot(pr_counts$temperature, residuals_long_term, main = "Residuals of Long-Term Trend Model", xlab = "Date", ylab = "Residuals")

ggplot(pr_counts, aes(x = temperature, y = residuals_long_term)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.25) +
  labs(title = "Residuals of Long-Term Trend Model", x = "Date", y = "Residuals")

plot(pr_counts$date, pr_counts$deaths, main = " Long-Term Trend Model", xlab = "Date", ylab = "Residuals")
# add 
lines(pr_counts$date, predict(fit_long_term, type = "response"), col = "red")


fit <- fit_hot(
  city_data          = pr_counts,
  df_per_year_spline = 2,
  L                  = 1, # single previous day lag
  verbose            = TRUE, 
  return_data_mat            = TRUE # for plots info
)

fit$c_hat
fit$ci_c_wald
fit$coef['g']
fit$se_wald['se_beta']

fit$coef['g'] + c(-1, 1) * qnorm(0.975) * fit$se_wald['se_beta']


hotR::date_effect_plot(
  data = pr_counts,
  model_fit = fit
)

hotR::effect_plot(
  data = pr_counts,
  model_fit = fit
)


fit_flexible_pr <- fit_flexible(city_data = pr_counts,
                         df_per_year_spline = 2, lag = 5, 
                         return_data_mat = TRUE, sim_boot = !TRUE, n_boot = 100)

hotR::exposure_response(
   fit_flexible_pr
)
