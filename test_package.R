# Run this integration example from the hotR package root.
if (!file.exists("DESCRIPTION")) {
  stop("Run test_package.R from the hotR package root.")
}

devtools::load_all(quiet = TRUE)
pr_counts <- puerto_rico_counts_tmax

pr_counts$city <- "Puerto Rico"
pr_counts$year <- as.integer(format(pr_counts$date, "%Y"))

fit_long_term <- MASS::glm.nb(
  deaths ~ splines::ns(as.numeric(date),
                       df = 2 * length(unique(year))) +
    offset(log(population)),
  data = pr_counts
)

# Residual checks for the seasonal background model.
residuals_long_term <- pr_counts$deaths -
  predict(fit_long_term, type = "response")

plot(
  pr_counts$temperature,
  residuals_long_term,
  main = "Residuals of Long-Term Trend Model",
  xlab = "Temperature (°C)",
  ylab = "Residuals"
)

ggplot2::ggplot(
  pr_counts,
  ggplot2::aes(x = temperature, y = residuals_long_term)
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(
    method = "loess", se = FALSE, color = "red", span = 0.25
  ) +
  ggplot2::labs(
    title = "Residuals of Long-Term Trend Model",
    x = "Temperature (°C)",
    y = "Residuals"
  )

plot(
  pr_counts$date,
  pr_counts$deaths,
  main = "Long-Term Trend Model",
  xlab = "Date",
  ylab = "Deaths"
)
lines(pr_counts$date, predict(fit_long_term, type = "response"), col = "red")

# The full-series model fit may take several minutes.
fit <- fit_hot(
  city_data          = pr_counts,
  df_per_year_spline = 2,
  L                  = 1, # previous-day temperature
  verbose            = TRUE, 
  return_data_mat    = TRUE # needed by the plotting functions
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
                                df_per_year_spline = 2,
                                lag = 5,
                                return_data_mat = TRUE,
                                sim_boot = FALSE,
                                n_boot = 100)
fit_flexible_pr$mmt

hotR::exposure_response(
   fit_flexible_pr
)
