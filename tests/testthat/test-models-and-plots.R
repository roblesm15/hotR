test_that("fit_hot returns coherent estimates and plotting inputs", {
  city_data <- puerto_rico_counts_tmax[1:365, ]

  fit <- suppressWarnings(fit_hot(
    city_data = city_data,
    df_per_year_spline = 2,
    L = 1,
    maxit = 30,
    return_data_mat = TRUE
  ))

  expect_s3_class(fit, "hot_fit")
  expect_true(fit$converged)
  expect_identical(fit$optimizer_convergence, 0L)
  expect_true(is.finite(fit$c_hat))
  expect_true(is.finite(fit$theta_hat) && fit$theta_hat > 0)
  expect_true(is.finite(fit$coef["g"]) && fit$coef["g"] > 0)
  expect_named(fit$ci_c_wald, c("lower", "upper"))
  expect_named(fit$ci_c_sandwich, c("lower", "upper"))
  expect_equal(nrow(fit$X_full), nrow(city_data))
  expect_identical(fit$dates, city_data$date)

  date_plot <- date_effect_plot(city_data, fit)
  temperature_plot <- effect_plot(city_data, fit)
  expect_s3_class(date_plot, "ggplot")
  expect_s3_class(temperature_plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(date_plot))
  expect_silent(suppressWarnings(ggplot2::ggplot_build(temperature_plot)))
})

test_that("fit_flexible returns curves that can be plotted", {
  city_data <- puerto_rico_counts_tmax[1:365, ]

  fit <- suppressWarnings(fit_flexible(
    city_data = city_data,
    df_per_year_spline = 2,
    lag = 2,
    return_data_mat = FALSE,
    sim_boot = FALSE
  ))

  expect_s3_class(fit, "flexible_fit")
  expect_true(is.finite(fit$mmt))
  expect_true(is.data.frame(fit$exposure))
  expect_true(is.data.frame(fit$lag_effects))
  expect_true(is.data.frame(fit$lag_importance))
  expect_null(fit$cb_temp)
  expect_null(fit$X_spline)
  expect_true(all(is.na(fit$ci)))

  response <- exposure_response(fit)
  expect_named(response, c("plot", "rr_df"))
  expect_s3_class(response$plot, "ggplot")
  expect_named(
    response$rr_df,
    c("temperature", "rr", "rr_lower", "rr_upper")
  )
  expect_silent(ggplot2::ggplot_build(response$plot))
})

test_that("model inputs fail early with informative errors", {
  city_data <- puerto_rico_counts_tmax[1:30, ]

  missing_temperature <- city_data
  missing_temperature$temperature[1] <- NA_real_
  expect_error(
    fit_hot(missing_temperature),
    "missing values in required columns"
  )

  unordered <- city_data[c(2, 1, 3:nrow(city_data)), ]
  expect_error(fit_hot(unordered), "strictly increasing")

  gap <- city_data[-10, ]
  expect_error(fit_flexible(gap, sim_boot = FALSE), "no gaps")

  expect_error(fit_hot(city_data, L = nrow(city_data)),
               "smaller than the number")
  expect_error(fit_hot(city_data, lower_q = 0.9, upper_q = 0.1),
               "lower_q must be smaller")
  expect_error(fit_flexible(city_data, lag = -1, sim_boot = FALSE),
               "greater than or equal to 0")
})
