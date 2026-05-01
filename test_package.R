devtools::load_all()

library(tidyverse)
# Load datasets
load(file = "../india_heat/data/india_mort_temp_ah.rda")
load(file = "../india_heat/data/exclusion_dates_india.rda")
load(file = "../india_heat/data/full_exclusion_dates_india.rda")

# Prepare India-specific dataset with exclusions for known PEMs
data_india <- india_mort_weather |>
  left_join(exclusion_dates_india, by = c("date", "city")) |>
  mutate(excl_period = ifelse(is.na(excl_period), "no", excl_period)) |>
  filter(excl_period == "no") |>
  mutate(day_of_year = yday(date)) |>
  dplyr::select(-excl_period) |>
  dplyr::rename(population = pop)

city_data <- data_india |>
    filter(city == "Delhi")  |>
    dplyr::select(date, max_temp, deaths, population) |>
    dplyr::rename(temperature = max_temp)
colnames(city_data)


fit <- fit_hot(
  city_data          = city_data,
  df_per_year_spline = 3,
  L                  = 1, # single previous day lag
  verbose            = TRUE
)

str(fit)

load(file = "../heat/usa/rdas/us_cities_mort_temp.rda")

colnames(us_cities_mort_temp)

us_city_data <- us_cities_mort_temp |>
  filter(city == "Chicago") |>
  dplyr::select(date, TEMP, tot, population) |>
  dplyr::rename(temperature = TEMP, deaths = tot) |>
  # make temperature in celsius
  mutate(temperature = (temperature - 32) * 5/9) |>
  group_by(date) |>
  summarize(temperature = mean(temperature, na.rm = TRUE),
            deaths = sum(deaths, na.rm = TRUE),
            population = sum(population, na.rm = TRUE))


fit_us <- fit_hot(
  city_data          = us_city_data,
  df_per_year_spline = 4,
  L                  = 1, # single day lag
  lower_q            = 0.01,
  verbose            = TRUE
)


fit_flexible_us <- fit_flexible(city_data = us_city_data,
                         df_per_year_spline = 4, lag = 7, 
                         return_data_mat = FALSE, sim_boot = FALSE)

fit_flexible_us$mmt

# fit simple long-term trend model
library(splines)

fit_long_term <- glm.nb(deaths ~ ns(date, df = 4 * length(unique(year(date)))) + offset(log(population)),
          data = us_city_data)

# residuals

residuals_long_term <- us_city_data$deaths- predict(fit_long_term, type = "response")

plot(us_city_data$temperature, residuals_long_term, main = "Residuals of Long-Term Trend Model", xlab = "Date", ylab = "Residuals")

plot(us_city_data$date, us_city_data$deaths, main = "Residuals of Long-Term Trend Model", xlab = "Date", ylab = "Residuals")
# add 
lines(us_city_data$date, predict(fit_long_term, type = "response"), col = "red")

