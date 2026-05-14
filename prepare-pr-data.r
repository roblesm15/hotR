# First get daily deaths and population estimates from excessmort

library(tidyverse)
library(excessmort)

# Load datasets

pr_counts <- puerto_rico_counts |> 
    group_by(date) |>
    summarize(deaths = sum(outcome, na.rm = TRUE),
              population = sum(population, na.rm = TRUE))

# 2023-12-31
# Now lets call temperature data from nasapower

# ── Packages ──────────────────────────────────────────────────────────────────
if (!require("httr"))     install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("dplyr"))    install.packages("dplyr")
if (!require("ggplot2"))  install.packages("ggplot2")

library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)

# ── Location: San Juan, Puerto Rico ───────────────────────────────────────────
lat <- 18.4655
lon <- -66.1057

# ── API Parameters ────────────────────────────────────────────────────────────
# T2M_MAX  = Maximum temperature at 2 meters above surface (°C)
# community options: RE (Renewable Energy) | AG (Agroclimatology) | SB (Buildings)

base_url <- "https://power.larc.nasa.gov/api/temporal/daily/point"

query_params <- list(
  parameters = "T2M_MAX",
  community  = "RE",
  latitude   = lat,
  longitude  = lon,
  start      = "19850101",   # YYYYMMDD
  end        = "20221231",
  format     = "JSON"
)

# ── API Call ──────────────────────────────────────────────────────────────────
cat("📡 Calling NASA POWER API...\n")

response <- GET(
  url   = base_url,
  query = query_params,
  timeout(120)
)

# ── Status Check ──────────────────────────────────────────────────────────────
cat("HTTP Status:", status_code(response), "\n")

if (http_error(response)) {
  stop("❌ Request failed: ", status_code(response),
       "\n", content(response, as = "text"))
}

# ── Parse JSON ────────────────────────────────────────────────────────────────
raw_text  <- content(response, as = "text", encoding = "UTF-8")
data_json <- fromJSON(raw_text, flatten = TRUE)

# Metadata peek
cat("Parameter unit:",
    data_json$parameters$T2M_MAX$units, "\n")   # should be "C"
cat("Location confirmed — lat:",
    data_json$geometry$coordinates[2],
    " lon:", data_json$geometry$coordinates[1], "\n")

# ── Build Data Frame ──────────────────────────────────────────────────────────
tmax_list <- data_json$properties$parameter$T2M_MAX

df <- data.frame(
  date    = names(tmax_list),
  T2M_MAX = as.numeric(unlist(tmax_list)),
  stringsAsFactors = FALSE
) |>
  mutate(
    date    = as.Date(date, format = "%Y%m%d"),
    temperature = if_else(T2M_MAX == -999, NA_real_, T2M_MAX),  # fill value → NA
    month   = format(date, "%b")
  )


# complete dataset

pr_counts <- pr_counts |> 
    left_join(df, by = "date") |> 
    arrange(date) |> 
    select(date, deaths, population, temperature)

