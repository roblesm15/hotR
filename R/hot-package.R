#' hotR: Heat-Mortality Analysis Tools
#'
#' Provides simulation, modelling, and visualisation utilities for
#' analysing the relationship between extreme heat and mortality.
#'
#' @keywords internal
"_PACKAGE"

## Suppress R CMD CHECK notes for NSE variables used across the package
utils::globalVariables(c(
  # data columns
  "date", "deaths", "temperature", "population", "city",
  "day", "ah",
  # computed columns
  "rr", "rr_lower", "rr_upper", "rr_flexible",
  "rr_flexible_lower", "rr_flexible_upper",
  "effect_rate_ref_fit", "effect_rate_ref_lower", "effect_rate_ref_upper",
  "effect_rate_ref_flex", "effect_rate_ref_flex_low", "effect_rate_ref_flex_up",
  "resid_rate", "temperature_avg",
  "weight", "se_weights", "lag",
  "effect", "effect_lower", "effect_upper",
  "pearson_residuals", "residuals",
  # annotation helpers
  "x", "xlow", "xup", "lab",
  "xintercept", "label",
  # lag helpers
  "RR", "low", "high", "temp", "se_effect"
))