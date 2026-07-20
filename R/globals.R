# R/globals.R

# Suppress R CMD check notes for variables used in non-standard evaluation.
utils::globalVariables(c(
  "date", "deaths", "fitted", "hi", "lo",
  "temperature", "temperature_avg",
  "population", "resid_rate",
  "effect_rate_ref_fit", "effect_rate_ref_lower", "effect_rate_ref_upper",
  "rr", "rr_lower", "rr_upper",
  "x", "xlow", "xup", "lab", "xintercept", "label"
))
