
# 1. Create the package structure locally
usethis::create_package("hotR")

# 2. Initialize git inside the package
usethis::use_git()

# 3. Push to GitHub and create the repo simultaneously
usethis::use_github()

# 4. Add a description file

usethis::use_description(fields = list(
  Title = "Structured Threshold Model for Heat-Onset Risk Threshold Estimation",
  Description = "Implements the structured threshold model (STM) for 
    estimating the heat-onset risk threshold (HOT) in mortality time series.",
  `Authors@R` = 'person("Monica", "Robles Fontan", email = "mroblesfontan@g.harvard.edu", 
    role = c("aut", "cre"))'
))

usethis::use_mit_license()

usethis::use_readme_rmd()

# ── Inside your package project (hotR/) ───────────────────────────────────────
# Run this ONCE to create the .rda file

library(usethis)  # install if needed

# Assume you already have your cleaned df from the NASA POWER call
# (your San Juan max temp data frame)

puerto_rico_counts_tmax <- pr_counts 

# ── Save to data/ folder ──────────────────────────────────────────────────────
usethis::use_data(
  puerto_rico_counts_tmax,       # object name = what users will call with data()
  overwrite  = TRUE,  # overwrite if rebuilding
  compress   = "xz"   # best compression for CRAN; "gzip" also fine
)

# ✅ Creates:  hotR/data/puerto_rico_counts_tmax.rda


usethis::use_package("dplyr")
usethis::use_package("ggplot2")
usethis::use_package("MASS")

# ── Run these in order after editing the files ────────────────────────────────

devtools::document()    # generate NAMESPACE + man/*.Rd from roxygen2 tags
devtools::check()       # R CMD CHECK — look for NOTEs / WARNINGs


