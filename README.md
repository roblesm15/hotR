
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# hotR

`hotR` estimates the heat-onset risk threshold (HOT): the temperature at
which mortality begins to increase above its seasonal background. It
fits a negative binomial structured threshold model with a smooth heat
effect and a natural-spline time trend. The package also includes a
flexible distributed lag nonlinear model (DLNM) for comparison.

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github(
  "roblesm15/hotR",
  build_vignettes = TRUE
)
```

## Example

Input data must have columns named `date`, `temperature`, `deaths`, and
`population`.

``` r
library(hotR)

pr_counts <- puerto_rico_counts_tmax

# Estimate the HOT using previous-day temperature. Retaining the model
# matrices enables the plotting helpers.
fit <- fit_hot(
  city_data = pr_counts,
  df_per_year_spline = 2,
  L = 1,
  return_data_mat = TRUE
)

fit$c_hat             # estimated heat-onset risk threshold (°C)
fit$ci_c_wald         # model-based 95% confidence interval
fit$ci_c_sandwich     # observation-level sandwich 95% interval
fit$coef["g"]         # heat-effect coefficient

date_effect_plot(pr_counts, fit)
effect_plot(pr_counts, fit)

# Fit a less structured DLNM comparison. Set sim_boot = TRUE to simulate
# an uncertainty interval for the minimum mortality temperature (MMT).
flex_fit <- fit_flexible(
  city_data = pr_counts,
  df_per_year_spline = 2,
  lag = 5,
  sim_boot = FALSE
)

flex_fit$mmt
exposure_response(flex_fit)$plot
```

The full Puerto Rico series spans 1985–2022, so the two fits can take
several minutes. For a quick code check, use a contiguous subset such as
`puerto_rico_counts_tmax[1:365, ]`; estimates from such a short subset
should not be treated as a substantive mortality analysis.

After installing with vignette building enabled, open the detailed
walkthrough with:

``` r
vignette("using-hotR", package = "hotR")
```
