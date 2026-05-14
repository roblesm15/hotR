#' Effect plot: excess mortality rate vs temperature
#'
#' Shows the fitted quadratic threshold effect, residual scatter, and key reference line for 97.5th percentile, and \eqn{\hat{\tau}}).
#'
#' @param data            Data frame with city-level daily data.
#' @param model_fit       Threshold model fit.
#' @param scale_pop       Numeric. Reference population for rate scaling
#'                        (default 1 000 000).
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom dplyr filter pull mutate arrange left_join
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon
#'   geom_vline geom_hline geom_rect geom_label scale_linetype_manual
#'   scale_y_continuous xlim labs
effect_plot <- function(data,
                        model_fit,
                        scale_pop = 1e6) {


  tau_hat <- model_fit$c_hat
  tau_hat_lower <- model_fit$ci_c_wald[1]
  tau_hat_upper <- model_fit$ci_c_wald[2]

  beta_hat <- model_fit$coef['g']
  beta_lower <- beta_hat - qnorm(0.975) * model_fit$se_wald['se_beta']
  beta_upper <- beta_hat + qnorm(0.975) * model_fit$se_wald['se_beta']
  
  fit_obj    <- model_fit
  X_spline   <- cbind(1, fit_obj$X_spline)
  fit_coef   <- fit_obj$coef[-length(fit_obj$coef)]

  eta_base   <- as.vector(X_spline %*% fit_coef + log(data$population))
  mu_base    <- exp(eta_base)

  mu_ref     <- mean(mu_base,            na.rm = TRUE)
  pop_ref    <- mean(data$population, na.rm = TRUE)
  c_ref      <- scale_pop * mu_ref / pop_ref        # secondary-axis constant


  lag_weights <- fit_obj$lag_weights
  temp_lags   <- sapply(seq(0, length(lag_weights) - 1), function(lag) {
    dplyr::lag(data$temperature, n = lag)
  })

  temp <- as.vector(temp_lags %*% lag_weights)
  
  temp_effect <- function(temp, beta, tau) {
    beta * smooth_hinge(temp, c = tau, k = 5)
  }

  rr_fit   <- exp(temp_effect(temp, beta_hat,   tau_hat))
  rr_lower <- exp(temp_effect(temp, beta_lower, tau_hat))
  rr_upper <- exp(temp_effect(temp, beta_upper, tau_hat))

  
  city_effect_df <- dplyr::arrange(
    dplyr::mutate(
      data.frame(
        temperature_avg = temp,
        temperature     = data$temperature,
        population      = data$population,
        deaths          = data$deaths,
        mu_base         = mu_base,
        rr_fit          = rr_fit,
        rr_lower        = rr_lower,
        rr_upper        = rr_upper
      )
    ),
    temperature_avg
  ) |>
    dplyr::mutate(
      resid_rate               = scale_pop * (deaths - mu_base) / population,
      effect_rate_ref_fit      = c_ref * (rr_fit - 1),
      effect_rate_ref_lower    = c_ref * (rr_lower - 1),
      effect_rate_ref_upper    = c_ref * (rr_upper - 1)
    )


  percentile <- stats::quantile(data$temperature, probs = 0.975,
                                na.rm = TRUE)

  ann_df <- data.frame(
    x    = tau_hat,
    xlow = tau_hat_lower,
    xup  = tau_hat_upper,
    lab  = paste0("hat(tau)=='", round(tau_hat, 2), "\u00b0C'")
  )

  vline_df <- data.frame(
    xintercept = c(percentile),
    label      = c("97.5th Percentile")
  )

  
  ggplot2::ggplot(city_effect_df,
                  ggplot2::aes(x = temperature_avg)) +
    ggplot2::geom_point(ggplot2::aes(y = resid_rate),
                        colour = "grey", size = 0.5, alpha = 0.25) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = effect_rate_ref_lower,
                   ymax = effect_rate_ref_upper),
      alpha = 0.25, fill = "steelblue"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = effect_rate_ref_fit), color = "steelblue") +
    ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = xintercept, linetype = label),
      colour = "black"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(data = ann_df,
                        ggplot2::aes(xintercept = x)) +
    ggplot2::geom_rect(
      data = ann_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = xlow, xmax = xup, ymin = -Inf, ymax = Inf), alpha = 0.20
    ) +
    ggplot2::geom_label(
      data = ann_df,
      ggplot2::aes(x = x, y = -1.05 * 10, label = lab), size = 3, parse = TRUE
    ) +
    ggplot2::scale_linetype_manual(
      name   = NULL,
      values = c("97.5th Percentile" = "dotted", "IMD" = "dotdash")
    ) +
    ggplot2::scale_y_continuous(
      name     = "Excess mortality rate (per million)",
      sec.axis = ggplot2::sec_axis(~ 1 + . / c_ref,
                                   name = "Relative Risk (RR)")
    ) +
    ggplot2::xlim(min(data$temperature, na.rm = TRUE), max(data$temperature, na.rm = TRUE)) +
    ggplot2::labs(
      title = "Temperature Effect on Excess Mortality Rate",
      x     = "Maximum Temperature (lag 1) (\u00b0C)"
    )+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

