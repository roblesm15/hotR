#' Effect plot: excess mortality rate vs temperature
#'
#' Shows the fitted smooth threshold effect, residual scatter, the
#' estimated heat-onset risk threshold, its Wald confidence interval,
#' and the 97.5th percentile of the lagged temperature exposure.
#'
#' @param data Data frame used to fit \code{model_fit}.
#' @param model_fit A \code{"hot_fit"} object fitted with
#'   \code{return_data_mat = TRUE}.
#' @param scale_pop Numeric. Reference population for rate scaling.
#'   Default is 1,000,000.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom dplyr mutate arrange
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon
#'   geom_vline geom_hline geom_rect geom_label scale_linetype_manual
#'   scale_y_continuous xlim labs
#'
#' @examples
#' \dontrun{
#' city_data <- puerto_rico_counts_tmax[1:365, ]
#' fit <- fit_hot(city_data, L = 1, return_data_mat = TRUE)
#' effect_plot(city_data, fit)
#' }
effect_plot <- function(data,
                        model_fit,
                        scale_pop = 1e6) {
  validate_city_data(data)
  validate_positive_scalar(scale_pop, "scale_pop")
  if (!inherits(model_fit, "hot_fit")) {
    stop("model_fit must be an object returned by fit_hot().",
         call. = FALSE)
  }
  if (is.null(model_fit$X_full) || is.null(model_fit$dates)) {
    stop(
      "Refit with return_data_mat = TRUE before calling effect_plot().",
      call. = FALSE
    )
  }
  if (nrow(data) != nrow(model_fit$X_full) ||
      !identical(data$date, model_fit$dates)) {
    stop("data must be the same ordered daily data used to fit model_fit.",
         call. = FALSE)
  }

  tau_hat <- model_fit$c_hat
  tau_hat_lower <- model_fit$ci_c_wald[1]
  tau_hat_upper <- model_fit$ci_c_wald[2]

  beta_hat <- as.numeric(model_fit$coef["g"])
  X_full <- model_fit$X_full
  p_gamma <- ncol(X_full)
  fit_coef <- model_fit$coef[seq_len(p_gamma)]

  eta_base <- as.vector(X_full %*% fit_coef + log(data$population))
  mu_base <- exp(eta_base)

  mu_ref  <- mean(mu_base)
  pop_ref <- mean(data$population)
  c_ref   <- scale_pop * mu_ref / pop_ref

  # Match the padding and lag definition used during model fitting.
  temp <- extract_single_lag(data$temperature, model_fit$L)
  hinge <- smooth_hinge(temp, c = tau_hat, k = model_fit$k)
  log_rr <- beta_hat * hinge

  # Pointwise delta-method interval using the joint covariance of b and c.
  idx_b <- p_gamma + 1L
  idx_c <- p_gamma + 2L
  V_bc <- model_fit$V_wald[c(idx_b, idx_c), c(idx_b, idx_c), drop = FALSE]
  effect_gradient <- cbind(
    beta_hat * hinge,
    beta_hat * dg_dc(temp, c = tau_hat, k = model_fit$k)
  )
  var_log_rr <- rowSums((effect_gradient %*% V_bc) * effect_gradient)
  se_log_rr <- sqrt(pmax(var_log_rr, 0))

  rr_fit <- exp(log_rr)
  rr_lower <- exp(log_rr - qnorm(0.975) * se_log_rr)
  rr_upper <- exp(log_rr + qnorm(0.975) * se_log_rr)

  city_effect_df <- dplyr::arrange(
    dplyr::mutate(
      data.frame(
        temperature_avg = temp,
        population      = data$population,
        deaths          = data$deaths,
        mu_base         = mu_base,
        rr_fit          = rr_fit,
        rr_lower        = rr_lower,
        rr_upper        = rr_upper
      ),
      resid_rate            = scale_pop * (deaths - mu_base) / population,
      effect_rate_ref_fit   = c_ref * (rr_fit - 1),
      effect_rate_ref_lower = c_ref * (rr_lower - 1),
      effect_rate_ref_upper = c_ref * (rr_upper - 1)
    ),
    temperature_avg
  )

  percentile <- stats::quantile(temp, probs = 0.975, names = FALSE)

  ann_df <- data.frame(
    x    = tau_hat,
    xlow = tau_hat_lower,
    xup  = tau_hat_upper,
    lab  = paste0("hat(tau)==", round(tau_hat, 2), "*degree*C")
  )

  vline_df <- data.frame(
    xintercept = percentile,
    label      = "97.5th Percentile"
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
    ggplot2::geom_line(
      ggplot2::aes(y = effect_rate_ref_fit), color = "steelblue"
    ) +
    ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = xintercept, linetype = label),
      colour = "black"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(
      data = ann_df,
      ggplot2::aes(xintercept = x)
    ) +
    ggplot2::geom_rect(
      data = ann_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = xlow, xmax = xup, ymin = -Inf, ymax = Inf),
      alpha = 0.20
    ) +
    ggplot2::geom_label(
      data = ann_df,
      ggplot2::aes(x = x, y = Inf, label = lab),
      inherit.aes = FALSE, size = 3, parse = TRUE, vjust = 1.5
    ) +
    ggplot2::scale_linetype_manual(
      name   = NULL,
      values = c("97.5th Percentile" = "dotted")
    ) +
    ggplot2::scale_y_continuous(
      name     = "Excess mortality rate (per million)",
      sec.axis = ggplot2::sec_axis(
        ~ 1 + . / c_ref,
        name = "Relative Risk (RR)"
      )
    ) +
    ggplot2::xlim(min(temp), max(temp)) +
    ggplot2::labs(
      title = "Temperature Effect on Excess Mortality Rate",
      x = paste0(
        "Temperature (lag ", model_fit$L, " day",
        if (model_fit$L == 1L) "" else "s", ") (\u00b0C)"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
