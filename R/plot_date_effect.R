#' Plot the estimated date (time-trend) effect
#'
#' Overlays fitted baseline deaths (spline date trend) with pointwise 95%
#' confidence ribbons on the raw daily death counts.
#' @param data Data frame used to fit \code{model_fit}.
#' @param model_fit A \code{"hot_fit"} object fitted with
#'   \code{return_data_mat = TRUE}.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon labs
#'
#' @examples
#' \dontrun{
#' city_data <- puerto_rico_counts_tmax[1:365, ]
#' fit <- fit_hot(city_data, L = 1, return_data_mat = TRUE)
#' date_effect_plot(city_data, fit)
#' }
date_effect_plot <- function(data, model_fit) {
  validate_city_data(data)
  if (!inherits(model_fit, "hot_fit")) {
    stop("model_fit must be an object returned by fit_hot().",
         call. = FALSE)
  }
  if (is.null(model_fit$X_full) || is.null(model_fit$dates)) {
    stop(
      "Refit with return_data_mat = TRUE before calling date_effect_plot().",
      call. = FALSE
    )
  }
  if (nrow(data) != nrow(model_fit$X_full) ||
      !identical(data$date, model_fit$dates)) {
    stop("data must be the same ordered daily data used to fit model_fit.",
         call. = FALSE)
  }

  city_data <- data
  p_gamma <- ncol(model_fit$X_full)

  date_spline <- model_fit$X_full
  date_effect <- date_spline %*% model_fit$coef[seq_len(p_gamma)] +
    log(city_data$population)

  vcov_sub <- model_fit$V_wald[seq_len(p_gamma), seq_len(p_gamma), drop = FALSE]
  se_date  <- sqrt(rowSums((date_spline %*% vcov_sub) * date_spline))

  plot_df <- dplyr::mutate(
    city_data,
    fitted = exp(date_effect),
    lo     = exp(date_effect - 1.96 * se_date),
    hi     = exp(date_effect + 1.96 * se_date)
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = date)) +
    ggplot2::geom_point(ggplot2::aes(y = deaths),
                        alpha = 0.25, colour = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "steelblue", alpha = 0.20) +
    ggplot2::geom_line(ggplot2::aes(y = fitted), colour = "steelblue") +
    ggplot2::labs(title = "Date Effect", x = "Date", y = "Deaths") + 
    ggplot2::theme_bw()
}
