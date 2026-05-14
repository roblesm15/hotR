#' Plot the estimated date (time-trend) effect
#'
#' Overlays fitted baseline deaths (spline date trend) with pointwise 95 \%
#' confidence ribbons on the raw daily death counts.
#' @param data        Data frame containing at least \code{date},
#'                    \code{deaths}, and \code{population} columns.
#' @param model_fit  Model fit object.  Each
#'                    element must contain \code{$fit$X_spline},
#'                    \code{$fit$coef}, and \code{$fit$V_wald}.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon labs
date_effect_plot <- function( data, model_fit) {

  city_data <- data
  n_coef    <- length(model_fit$coef)

  date_spline <- cbind(1, model_fit$X_spline)
  date_effect <- date_spline %*% model_fit$coef[-n_coef] + log(city_data$population)

  vcov_sub <- model_fit$V_wald[seq_len(n_coef - 1), seq_len(n_coef - 1)]
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