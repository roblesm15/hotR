#' Exposure response (temperature mortality) curve
#'
#' Plots the flexible spline-estimated relative risk curve with a vertical
#' line at the minimum mortality temperature (MMT).
#'
#' @param flex_model_fit Flexible model fit object; must
#'                        contain \code{$fit$exposure} (matrix / data frame
#'                        with columns temperature, rr, rr_lower, rr_upper)
#'                        and \code{$fit$mmt}.
#'
#' @return Named list: \code{plot} (ggplot) and \code{rr_df} (data frame).
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_vline labs
exposure_response <- function(flex_model_fit) {

  res     <- flex_model_fit
  rr_df   <- as.data.frame(res$exposure)
  colnames(rr_df) <- c("temperature", "rr", "rr_lower", "rr_upper")


  mmt_city <- res$mmt

  p <- ggplot2::ggplot(rr_df, ggplot2::aes(x = temperature, y = rr)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = rr_lower, ymax = rr_upper),
                         alpha = 0.2, fill = "steelblue") +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_vline(xintercept = mmt_city,
                        linetype = "dashed") +
    ggplot2::labs(
      title = "Temperature mortality curve",
      x     = "Temperature (\u00b0C)",
      y     = "Relative Risk (RR)"
    ) + ggplot2::theme_bw()

  list(plot = p, rr_df = rr_df)
}