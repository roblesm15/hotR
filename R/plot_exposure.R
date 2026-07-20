#' Exposure response (temperature mortality) curve
#'
#' Plots the flexible spline-estimated relative risk curve with a vertical
#' line at the minimum mortality temperature (MMT).
#'
#' @param flex_model_fit A \code{"flexible_fit"} object returned by
#'   \code{fit_flexible()}. It must contain \code{$exposure} and
#'   \code{$mmt}.
#'
#' @return Named list: \code{plot} (ggplot) and \code{rr_df} (data frame).
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_vline labs
#'
#' @examples
#' \dontrun{
#' city_data <- puerto_rico_counts_tmax[1:365, ]
#' fit <- fit_flexible(city_data, lag = 2, sim_boot = FALSE)
#' exposure_response(fit)$plot
#' }
exposure_response <- function(flex_model_fit) {
  if (!inherits(flex_model_fit, "flexible_fit")) {
    stop(
      "flex_model_fit must be an object returned by fit_flexible().",
      call. = FALSE
    )
  }
  if (is.null(flex_model_fit$exposure)) {
    stop("flex_model_fit does not contain an exposure-response curve.",
         call. = FALSE)
  }

  rr_df <- as.data.frame(flex_model_fit$exposure)
  required_cols <- c("temp", "RR", "low", "high")
  if (!all(required_cols %in% names(rr_df))) {
    stop("flex_model_fit$exposure has an unexpected structure.",
         call. = FALSE)
  }
  rr_df <- rr_df[required_cols]
  names(rr_df) <- c("temperature", "rr", "rr_lower", "rr_upper")

  mmt_city <- flex_model_fit$mmt
  if (!is.numeric(mmt_city) || length(mmt_city) != 1L ||
      !is.finite(mmt_city)) {
    stop("flex_model_fit$mmt must be a finite numeric scalar.",
         call. = FALSE)
  }

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
