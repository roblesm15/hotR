# flexible.R
# User-facing interface for the flexible distributed lag model.
# -------------------------------------------------------------------
# Minimum mortality temperature estimation
# -------------------------------------------------------------------
#' Find the minimum mortality temperature from a DLNM fit
#'
#' @noRd
findmin <- function(basis,
                    model = NULL,
                    coef = NULL,
                    vcov = NULL,
                    at = NULL,
                    from = NULL,
                    to = NULL,
                    by = NULL,
                    sim = FALSE,
                    nsim = 5000,
                    basis_name = "cb_temp") {

  if (!any(class(basis) %in% c("crossbasis", "onebasis"))) {
    stop(
      "the first argument must be an object of class ",
      "'crossbasis' or 'onebasis'"
    )
  }

  if (is.null(at)) {
    if (is.null(from)) {
      from <- attr(basis, "range")[1]
    }
    if (is.null(to)) {
      to <- attr(basis, "range")[2]
    }
    if (is.null(by)) {
      by <- 0.1
    }

    at <- seq(from, to, by = by)
  }

  at <- sort(unique(at[is.finite(at)]))

  if (length(at) < 2L) {
    stop("'at' must contain at least two finite values.")
  }

  if (!is.null(model)) {
    coef_all <- stats::coef(model)
    vcov_all <- stats::vcov(model)

    cb_ind <- grep(paste0("^", basis_name), names(coef_all))

    if (length(cb_ind) != ncol(basis)) {
      stop("Could not identify cross-basis coefficients in model.")
    }

    coef <- coef_all[cb_ind]
    vcov <- vcov_all[cb_ind, cb_ind, drop = FALSE]
  } else {
    if (is.null(coef) || is.null(vcov)) {
      stop("At least 'model' or 'coef'-'vcov' must be provided.")
    }
  }

  if (length(coef) != ncol(basis) ||
      length(coef) != nrow(vcov) ||
      length(coef) != ncol(vcov) ||
      anyNA(coef) ||
      anyNA(vcov)) {
    stop("model or coef/vcov not consistent with basis.")
  }

  cp <- dlnm::crosspred(
    basis,
    coef = coef,
    vcov = vcov,
    model.link = "log",
    at = at,
    cumul = TRUE,
    cen = at[1]
  )

  pred <- cp$allfit
  min <- at[which.min(pred)]

  if (!sim) {
    return(min)
  }

  eig <- eigen(vcov, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0

  z <- matrix(
    stats::rnorm(length(coef) * nsim),
    nrow = length(coef),
    ncol = nsim
  )

  coefsim <- matrix(coef, nrow = length(coef), ncol = nsim) +
    eig$vectors %*%
    diag(sqrt(eig$values), nrow = length(eig$values)) %*%
    z

  minsim <- numeric(nsim)

  for (i in seq_len(nsim)) {
    cp_i <- dlnm::crosspred(
      basis,
      coef = coefsim[, i],
      vcov = vcov,
      model.link = "log",
      at = at,
      cumul = TRUE,
      cen = at[1]
    )

    minsim[i] <- at[which.min(cp_i$allfit)]
  }

  minsim
}

# -------------------------------------------------------------------
# Main exported function
# -------------------------------------------------------------------

#' Fit a flexible distributed lag model for temperature-mortality analysis
#'
#' Estimates the full temperature-mortality exposure-response
#' relationship using a distributed lag nonlinear model (DLNM).
#' Returns cumulative and lag-specific risk curves relative to the
#' minimum mortality temperature (MMT), along with a summary of lag
#' importance across the temperature distribution.
#'
#' @details
#' Mortality counts \eqn{y_t} are modeled as:
#' \deqn{y_t ~ NegBin(mu_t, r)}
#' \deqn{log(mu_t) = x_t^T gamma + f(T_t, l) + log(pop_t)}
#'
#' where \eqn{f(T_t, l)} is a cross-basis function constructed by
#' \code{dlnm::crossbasis} using natural splines in the temperature
#' direction and integer lag weights up to \code{lag} days. By default,
#' interior temperature knots are placed at the 20th, 40th, 60th, and
#' 80th percentiles of the observed temperature distribution.
#'
#' The MMT is estimated as the temperature minimizing the cumulative
#' exposure-response function. A 95 percent confidence interval is
#' obtained by simulation from the multivariate normal approximation
#' to the distribution of the cross-basis coefficients when
#' \code{sim_boot = TRUE}.
#'
#' @param city_data Data frame with columns \code{date},
#'   \code{temperature}, \code{deaths}, and \code{population}, ordered
#'   by date with one row per day and no gaps.
#' @param df_per_year_spline Positive numeric scalar. Degrees of
#'   freedom per year for the natural spline baseline. Default is
#'   \code{3}.
#' @param ktemp Optional numeric vector. Interior knot positions for
#'   the natural spline in the temperature direction of the cross-basis.
#'   If \code{NULL}, knots are placed at the 20th, 40th, 60th, and
#'   80th percentiles of observed temperature.
#' @param lag Non-negative integer. Maximum lag in days. Default is
#'   \code{7}.
#' @param at Optional numeric vector. Temperature values at which the
#'   exposure-response curve is evaluated. If \code{NULL}, all unique
#'   observed finite temperatures are used.
#' @param return_data_mat Logical. If \code{TRUE}, returns the
#'   cross-basis object and spline design matrix. Exposure-response
#'   summaries are returned regardless of this setting. Default is
#'   \code{TRUE}.
#' @param sim_boot Logical. If \code{TRUE}, simulates uncertainty in
#'   the MMT. Default is \code{TRUE}.
#' @param n_boot Positive integer. Number of simulations for the MMT
#'   confidence interval when \code{sim_boot = TRUE}. Default is
#'   \code{1000}.
#'
#' @return An object of class \code{"flexible_fit"}, a named list with:
#' \describe{
#'   \item{df_loc}{Degrees of freedom used for the time spline.}
#'   \item{lag}{Maximum temperature lag in days.}
#'   \item{mmt}{Estimated minimum mortality temperature.}
#'   \item{ci}{Simulation-based 95% interval for the MMT, or two
#'     \code{NA} values when \code{sim_boot = FALSE}.}
#'   \item{coef}{Fitted negative binomial regression coefficients.}
#'   \item{r}{Estimated negative binomial dispersion parameter.}
#'   \item{cb_temp, X_spline}{Cross-basis and time-spline matrices when
#'     \code{return_data_mat = TRUE}; otherwise \code{NULL}.}
#'   \item{exposure}{Cumulative exposure-response curve.}
#'   \item{lag_effects}{Lag-specific exposure-response curves.}
#'   \item{lag_importance}{Summaries of absolute and squared effects by lag.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use one year here for a quick example; use the full series for analysis.
#' city_data <- puerto_rico_counts_tmax[1:365, ]
#'
#' fit <- fit_flexible(
#'   city_data = city_data,
#'   df_per_year_spline = 2,
#'   lag = 2,
#'   sim_boot = FALSE
#' )
#'
#' fit$mmt
#' exposure_response(fit)$plot
#' }
fit_flexible <- function(
    city_data,
    df_per_year_spline = 3,
    ktemp              = NULL,
    lag                = 7,
    at                 = NULL,
    return_data_mat    = TRUE,
    sim_boot           = TRUE, 
    n_boot             = 1000
) {
  validate_city_data(city_data)
  validate_positive_scalar(df_per_year_spline, "df_per_year_spline")
  validate_integer_scalar(lag, "lag")
  validate_logical_scalar(return_data_mat, "return_data_mat")
  validate_logical_scalar(sim_boot, "sim_boot")

  if (lag >= nrow(city_data)) {
    stop("lag must be smaller than the number of observations.",
         call. = FALSE)
  }
  if (sim_boot) {
    validate_integer_scalar(n_boot, "n_boot", minimum = 1L)
  }

  tnum <- as.numeric(city_data$date)
  temp <- city_data$temperature
  y    <- city_data$deaths
  off  <- log(city_data$population)

  no_years <- as.numeric(max(city_data$date) - min(city_data$date) + 1) / 365.25
  df_loc <- max(4L, round(df_per_year_spline * no_years))

  if (df_loc + 2L >= nrow(city_data) - lag) {
    stop(
      "The requested spline and lag have too many parameters for city_data.",
      call. = FALSE
    )
  }

  if (is.null(ktemp)) {
    ktemp <- stats::quantile(
      temp,
      probs = c(0.2, 0.4, 0.6, 0.8),
      na.rm = TRUE,
      names = FALSE
    )
  }

  if (!is.numeric(ktemp) || length(ktemp) == 0L ||
      any(!is.finite(ktemp))) {
    stop("ktemp must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (is.unsorted(ktemp, strictly = TRUE) || anyDuplicated(ktemp)) {
    stop("ktemp must contain unique values in increasing order.",
         call. = FALSE)
  }
  temp_range <- range(temp)
  if (any(ktemp <= temp_range[1] | ktemp >= temp_range[2])) {
    stop("ktemp values must lie strictly inside the temperature range.",
         call. = FALSE)
  }

  cb_temp <- dlnm::crossbasis(
    x      = temp,
    lag    = lag,
    argvar = list(fun = "ns", knots = ktemp),
    arglag = list(fun = "integer")
  )

  X_spline <- stats::model.matrix(~ splines::ns(tnum, df = df_loc))
  X_spline <- X_spline[, -1, drop = FALSE]

  fit <- MASS::glm.nb(
    y ~ X_spline + cb_temp + offset(off),
    link = log
  )

  if (is.null(at)) {
    at <- sort(unique(temp))
  } else {
    if (!is.numeric(at)) {
      stop("at must be NULL or a numeric vector.", call. = FALSE)
    }
    at <- sort(unique(at[is.finite(at)]))
  }

  if (length(at) < 2L) {
    stop("'at' must contain at least two finite temperature values.")
  }

  mmt <- findmin(
    basis = cb_temp,
    model = fit,
    at    = at,
    sim   = FALSE
  )

  ci <- if (sim_boot) {
    boots <- findmin(
      basis = cb_temp,
      model = fit,
      at    = at,
      sim   = TRUE,
      nsim  = n_boot
    )

    stats::quantile(
      boots,
      probs = c(0.025, 0.975),
      names = TRUE,
      na.rm = TRUE
    )
  } else {
    stats::setNames(c(NA_real_, NA_real_), c("2.5%", "97.5%"))
  }

  cp <- dlnm::crosspred(
    cb_temp,
    fit,
    at     = at,
    bylag  = 1,
    cumul  = TRUE,
    cen    = mmt
  )

  cum_curve <- data.frame(
    temp = as.numeric(names(cp$allRRfit)),
    RR   = as.numeric(cp$allRRfit),
    low  = as.numeric(cp$allRRlow),
    high = as.numeric(cp$allRRhigh),
    row.names = NULL
  )

  lag_curve <- local({
    rr <- cp$matRRfit
    low <- cp$matRRlow
    hi <- cp$matRRhigh

    temps <- as.numeric(rownames(rr))
    lags <- as.integer(gsub("[^0-9]", "", colnames(rr)))

    data.frame(
      temp = rep(temps, times = ncol(rr)),
      lag  = rep(lags, each = nrow(rr)),
      RR   = as.vector(rr),
      low  = as.vector(low),
      high = as.vector(hi),
      row.names = NULL
    )
  })

  lag_importance <- local({
    eta <- cp$matfit

    temps <- as.numeric(rownames(eta))
    lags <- as.integer(gsub("[^0-9]", "", colnames(eta)))

    w <- if (length(temps) > 1L) {
      c(diff(temps), tail(diff(temps), 1))
    } else {
      1
    }

    data.frame(
      lag = lags,
      importance_abs = as.numeric(apply(abs(eta), 2, function(z) sum(z * w))),
      importance_sq  = as.numeric(apply(eta^2, 2, function(z) sum(z * w))),
      mean_abs_eta   = as.numeric(colMeans(abs(eta))),
      max_abs_eta    = as.numeric(apply(abs(eta), 2, max)),
      row.names = NULL
    )
  })

  out <- list(
    df_loc         = df_loc,
    lag            = lag,
    mmt            = mmt,
    ci             = ci,
    coef           = stats::coef(fit),
    r              = fit$theta,
    cb_temp        = if (return_data_mat) cb_temp else NULL,
    X_spline       = if (return_data_mat) X_spline else NULL,
    exposure       = cum_curve,
    lag_effects    = lag_curve,
    lag_importance = lag_importance
  )

  class(out) <- "flexible_fit"
  out
}
