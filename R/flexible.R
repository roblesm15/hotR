# flexible.R
# User-facing interface for the flexible distributed lag model.
#
# This file provides fit_flexible(), a wrapper around a distributed
# lag nonlinear model (DLNM) for estimating the full temperature-
# mortality exposure-response relationship. Unlike fit_hot(), which
# estimates a single threshold parameter under a structured model,
# fit_flexible() makes minimal assumptions about the shape of the
# relationship and the lag structure.
#
# The model assumes:
#
#   log(mu_t) = X_t gamma + f(T_t, l) + log(pop_t)
#
# where f(T_t, l) is a cross-basis function spanning both the
# temperature and lag dimensions, estimated using natural splines
# in the temperature direction and integer lag weights.
#
# The minimum mortality temperature (MMT) is estimated as the
# temperature at which the cumulative exposure-response function
# is minimized. Uncertainty in the MMT is quantified by simulation
# from the joint distribution of the cross-basis coefficients.
#
# This model is intended as a flexible reference against which the
# structured threshold model can be compared. It is better suited
# to longer and less noisy series where the full shape of the
# temperature-mortality relationship can be identified.
#
# Note: fit_flexible() depends on internal functions from the dlnm
# package accessed via :::. This is necessary because dlnm does not
# export the prediction utilities required by findmin(). This
# dependency may break if dlnm changes its internal interface.

# -------------------------------------------------------------------
# Minimum mortality temperature estimation
# -------------------------------------------------------------------

#' Find the minimum mortality temperature from a DLNM fit
#'
#' Estimates the minimum mortality temperature (MMT) as the value of
#' the exposure variable at which the cumulative exposure-response
#' function from a distributed lag nonlinear model is minimized.
#' Optionally quantifies uncertainty by simulation from the joint
#' distribution of cross-basis coefficients.
#'
#' @details
#' The MMT is found by evaluating the cumulative linear predictor
#' \eqn{\hat\eta(x) = \mathbf{X}_\text{pred}(x) \hat\beta} over a
#' fine grid of exposure values and returning the grid point at which
#' \eqn{\hat\eta} is smallest. When \code{sim = TRUE}, the same
#' procedure is repeated \code{nsim} times using coefficient vectors
#' drawn from the multivariate normal approximation to the posterior
#' of the cross-basis coefficients, yielding a distribution of MMT
#' estimates from which confidence intervals can be computed.
#'
#' This function is adapted from the \code{findmin} utility in the
#' \pkg{dlnm} package and relies on several internal \pkg{dlnm}
#' functions via \code{:::}. It is included here because \pkg{dlnm}
#' does not export these utilities.
#'
#' @param basis An object of class \code{"crossbasis"} or
#'   \code{"onebasis"} as produced by \code{dlnm::crossbasis}.
#' @param model A fitted model object (e.g., from \code{MASS::glm.nb})
#'   containing the cross-basis term. Either \code{model} or both
#'   \code{coef} and \code{vcov} must be provided.
#' @param coef Numeric vector of cross-basis coefficients. Used only
#'   if \code{model} is \code{NULL}.
#' @param vcov Numeric matrix. Covariance matrix for \code{coef}.
#'   Used only if \code{model} is \code{NULL}.
#' @param at Numeric vector. Grid of exposure values over which to
#'   search for the minimum. If \code{NULL}, a regular grid from
#'   \code{from} to \code{to} with step \code{by} is used.
#' @param from Numeric scalar. Left endpoint of the search grid.
#'   Defaults to the minimum of the observed exposure range.
#' @param to Numeric scalar. Right endpoint of the search grid.
#'   Defaults to the maximum of the observed exposure range.
#' @param by Numeric scalar. Step size for the search grid. Default
#'   is \code{0.1}.
#' @param sim Logical. If \code{TRUE}, returns a vector of MMT
#'   estimates from \code{nsim} simulated coefficient draws rather
#'   than the point estimate. Default is \code{FALSE}.
#' @param nsim Positive integer. Number of simulated draws used when
#'   \code{sim = TRUE}. Default is \code{5000}.
#'
#' @return If \code{sim = FALSE}, a numeric scalar giving the point
#'   estimate of the MMT. If \code{sim = TRUE}, a numeric vector of
#'   length \code{nsim} giving simulated MMT estimates.
#'
#' @noRd
findmin <- function(basis, model = NULL, coef = NULL, vcov = NULL,
                    at = NULL, from = NULL, to = NULL, by = NULL,
                    sim = FALSE, nsim = 5000) {

  if (!any(class(basis) %in% c("crossbasis", "onebasis"))) {
    stop(
      "the first argument must be an object of class ",
      "'crossbasis' or 'onebasis'"
    )
  }

  one   <- any(class(basis) %in% c("onebasis"))
  range <- attr(basis, "range")

  if (is.null(by)) by <- 0.1
  lag <- if (one) c(0, 0) else attr(basis, "lag")

  if (is.null(model) && (is.null(coef) || is.null(vcov))) {
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  }

  name <- "cb_temp"
  cond <- if (one) {
    paste(name, "[[:print:]]*b[0-9]{1,2}", sep = "")
  } else {
    paste(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}", sep = "")
  }

  if (!is.null(model)) {
    model.class <- class(model)
    coef        <- dlnm:::getcoef(model, model.class)
    ind         <- grep(cond, names(coef))
    coef        <- coef[ind]
    vcov        <- dlnm:::getvcov(model, model.class)[ind, ind, drop = FALSE]
  } else {
    model.class <- NA
  }

  if (length(coef) != ncol(basis) ||
      length(coef) != dim(vcov)[1] ||
      any(is.na(coef)) ||
      any(is.na(vcov))) {
    stop("model or coef/vcov not consistent with basis")
  }

  at      <- dlnm:::mkat(at, from, to, by, range, lag, bylag = 1)
  predvar <- if (is.matrix(at)) rownames(at) else at
  predlag <- dlnm:::seqlag(lag, by = 1)

  type   <- if (one) "one" else "cb"
  Xpred  <- dlnm:::mkXpred(type, basis, at, predvar, predlag, cen = NULL)

  Xpredall <- 0
  for (i in seq(length(predlag))) {
    ind      <- seq(length(predvar)) + length(predvar) * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
  }

  pred <- drop(Xpredall %*% coef)
  ind  <- which.min(pred)
  min  <- predvar[ind]

  if (sim) {
    k      <- length(coef)
    eigen  <- eigen(vcov)
    X      <- matrix(rnorm(length(coef) * nsim), nsim)
    coefsim <- coef +
      eigen$vectors %*% diag(sqrt(eigen$values), k) %*% t(X)

    minsim <- apply(coefsim, 2, function(coefi) {
      pred <- drop(Xpredall %*% coefi)
      predvar[which.min(pred)]
    })
  }

  if (sim) minsim else min
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
#' ## Model
#' Mortality counts \eqn{y_t} are modeled as:
#' \deqn{y_t \sim \mathrm{NegBin}(\mu_t,\, r)}
#' \deqn{\log\mu_t = \mathbf{x}_t^\top\gamma +
#'   f(T_t, l) + \log(\mathrm{pop}_t)}
#'
#' where \eqn{f(T_t, l)} is a cross-basis function constructed by
#' \code{dlnm::crossbasis} using natural splines in the temperature
#' direction (with knots at the 20th, 40th, 60th, and 80th
#' percentiles of temperature by default) and integer lag weights
#' up to \code{lag} days. The baseline \eqn{\mathbf{x}_t^\top\gamma}
#' captures long-term trend and seasonality via a natural spline with
#' \code{df_per_year_spline} degrees of freedom per year.
#'
#' ## Minimum mortality temperature
#' The MMT is estimated as the temperature minimizing the cumulative
#' exposure-response function. A 95\% confidence interval is obtained
#' by simulation from the multivariate normal approximation to the
#' distribution of the cross-basis coefficients when
#' \code{sim_boot = TRUE}.
#'
#' ## Lag importance
#' A data frame summarizing the contribution of each lag to the
#' overall exposure-response is returned as \code{lag_importance}
#' when \code{return_data_mat = TRUE}. Four summaries are provided:
#' integrated absolute log-RR, integrated squared log-RR, mean
#' absolute log-RR, and maximum absolute log-RR across the temperature
#' distribution.
#'
#' ## Comparison with fit_hot
#' \code{fit_flexible} imposes no structural assumption on the shape
#' of the temperature-mortality relationship or the lag structure. It
#' is well suited to longer series (five or more years) where the full
#' shape of the relationship can be identified. For short or noisy
#' series, \code{fit_hot} is preferred.
#'
#' @param city_data Data frame with columns:
#'   \describe{
#'     \item{date}{Date vector ordered chronologically with no gaps.}
#'     \item{temperature}{Numeric. Daily mean or maximum temperature
#'       in degrees C.}
#'     \item{deaths}{Non-negative integer. Daily all-cause or
#'       cause-specific mortality counts.}
#'     \item{population}{Positive numeric. Population at risk, used
#'       as a multiplicative offset on the log scale.}
#'   }
#' @param df_per_year_spline Positive numeric scalar. Degrees of
#'   freedom per year for the natural spline baseline. Default is
#'   \code{3}.
#' @param ktemp Optional numeric vector. Interior knot positions for
#'   the natural spline in the temperature direction of the
#'   cross-basis. If \code{NULL} (default), four knots are placed at
#'   the 20th, 40th, 60th, and 80th percentiles of the observed
#'   temperature distribution.
#' @param lag Non-negative integer. Maximum lag in days included in
#'   the cross-basis. Default is \code{7}.
#' @param at Optional numeric vector. Temperature values at which the
#'   exposure-response curve is evaluated. If \code{NULL} (default),
#'   all unique observed temperatures are used.
#' @param return_data_mat Logical. If \code{TRUE}, the cross-basis
#'   object, design matrix, and detailed exposure-response curves are
#'   included in the returned list. Default is \code{FALSE}.
#' @param sim_boot Logical. If \code{TRUE}, uncertainty in the MMT is
#'   quantified by simulation from the coefficient distribution,
#'   yielding a 95\% confidence interval. Default is \code{TRUE}.
#'
#' @return An object of class \code{"flexible_fit"}: a named list
#'   with components:
#' \describe{
#'   \item{df_loc}{Integer. Degrees of freedom used for the spline
#'     baseline.}
#'   \item{lag}{Integer. Maximum lag used in the cross-basis.}
#'   \item{mmt}{Numeric scalar. Point estimate of the MMT in degrees
#'     C.}
#'   \item{ci}{Named numeric vector of length 2. 95\% simulation-based
#'     confidence interval for the MMT (\code{c(lower, upper)}).
#'     \code{NA} if \code{sim_boot = FALSE}.}
#'   \item{coef}{Named numeric vector of all model coefficients.}
#'   \item{r}{Numeric scalar. Estimated negative binomial dispersion
#'     parameter.}
#'   \item{cb_temp}{The \code{crossbasis} object (if
#'     \code{return_data_mat = TRUE}, else \code{NULL}).}
#'   \item{X_spline}{Spline design matrix without intercept (if
#'     \code{return_data_mat = TRUE}, else \code{NULL}).}
#'   \item{exposure}{Data frame of cumulative exposure-response
#'     estimates with columns \code{temp}, \code{RR}, \code{low},
#'     \code{high} (if \code{return_data_mat = TRUE}, else
#'     \code{NULL}).}
#'   \item{lag_effects}{Data frame of lag-specific exposure-response
#'     estimates with columns \code{temp}, \code{lag}, \code{RR},
#'     \code{low}, \code{high} (if \code{return_data_mat = TRUE},
#'     else \code{NULL}).}
#'   \item{lag_importance}{Data frame summarizing lag contributions
#'     with columns \code{lag}, \code{importance_abs},
#'     \code{importance_sq}, \code{mean_abs_eta},
#'     \code{max_abs_eta} (if \code{return_data_mat = TRUE}, else
#'     \code{NULL}).}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(ahmedabad)
#'
#' fit <- fit_flexible(
#'   city_data          = ahmedabad,
#'   df_per_year_spline = 3,
#'   lag                = 7,
#'   return_data_mat    = TRUE
#' )
#'
#' # Minimum mortality temperature
#' fit$mmt
#'
#' # 95% confidence interval for the MMT
#' fit$ci
#' }
fit_flexible <- function(
    city_data,
    df_per_year_spline = 3,
    ktemp              = NULL,
    lag                = 7,
    at                 = NULL,
    return_data_mat    = FALSE,
    sim_boot           = TRUE
) {

  # --- input checks ---
  required_cols <- c("date", "temperature", "deaths", "population")
  missing_cols  <- setdiff(required_cols, names(city_data))
  if (length(missing_cols) > 0) {
    stop(
      "city_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!inherits(city_data$date, "Date")) {
    stop("city_data$date must be of class Date.")
  }
  if (any(city_data$deaths < 0, na.rm = TRUE)) {
    stop("city_data$deaths must be non-negative.")
  }
  if (any(city_data$population <= 0, na.rm = TRUE)) {
    stop("city_data$population must be strictly positive.")
  }

  # --- prepare data ---
  tnum <- as.numeric(city_data$date)
  temp <- city_data$temperature
  y    <- city_data$deaths
  off  <- log(city_data$population)

  no_years <- (max(city_data$date) - min(city_data$date) + 1) / 365.25
  df_loc   <- max(4, round(df_per_year_spline * no_years))

  # --- temperature knots ---
  if (is.null(ktemp)) {
    r     <- range(temp, na.rm = TRUE)
    ktemp <- r[1] + diff(r) / 5 * 1:4
  }

  # --- cross-basis ---
  cb_temp <- dlnm::crossbasis(
    x      = temp,
    lag    = lag,
    argvar = list(fun = "ns", knots = ktemp),
    arglag = list(fun = "integer")
  )

  # --- baseline spline ---
  X_spline <- model.matrix(~ splines::ns(tnum, df = df_loc))
  X_spline <- X_spline[, -1, drop = FALSE]

  # --- fit ---
  fit <- MASS::glm.nb(
    y ~ X_spline + cb_temp + offset(off),
    link = log
  )

  # --- MMT ---
  if (is.null(at)) at <- sort(unique(temp))

  mmt <- findmin(
    cb_temp,
    model = fit,
    at    = at,
    by    = 0.001,
    sim   = FALSE
  )

  ci <- if (sim_boot) {
    boots <- findmin(
      cb_temp,
      model = fit,
      at    = at,
      by    = 0.001,
      sim   = TRUE,
      nsim  = 10000
    )
    quantile(boots, probs = c(0.025, 0.975))
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }

  # --- exposure-response curves ---
  cp <- dlnm::crosspred(
    cb_temp, fit,
    at     = at,
    bylag  = 1,
    cumul  = TRUE,
    cen    = mmt
  )

  cum_curve <- data.frame(
    temp = as.numeric(names(cp$allRRfit)),
    RR   = cp$allRRfit,
    low  = cp$allRRlow,
    high = cp$allRRhigh
  )

  lag_curve <- local({
    rr    <- cp$matRRfit
    low   <- cp$matRRlow
    hi    <- cp$matRRhigh
    temps <- as.numeric(rownames(rr))
    lags  <- as.integer(gsub("[^0-9]", "", colnames(rr)))
    data.frame(
      temp = rep(temps, times = ncol(rr)),
      lag  = rep(lags,  each  = nrow(rr)),
      RR   = as.vector(rr),
      low  = as.vector(low),
      high = as.vector(hi),
      row.names = NULL
    )
  })

  lag_importance <- local({
    eta   <- cp$matfit
    temps <- as.numeric(rownames(eta))
    lags  <- as.integer(gsub("[^0-9]", "", colnames(eta)))
    w     <- c(diff(temps), tail(diff(temps), 1))
    data.frame(
      lag            = lags,
      importance_abs = as.numeric(apply(abs(eta), 2,
                         function(z) sum(z * w))),
      importance_sq  = as.numeric(apply(eta^2,    2,
                         function(z) sum(z * w))),
      mean_abs_eta   = as.numeric(colMeans(abs(eta))),
      max_abs_eta    = as.numeric(apply(abs(eta), 2, max)),
      row.names      = NULL
    )
  })

  # --- assemble output ---
  out <- list(
    df_loc         = df_loc,
    lag            = lag,
    mmt            = mmt,
    ci             = ci,
    coef           = stats::coef(fit),
    r              = fit$theta,
    cb_temp        = if (return_data_mat) cb_temp    else NULL,
    X_spline       = if (return_data_mat) X_spline   else NULL,
    exposure       = if (return_data_mat) cum_curve  else NULL,
    lag_effects    = if (return_data_mat) lag_curve  else NULL,
    lag_importance = if (return_data_mat) lag_importance else NULL
  )

  class(out) <- "flexible_fit"
  out
}