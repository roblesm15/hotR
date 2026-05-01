# lag.R
# Utilities for constructing lagged temperature exposures.
# The structured threshold model operates on a single chosen lag L,
# meaning mortality at time t is linked to temperature at time t - L.
# These functions build the full lag matrix and extract the relevant
# column. Compatibility stubs for the softmax lag-weight parameterization
# are retained for internal consistency but are not part of the
# user-facing API.

# -------------------------------------------------------------------
# Lag matrix construction
# -------------------------------------------------------------------

#' Build a matrix of lagged temperature values
#'
#' Constructs an \eqn{n \times (L + 1)} matrix where column \eqn{\ell + 1}
#' contains the temperature series lagged by \eqn{\ell} days, for
#' \eqn{\ell = 0, 1, \ldots, L}. Missing values at the start of each
#' lagged series (before enough observations are available) are filled
#' with the first observed temperature value.
#'
#' @details
#' The padding strategy — filling early rows with \code{temp[1]} rather
#' than \code{NA} — ensures the lag matrix has no missing values and the
#' model can be fit on the full observed series. This is appropriate when
#' the time series is short and dropping initial observations would
#' meaningfully reduce the sample size, as is typical for the Indian city
#' mortality series that motivate this package.
#'
#' @param temp Numeric vector of daily temperature values of length
#'   \eqn{n}, ordered chronologically.
#' @param L Non-negative integer. Maximum lag to include. \code{L = 0}
#'   returns a single-column matrix containing \code{temp} itself.
#'
#' @return An \eqn{n \times (L + 1)} numeric matrix. Column \eqn{j}
#'   (1-indexed) contains \code{temp} lagged by \eqn{j - 1} days.
#'
#' @examples
#' \dontrun{
#' temp <- c(28, 30, 33, 31, 29, 35, 36)
#' build_temp_lag_matrix(temp, L = 2)
#' }
#'
#' @noRd
build_temp_lag_matrix <- function(temp, L) {
  n   <- length(temp)
  out <- matrix(temp[1], nrow = n, ncol = L + 1)
  out[, 1] <- temp
  if (L >= 1) {
    for (ell in 1:L) {
      out[(ell + 1):n, ell + 1] <- temp[1:(n - ell)]
    }
  }
  out
}

#' Extract a single lagged temperature series
#'
#' A convenience wrapper around \code{build_temp_lag_matrix} that
#' returns only the column corresponding to lag \code{L}. This is
#' the exposure vector \eqn{z_t = T_{t-L}} that enters the smooth
#' hinge function in the structured threshold model.
#'
#' @inheritParams build_temp_lag_matrix
#'
#' @return Numeric vector of length \eqn{n} containing temperature
#'   lagged by \code{L} days.
#'
#' @examples
#' \dontrun{
#' temp <- c(28, 30, 33, 31, 29, 35, 36)
#' extract_single_lag(temp, L = 1)
#' }
#'
#' @noRd
extract_single_lag <- function(temp, L) {
  build_temp_lag_matrix(temp, L)[, L + 1]
}

# -------------------------------------------------------------------
# Softmax lag-weight compatibility stubs
# -------------------------------------------------------------------
# The model previously supported estimation of a full distribution of
# lag weights via a softmax parameterization. The current implementation
# uses a single fixed lag L selected by the user (or by cross-validation
# in the analysis scripts). These stubs preserve the internal interface
# so that the fitting and covariance routines do not need branching logic.

#' Softmax lag weight (stub)
#'
#' Returns 1, reflecting that the current model assigns all weight to
#' a single lag. Retained for interface compatibility with the fitting
#' engine.
#'
#' @param eta Ignored. Numeric vector of softmax parameters (empty in
#'   the single-lag model).
#'
#' @return Numeric scalar 1.
#'
#' @noRd
eta_to_w_softmax <- function(eta) 1

#' Initialize softmax parameters (stub)
#'
#' Returns an empty numeric vector, reflecting that no lag-weight
#' parameters are estimated in the single-lag model.
#'
#' @param L Non-negative integer. Maximum lag (unused).
#' @param sd Numeric scalar. Standard deviation for random
#'   initialization (unused).
#'
#' @return Empty numeric vector.
#'
#' @noRd
init_eta_softmax_random <- function(L, sd = 1) numeric(0)

#' Compute the effective temperature exposure (stub)
#'
#' Extracts the single-lag temperature exposure from the lag matrix.
#' In the softmax model this would compute a weighted average across
#' lags; here it simply returns the last column of the lag matrix,
#' which corresponds to lag \code{L}.
#'
#' @param temp_lag_mat Numeric matrix produced by
#'   \code{build_temp_lag_matrix}, or a numeric vector when \code{L = 0}.
#' @param eta Ignored. Included for interface compatibility.
#'
#' @return Numeric vector of length \eqn{n}.
#'
#' @noRd
compute_z_softmax <- function(temp_lag_mat, eta = numeric(0)) {
  if (is.null(dim(temp_lag_mat))) return(as.numeric(temp_lag_mat))
  as.numeric(temp_lag_mat[, ncol(temp_lag_mat)])
}

#' Jacobian of softmax weights with respect to eta (stub)
#'
#' Returns a \eqn{1 \times 0} matrix, reflecting that there are no
#' lag-weight parameters in the single-lag model. Used in the delta
#' method transformation of the covariance matrix.
#'
#' @param eta Ignored.
#'
#' @return A \eqn{1 \times 0} numeric matrix.
#'
#' @noRd
jacobian_w_wrt_eta_softmax <- function(eta) matrix(0, nrow = 1, ncol = 0)

#' Hessian of softmax weights with respect to eta (stub)
#'
#' Returns a \eqn{1 \times 0 \times 0} array, reflecting that there
#' are no lag-weight parameters in the single-lag model.
#'
#' @param eta Ignored.
#'
#' @return A \eqn{1 \times 0 \times 0} numeric array.
#'
#' @noRd
hessian_w_wrt_eta_softmax <- function(eta) array(0, dim = c(1, 0, 0))