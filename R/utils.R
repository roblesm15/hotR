# utils.R
# Core mathematical utilities for the structured threshold model.
# These functions implement the smooth hinge basis function and its
# derivatives, along with negative binomial likelihood components
# used throughout the fitting and inference routines.

# -------------------------------------------------------------------
# Smooth hinge function and derivatives
# -------------------------------------------------------------------

#' Smooth hinge function
#'
#' Computes a smooth approximation to the one-sided quadratic hinge
#' \eqn{(x - c)_+^2 = \max(x - c, 0)^2} using a logistic transition.
#' The parameter \code{k} controls the sharpness of the transition around
#' the threshold \code{c}: larger values of \code{k} make the function
#' approach the hard hinge, while smaller values produce a smoother bend.
#'
#' The function is defined as:
#' \deqn{g(x; c, k) = (x - c)^2 \cdot \sigma(k(x - c))}
#' where \eqn{\sigma} is the logistic function.
#'
#' @param x Numeric vector of observed values (e.g., lagged temperature).
#' @param c Numeric scalar. Threshold parameter (the heat-onset risk
#'   threshold, HOT).
#' @param k Numeric scalar. Smoothness parameter controlling the sharpness
#'   of the hinge transition. Default in the model is \code{k = 5}.
#'
#' @return Numeric vector of the same length as \code{x}.
#'
#' @noRd
smooth_hinge <- function(x, c, k) {
  z <- x - c
  s <- plogis(k * z)
  z^2 * s
}

#' First derivative of the smooth hinge with respect to the threshold
#'
#' Computes \eqn{\partial g / \partial c} where \eqn{g} is the smooth
#' hinge function defined in \code{smooth_hinge}. Used in gradient
#' computations for the threshold parameter during optimization and
#' in the score and Hessian matrices for inference.
#'
#' @inheritParams smooth_hinge
#'
#' @return Numeric vector of the same length as \code{x}.
#'
#' @noRd
dg_dc <- function(x, c, k) {
  z  <- x - c
  s  <- plogis(k * z)
  s1 <- k * s * (1 - s)
  -(2 * z * s + z^2 * s1)
}

#' Second derivative of the smooth hinge with respect to the threshold
#'
#' Computes \eqn{\partial^2 g / \partial c^2} where \eqn{g} is the smooth
#' hinge function defined in \code{smooth_hinge}. Used in the Hessian
#' matrix for the threshold parameter, which enters the observed
#' information matrix used for Wald-type inference on the HOT.
#'
#' @inheritParams smooth_hinge
#'
#' @return Numeric vector of the same length as \code{x}.
#'
#' @noRd
d2g_dc2 <- function(x, c, k) {
  z  <- x - c
  s  <- plogis(k * z)
  s1 <- k * s * (1 - s)
  s2 <- k^2 * s * (1 - s) * (1 - 2 * s)
  2 * s + 4 * z * s1 + z^2 * s2
}

# -------------------------------------------------------------------
# Negative binomial likelihood components
# -------------------------------------------------------------------

#' Negative binomial log-likelihood on the linear predictor scale
#'
#' Evaluates the negative binomial log-likelihood given a vector of
#' linear predictors \code{eta}, observed counts \code{y}, and
#' dispersion parameter \code{r}. The mean is \eqn{\mu = \exp(\eta)}.
#'
#' The log-likelihood contribution for observation \eqn{i} follows the
#' NB2 parameterization used by \code{MASS::glm.nb}:
#' \deqn{\ell_i = \log \Gamma(y_i + r) - \log \Gamma(r) - \log(y_i!) +
#'   r \log\frac{r}{r + \mu_i} + y_i \log\frac{\mu_i}{r + \mu_i}}
#'
#' @param eta Numeric vector of linear predictors (log-scale means
#'   including offset).
#' @param y Non-negative integer vector of observed counts (deaths).
#' @param r Positive numeric scalar. Negative binomial dispersion
#'   (size) parameter. Larger values correspond to less overdispersion.
#'
#' @return Scalar. Total log-likelihood summed over all observations.
#'
#' @noRd
nb_loglik_full <- function(eta, y, r) {
  mu <- exp(eta)
  sum(dnbinom(y, size = r, mu = mu, log = TRUE))
}

#' Score of the negative binomial log-likelihood with respect to eta
#'
#' Computes the observation-level score (first derivative of the
#' log-likelihood with respect to the linear predictor \eqn{\eta})
#' under the NB2 parameterization. This is used in gradient-based
#' optimization of regression coefficients and in the construction
#' of the sandwich covariance estimator.
#'
#' The score for observation \eqn{i} is:
#' \deqn{\frac{\partial \ell_i}{\partial \eta_i} =
#'   y_i - \frac{(y_i + r)\mu_i}{r + \mu_i}}
#'
#' @param eta_lin Numeric vector of linear predictors.
#' @param y Non-negative integer vector of observed counts.
#' @param r Positive numeric scalar. Negative binomial dispersion
#'   parameter.
#'
#' @return Numeric vector of the same length as \code{y} containing
#'   the per-observation score contributions.
#'
#' @noRd
nb_score_eta <- function(eta_lin, y, r) {
  mu <- exp(eta_lin)
  y - (y + r) * mu / (r + mu)
}

#' @importFrom stats dnbinom fitted glm.fit model.matrix optim optimize
#'   plogis predict qnorm quantile rnorm
#' @importFrom utils tail
"_PACKAGE"