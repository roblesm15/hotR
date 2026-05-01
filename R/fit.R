# fit.R
# Fitting engine for the structured threshold model (STM).
#
# The model assumes that mortality counts follow a negative binomial
# distribution with mean:
#
#   log(mu_t) = X_t gamma + beta * g(z_t; c, k) + log(pop_t)
#
# where:
#   - X_t        is a row of the natural spline design matrix for time,
#                capturing long-term trend and seasonality
#   - gamma      is the vector of spline coefficients
#   - g(z_t;c,k) is the smooth hinge function evaluated at the lagged
#                temperature exposure z_t = T_{t-L}
#   - beta > 0   is the heat effect coefficient (constrained positive)
#   - c          is the heat-onset risk threshold (HOT), the target
#                parameter of primary interest
#   - k          controls the smoothness of the hinge transition
#   - log(pop_t) is the population offset
#
# Optimization proceeds by block coordinate descent, alternating between:
#   (1) updating (gamma, beta) with c fixed via BFGS
#   (2) updating c with (gamma, beta) fixed via golden-section search
#   (3) periodically updating the dispersion parameter r via ML
#
# beta is constrained to be positive by optimizing over b = log(beta)
# and transforming back. All internal fitting functions operate on b;
# the user-facing wrapper returns beta = exp(b).

# -------------------------------------------------------------------
# Negative log-likelihood and gradient for (gamma, b) with fixed theta
# -------------------------------------------------------------------

#' Negative log-likelihood for regression coefficients at fixed dispersion
#'
#' Evaluates the negative binomial negative log-likelihood as a function
#' of the regression parameter vector \eqn{(\gamma, b)}, where
#' \eqn{\beta = \exp(b)} is the heat effect coefficient. The dispersion
#' parameter \eqn{\theta} (size) and the smooth hinge values \code{g}
#' are treated as fixed.
#'
#' This function is passed to \code{optim} inside
#' \code{fit_nb_coef_fixed_theta} and is not intended to be called
#' directly.
#'
#' @param par Numeric vector of length \eqn{p_\gamma + 1} containing
#'   the spline coefficients \eqn{\gamma} followed by \eqn{b = \log\beta}.
#' @param y Non-negative integer vector of observed death counts.
#' @param off Numeric vector of offsets (log population).
#' @param X_base Numeric matrix of dimension \eqn{n \times p_\gamma}.
#'   Design matrix including the intercept and natural spline columns
#'   for time.
#' @param g Numeric vector of length \eqn{n}. Pre-computed smooth hinge
#'   values \eqn{g(z_t; c, k)} at the current threshold \eqn{c}.
#' @param theta Positive numeric scalar. Negative binomial dispersion
#'   parameter, held fixed during this optimization step.
#'
#' @return Scalar. Negative log-likelihood value.
#'
#' @noRd
negloglik_gamma_b_fixed_theta <- function(par, y, off, X_base, g, theta) {
  p_gamma <- ncol(X_base)
  gamma   <- par[1:p_gamma]
  b       <- par[p_gamma + 1]
  beta    <- exp(b)
  eta_lin <- as.vector(X_base %*% gamma + beta * g + off)
  -nb_loglik_full(eta_lin, y, theta)
}

#' Gradient of the negative log-likelihood for (gamma, b) at fixed dispersion
#'
#' Computes the analytic gradient of the negative log-likelihood with
#' respect to \eqn{(\gamma, b)} where \eqn{\beta = \exp(b)}. Used
#' alongside \code{negloglik_gamma_b_fixed_theta} in the BFGS step of
#' the block coordinate descent.
#'
#' The chain rule gives:
#' \deqn{\frac{\partial \ell}{\partial b} =
#'   \sum_t u_t \cdot \beta \cdot g(z_t; c, k)}
#' where \eqn{u_t} is the score of the log-likelihood with respect to
#' \eqn{\eta_t} as computed by \code{nb_score_eta}.
#'
#' @inheritParams negloglik_gamma_b_fixed_theta
#'
#' @return Numeric vector of length \eqn{p_\gamma + 1} containing the
#'   gradient with respect to \eqn{(\gamma, b)}.
#'
#' @noRd
grad_gamma_b_fixed_theta <- function(par, y, off, X_base, g, theta) {
  p_gamma <- ncol(X_base)
  gamma   <- par[1:p_gamma]
  b       <- par[p_gamma + 1]
  beta    <- exp(b)

  eta_lin    <- as.vector(X_base %*% gamma + beta * g + off)
  u          <- nb_score_eta(eta_lin, y, theta)

  grad_gamma <- -colSums(X_base * u)
  grad_b     <- -sum(u * beta * g)

  c(grad_gamma, grad_b)
}

# -------------------------------------------------------------------
# BFGS step: fit (gamma, b) for fixed theta and fixed c
# -------------------------------------------------------------------

#' Fit regression coefficients at fixed dispersion and threshold
#'
#' Estimates the spline coefficients \eqn{\gamma} and the log heat
#' effect \eqn{b = \log\beta} by maximizing the negative binomial
#' log-likelihood via BFGS, treating the dispersion parameter
#' \eqn{\theta} and the pre-computed hinge values \code{g} as fixed.
#'
#' @details
#' Starting values are obtained by fitting a \code{glm} with
#' \code{MASS::negative.binomial} family, treating \eqn{\beta} as an
#' unrestricted coefficient and then transforming to the log scale.
#' If \code{start} is provided it is used directly, which allows warm
#' starts across iterations of the block coordinate descent.
#'
#' The positivity constraint \eqn{\beta > 0} is enforced by
#' parameterizing the optimization over \eqn{b = \log\beta} and
#' transforming back after convergence.
#'
#' @inheritParams negloglik_gamma_b_fixed_theta
#' @param start Optional numeric vector of length \eqn{p_\gamma + 1}
#'   providing starting values for \eqn{(\gamma, b)}. If \code{NULL}
#'   (default), starting values are derived from a GLM fit.
#'
#' @return A named list with components:
#' \describe{
#'   \item{coef}{Numeric vector \eqn{(\hat\gamma, \hat b)} of length
#'     \eqn{p_\gamma + 1}.}
#'   \item{gamma}{Numeric vector \eqn{\hat\gamma} of spline coefficients.}
#'   \item{b}{Numeric scalar \eqn{\hat b = \log\hat\beta}.}
#'   \item{beta}{Numeric scalar \eqn{\hat\beta = \exp(\hat b)}.}
#'   \item{eta_lin}{Numeric vector of fitted linear predictors.}
#'   \item{mu}{Numeric vector of fitted means \eqn{\exp(\hat\eta)}.}
#' }
#'
#' @noRd
fit_nb_coef_fixed_theta <- function(y, off, X_base, g, theta, start = NULL) {
  p_gamma <- ncol(X_base)

  if (is.null(start)) {
    fam  <- MASS::negative.binomial(theta = theta, link = "log")
    fit0 <- glm.fit(
      x         = cbind(X_base, g),
      y         = y,
      family    = fam,
      offset    = off,
      intercept = FALSE
    )
    coef0  <- fit0$coefficients
    gamma0 <- coef0[1:p_gamma]
    beta0  <- max(1e-8, coef0[p_gamma + 1])
    start  <- c(gamma0, log(beta0))
  } else {
    start <- as.numeric(start)
    if (length(start) != p_gamma + 1) stop("start has wrong length")
  }

  opt <- optim(
    par     = start,
    fn      = negloglik_gamma_b_fixed_theta,
    gr      = grad_gamma_b_fixed_theta,
    y       = y,
    off     = off,
    X_base  = X_base,
    g       = g,
    theta   = theta,
    method  = "BFGS",
    control = list(maxit = 200, reltol = 1e-10)
  )

  par_hat   <- opt$par
  gamma_hat <- par_hat[1:p_gamma]
  b_hat     <- par_hat[p_gamma + 1]
  beta_hat  <- exp(b_hat)

  eta_lin <- as.vector(X_base %*% gamma_hat + beta_hat * g + off)
  mu      <- exp(eta_lin)

  list(
    coef  = c(gamma_hat, b_hat),
    gamma = gamma_hat,
    b     = b_hat,
    beta  = beta_hat,
    eta_lin = eta_lin,
    mu    = mu
  )
}

# -------------------------------------------------------------------
# Dispersion update
# -------------------------------------------------------------------

#' Update the negative binomial dispersion parameter by maximum likelihood
#'
#' Re-estimates the dispersion (size) parameter \eqn{r} of the negative
#' binomial distribution using \code{MASS::theta.ml}, treating the
#' current fitted means \code{mu} as fixed. This profile ML step is
#' called periodically within the block coordinate descent to avoid
#' the computational cost of updating \eqn{r} at every iteration.
#'
#' @details
#' If \code{MASS::theta.ml} throws an error or returns a non-finite or
#' non-positive value, the function returns \code{theta_start} unchanged.
#' This fallback ensures robustness in early iterations when the model
#' may not yet be close to the optimum.
#'
#' @param y Non-negative integer vector of observed death counts.
#' @param mu Positive numeric vector of current fitted means.
#' @param theta_start Positive numeric scalar. Current dispersion estimate,
#'   used as a fallback if the update fails.
#' @param limit Integer. Maximum number of iterations passed to
#'   \code{MASS::theta.ml}. Default is \code{25}.
#' @param eps Numeric scalar. Convergence tolerance for
#'   \code{MASS::theta.ml}. Default is \code{1e-8}.
#'
#' @return Positive numeric scalar. Updated dispersion estimate, or
#'   \code{theta_start} if the update failed.
#'
#' @noRd
refresh_theta_ml <- function(y, mu, theta_start, limit = 25, eps = 1e-8) {
  th <- tryCatch(
    MASS::theta.ml(
      y       = y,
      mu      = mu,
      n       = length(y),
      weights = rep(1, length(y)),
      limit   = limit,
      eps     = eps,
      trace   = FALSE
    ),
    error = function(e) theta_start
  )

  if (!is.finite(th) || th <= 0) theta_start else th
}

# -------------------------------------------------------------------
# Golden-section step: optimize c with (gamma, b) fixed
# -------------------------------------------------------------------

#' Conditional negative log-likelihood for the threshold parameter
#'
#' Evaluates the negative binomial negative log-likelihood as a function
#' of the threshold \eqn{c} alone, with all other parameters
#' \eqn{(\gamma, b, r)} held fixed. Used as the objective function
#' in the golden-section search over \eqn{c} within the block
#' coordinate descent.
#'
#' @param c Numeric scalar. Candidate threshold value.
#' @param gamma Numeric vector of spline coefficients, held fixed.
#' @param b Numeric scalar. Log heat effect \eqn{b = \log\beta},
#'   held fixed.
#' @param r Positive numeric scalar. Negative binomial dispersion,
#'   held fixed.
#' @param X_base Numeric matrix. Design matrix including intercept
#'   and spline columns for time.
#' @param z Numeric vector of lagged temperature exposures.
#' @param off Numeric vector of offsets (log population).
#' @param y Non-negative integer vector of observed death counts.
#' @param k Numeric scalar. Smoothness parameter for the hinge
#'   transition.
#'
#' @return Scalar. Negative log-likelihood value, or \code{Inf} if
#'   \code{c} is not finite.
#'
#' @noRd
negloglik_c_cond_softmax <- function(c, gamma, b, r, X_base, z, off, y, k) {
  if (!is.finite(c)) return(Inf)
  beta    <- exp(b)
  g       <- smooth_hinge(z, c, k)
  eta_lin <- as.vector(X_base %*% gamma + beta * g + off)
  -nb_loglik_full(eta_lin, y, r)
}

#' Conditional negative log-likelihood for lag weights (stub)
#'
#' Returns 0. Retained for interface compatibility with the softmax
#' lag-weight parameterization. In the single-lag model there are no
#' lag-weight parameters to optimize.
#'
#' @noRd
negloglik_eta_cond_softmax <- function(eta, gamma, beta, r, c,
                                       X_base, temp_lag_mat, off, y, k) 0

#' Gradient of the conditional log-likelihood for lag weights (stub)
#'
#' Returns an empty numeric vector. Retained for interface
#' compatibility with the softmax lag-weight parameterization.
#'
#' @noRd
grad_eta_cond_softmax <- function(eta, gamma, beta, r, c,
                                  X_base, temp_lag_mat, off, y, k) {
  numeric(0)
}

# -------------------------------------------------------------------
# Block coordinate descent: main fitting engine
# -------------------------------------------------------------------

#' Block coordinate descent for the structured threshold model
#'
#' Fits the structured threshold model (STM) by iterating between
#' three update steps until convergence:
#'
#' \enumerate{
#'   \item **Coefficient step**: maximize the log-likelihood over
#'     \eqn{(\gamma, b)} with \eqn{c} and \eqn{r} fixed, using BFGS
#'     via \code{fit_nb_coef_fixed_theta}.
#'   \item **Threshold step**: maximize the log-likelihood over
#'     \eqn{c} with \eqn{(\gamma, b, r)} fixed, using golden-section
#'     search via \code{stats::optimize}.
#'   \item **Dispersion step** (every \code{theta_update_every}
#'     iterations): update \eqn{r} by profile ML via
#'     \code{refresh_theta_ml}.
#' }
#'
#' Convergence is declared when either the absolute change in
#' log-likelihood falls below \code{tol} or the absolute change in
#' \eqn{c} falls below \code{tol_c}.
#'
#' @details
#' The search interval for \eqn{c} is restricted to
#' \eqn{[Q(\code{lower\_q}),\, Q(\code{upper\_q})]} of the observed
#' lagged temperature exposure \eqn{z}. This prevents the threshold
#' from being estimated in regions of the temperature distribution
#' where there is insufficient data to identify a heat effect, which
#' is a particular concern for the short Indian city mortality series.
#'
#' @param y Non-negative integer vector of observed death counts.
#' @param off Numeric vector of offsets (log population).
#' @param X_spline Numeric matrix of natural spline columns for time
#'   (without intercept). The intercept is added internally.
#' @param temp_lag_mat Numeric vector or matrix of lagged temperature
#'   exposures, as produced by \code{extract_single_lag} or
#'   \code{build_temp_lag_matrix}.
#' @param c0 Optional numeric scalar. Initial value for the threshold
#'   \eqn{c}. If \code{NULL} (default), the median of \eqn{z} is used.
#' @param eta0 Ignored. Included for interface compatibility.
#' @param k Numeric scalar. Smoothness parameter for the hinge
#'   transition. Default is \code{5}.
#' @param lower_q Numeric scalar in \eqn{(0, 1)}. Lower quantile of
#'   \eqn{z} defining the left boundary of the search interval for
#'   \eqn{c}. Default is \code{0.10}.
#' @param upper_q Numeric scalar in \eqn{(0, 1)}. Upper quantile of
#'   \eqn{z} defining the right boundary of the search interval for
#'   \eqn{c}. Default is \code{0.99}.
#' @param tol Positive numeric scalar. Convergence tolerance on the
#'   absolute change in log-likelihood. Default is \code{1e-8}.
#' @param tol_c Positive numeric scalar. Convergence tolerance on the
#'   absolute change in \eqn{c}. Default is \code{1e-5}.
#' @param tol_eta Ignored. Included for interface compatibility.
#' @param maxit Positive integer. Maximum number of block coordinate
#'   descent iterations. Default is \code{100}.
#' @param theta_update_every Positive integer. Number of iterations
#'   between dispersion updates. Default is \code{10}.
#' @param eta_maxit Ignored. Included for interface compatibility.
#' @param eta_reltol Ignored. Included for interface compatibility.
#' @param verbose Logical. If \code{TRUE}, prints iteration-level
#'   diagnostics including log-likelihood, \eqn{\hat c}, \eqn{\hat\beta},
#'   and \eqn{\hat r}. Default is \code{FALSE}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{converged}{Logical. Whether the algorithm converged within
#'     \code{maxit} iterations.}
#'   \item{it}{Integer. Number of iterations performed.}
#'   \item{c}{Numeric scalar. Estimated threshold \eqn{\hat c} (the HOT).}
#'   \item{gamma}{Numeric vector. Estimated spline coefficients
#'     \eqn{\hat\gamma}.}
#'   \item{b}{Numeric scalar. Estimated \eqn{\hat b = \log\hat\beta}.}
#'   \item{beta}{Numeric scalar. Estimated heat effect
#'     \eqn{\hat\beta = \exp(\hat b)}.}
#'   \item{r}{Numeric scalar. Estimated dispersion parameter \eqn{\hat r}.}
#'   \item{eta}{Empty numeric vector. Retained for compatibility.}
#'   \item{w}{Numeric scalar 1. Lag weight. Retained for compatibility.}
#'   \item{z}{Numeric vector. Lagged temperature exposure used in fitting.}
#'   \item{logLik}{Numeric scalar. Log-likelihood at convergence.}
#' }
#'
#' @noRd
joint_nb_threshold_avg_lag_softmax <- function(
    y, off, X_spline, temp_lag_mat,
    c0 = NULL,
    eta0 = NULL,
    k = 5,
    lower_q = 0.10,
    upper_q = 0.99,
    tol = 1e-8,
    tol_c = 1e-5,
    tol_eta = 1e-5,
    maxit = 100,
    theta_update_every = 10,
    eta_maxit = 10,
    eta_reltol = 1e-4,
    verbose = FALSE
) {
  X_base <- cbind(1, X_spline)
  storage.mode(X_base) <- "double"

  z       <- compute_z_softmax(temp_lag_mat)
  p_gamma <- ncol(X_base)

  c <- if (is.null(c0)) {
    as.numeric(stats::median(z, na.rm = TRUE))
  } else {
    as.numeric(c0)
  }

  c_lo_global <- as.numeric(stats::quantile(z, lower_q, na.rm = TRUE))
  c_hi_global <- as.numeric(stats::quantile(z, upper_q, na.rm = TRUE))
  c <- min(c_hi_global, max(c_lo_global, c))

  g    <- smooth_hinge(z, c, k)
  fit0 <- MASS::glm.nb(y ~ X_spline + g + offset(off), link = log)

  coefs  <- stats::coef(fit0)
  beta0  <- max(1e-8, as.numeric(coefs["g"]))
  gamma  <- as.numeric(coefs[setdiff(names(coefs), "g")])
  b      <- log(beta0)
  r      <- fit0$theta

  coef_start <- c(gamma, b)
  eta_lin    <- as.vector(X_base %*% gamma + exp(b) * g + off)
  mu         <- exp(eta_lin)
  ll         <- nb_loglik_full(eta_lin, y, r)

  converged <- FALSE
  it        <- 0L

  for (it in seq_len(maxit)) {
    ll_old <- ll
    c_old  <- c

    g <- smooth_hinge(z, c, k)

    fit_coef   <- fit_nb_coef_fixed_theta(
      y = y, off = off, X_base = X_base,
      g = g, theta = r, start = coef_start
    )
    coef_start <- fit_coef$coef
    gamma      <- fit_coef$gamma
    b          <- fit_coef$b
    beta       <- fit_coef$beta
    eta_lin    <- fit_coef$eta_lin
    mu         <- fit_coef$mu
    ll         <- nb_loglik_full(eta_lin, y, r)

    if (it %% theta_update_every == 0) {
      r  <- refresh_theta_ml(y = y, mu = mu, theta_start = r)
      ll <- nb_loglik_full(eta_lin, y, r)
    }

    opt_c <- optimize(
      f        = negloglik_c_cond_softmax,
      interval = c(c_lo_global, c_hi_global),
      gamma    = gamma,
      b        = b,
      r        = r,
      X_base   = X_base,
      z        = z,
      off      = off,
      y        = y,
      k        = k
    )
    c <- as.numeric(opt_c$minimum)

    g       <- smooth_hinge(z, c, k)
    eta_lin <- as.vector(X_base %*% gamma + beta * g + off)
    mu      <- exp(eta_lin)
    ll      <- nb_loglik_full(eta_lin, y, r)

    if (verbose) {
      cat(sprintf(
        "it=%d  ll=%.6f  c=%.4f  beta=%.6g  r=%.6g\n",
        it, ll, c, beta, r
      ))
    }

    if (abs(ll - ll_old) < tol || abs(c - c_old) < tol_c) {
      converged <- TRUE
      break
    }
  }

  g        <- smooth_hinge(z, c, k)
  fit_coef <- fit_nb_coef_fixed_theta(
    y = y, off = off, X_base = X_base,
    g = g, theta = r, start = coef_start
  )

  gamma   <- fit_coef$gamma
  b       <- fit_coef$b
  beta    <- fit_coef$beta
  eta_lin <- fit_coef$eta_lin
  mu      <- fit_coef$mu

  r  <- refresh_theta_ml(y = y, mu = mu, theta_start = r)
  ll <- nb_loglik_full(eta_lin, y, r)

  list(
    converged = converged,
    it        = it,
    c         = c,
    gamma     = gamma,
    b         = b,
    beta      = beta,
    r         = r,
    eta       = numeric(0),
    w         = 1,
    z         = z,
    logLik    = ll
  )
}