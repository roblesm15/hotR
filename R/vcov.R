# vcov.R
# Inference for the structured threshold model (STM).
#
# The internal parameter vector is (gamma, b, c, phi) where:
#   - gamma  : spline coefficients (length p_gamma)
#   - b      : log heat effect, b = log(beta)
#   - c      : heat-onset risk threshold (HOT), the target parameter
#   - phi    : log dispersion, phi = log(r)
#
# Two covariance estimators are computed:
#
#   Wald (model-based):
#     V_wald = [ -H(theta_hat) ]^{-1}
#     where H is the observed Hessian of the log-likelihood.
#     Valid under correct model specification.
#
#   Sandwich (robust):
#     V_sand = V_wald %*% J %*% V_wald
#     where J = S^T S is the outer product of the score matrix.
#     Provides robustness to overdispersion misspecification and
#     temporal dependence in the mortality series.
#
# Standard errors for beta = exp(b) are obtained by the delta method:
#   se(beta) = beta * se(b)
#
# Standard errors for c are read directly from V[idx_c, idx_c] and
# used to construct Wald confidence intervals for the HOT.

# -------------------------------------------------------------------
# Score matrix
# -------------------------------------------------------------------

#' Observation-level score matrix for the structured threshold model
#'
#' Computes the \eqn{n \times (p_\gamma + 3)} matrix of score
#' contributions, where each row contains the gradient of the
#' log-likelihood for a single observation with respect to the
#' internal parameter vector \eqn{(\gamma, b, c, \phi)}.
#'
#' @details
#' The columns of the score matrix correspond to:
#' \enumerate{
#'   \item Columns \eqn{1, \ldots, p_\gamma}: score with respect to
#'     the spline coefficients \eqn{\gamma}.
#'   \item Column \eqn{p_\gamma + 1}: score with respect to
#'     \eqn{b = \log\beta}. By the chain rule this equals
#'     \eqn{u_t \cdot \beta \cdot g(z_t; c, k)}, where \eqn{u_t}
#'     is the score with respect to \eqn{\eta_t}.
#'   \item Column \eqn{p_\gamma + 2}: score with respect to the
#'     threshold \eqn{c}. Equals
#'     \eqn{u_t \cdot \beta \cdot \partial g / \partial c}.
#'   \item Column \eqn{p_\gamma + 3}: score with respect to
#'     \eqn{\phi = \log r}. Derived from the negative binomial
#'     log-likelihood differentiated with respect to \eqn{r},
#'     scaled by \eqn{r} via the chain rule.
#' }
#'
#' The score matrix is the primary ingredient of the sandwich
#' covariance estimator. Each row being an independent observation
#' score contribution means the outer product \eqn{J = S^\top S}
#' is a consistent estimator of the expected Fisher information
#' under mild misspecification.
#'
#' @param gamma Numeric vector of length \eqn{p_\gamma}. Estimated
#'   spline coefficients.
#' @param beta Positive numeric scalar. Estimated heat effect
#'   \eqn{\hat\beta = \exp(\hat b)}.
#' @param c Numeric scalar. Estimated threshold \eqn{\hat c}.
#' @param phi Numeric scalar. Log dispersion \eqn{\phi = \log r}.
#' @param eta Ignored. Included for interface compatibility.
#' @param X_full Numeric matrix of dimension \eqn{n \times p_\gamma}.
#'   Full design matrix including intercept and spline columns.
#' @param temp_lag_mat Numeric vector or matrix of lagged temperature
#'   exposures. Only the last column is used (corresponding to lag
#'   \eqn{L}).
#' @param y Non-negative integer vector of observed death counts.
#' @param off Numeric vector of offsets (log population).
#' @param k Numeric scalar. Smoothness parameter for the hinge
#'   transition. Default is \code{5}.
#'
#' @return Numeric matrix of dimension \eqn{n \times (p_\gamma + 3)}.
#'   Columns are ordered as \eqn{(\gamma, b, c, \phi)}.
#'
#' @noRd
score_matrix_nb_threshold_lag_softmax <- function(gamma, beta, c, phi, eta,
                                                  X_full, temp_lag_mat, y, off,
                                                  k = 5) {
  r <- exp(phi)
  z <- as.numeric(temp_lag_mat)
  g <- smooth_hinge(z, c, k)

  eta_lin <- as.vector(X_full %*% gamma + beta * g + off)
  mu      <- exp(eta_lin)

  u <- y - (y + r) * mu / (r + mu)

  a_r <- digamma(y + r) - digamma(r) +
    log(r) + 1 - log(r + mu) - (y + r) / (r + mu)

  S_gamma <- X_full * u
  S_b     <- u * (beta * g)
  S_c     <- u * (beta * dg_dc(z, c, k))
  S_phi   <- a_r * r

  cbind(S_gamma, S_b, S_c, S_phi)
}

# -------------------------------------------------------------------
# Hessian matrix
# -------------------------------------------------------------------

#' Observed Hessian of the negative binomial log-likelihood
#'
#' Computes the \eqn{(p_\gamma + 3) \times (p_\gamma + 3)} observed
#' Hessian of the log-likelihood with respect to the internal
#' parameter vector \eqn{(\gamma, b, c, \phi)}.
#'
#' @details
#' The Hessian is structured into blocks corresponding to pairs of
#' parameter subvectors. Let \eqn{u_t} denote the score with respect
#' to \eqn{\eta_t}, and \eqn{u_t^{(2)}} its second derivative
#' component:
#' \deqn{u_t^{(2)} = -\frac{(y_t + r) r \mu_t}{(r + \mu_t)^2}}
#'
#' Key diagonal blocks:
#' \itemize{
#'   \item \eqn{H_{\gamma\gamma} = X^\top \mathrm{diag}(u^{(2)}) X}
#'   \item \eqn{H_{bb}} includes a term \eqn{u_t \cdot \beta \cdot g}
#'     from the second derivative of \eqn{\exp(b)}.
#'   \item \eqn{H_{cc}} includes a term involving the second derivative
#'     of the smooth hinge \eqn{\partial^2 g / \partial c^2}.
#'   \item \eqn{H_{\phi\phi}} uses trigamma functions from
#'     differentiating the negative binomial log-likelihood twice
#'     with respect to \eqn{r}.
#' }
#'
#' The matrix is symmetrized as \eqn{(H + H^\top) / 2} before
#' returning to guard against floating-point asymmetry.
#'
#' @inheritParams score_matrix_nb_threshold_lag_softmax
#'
#' @return Symmetric numeric matrix of dimension
#'   \eqn{(p_\gamma + 3) \times (p_\gamma + 3)}.
#'
#' @noRd
hessian_nb_threshold_lag_softmax <- function(gamma, beta, c, phi, eta,
                                             X_full, temp_lag_mat, y, off,
                                             k = 5) {
  r <- exp(phi)

  z   <- as.numeric(temp_lag_mat)
  g   <- smooth_hinge(z, c, k)
  gc  <- dg_dc(z, c, k)
  gcc <- d2g_dc2(z, c, k)

  eta_lin <- as.vector(X_full %*% gamma + beta * g + off)
  mu      <- exp(eta_lin)

  u1 <- y - (y + r) * mu / (r + mu)
  u2 <- -(y + r) * r * mu / (r + mu)^2

  q <- -mu / (r + mu) + (y + r) * mu / (r + mu)^2

  a_r <- digamma(y + r) - digamma(r) +
    log(r) + 1 - log(r + mu) - (y + r) / (r + mu)

  b_r <- trigamma(y + r) - trigamma(r) +
    1 / r - 1 / (r + mu) + (y - mu) / (r + mu)^2

  p_gamma <- ncol(X_full)

  H <- matrix(0, nrow = p_gamma + 3, ncol = p_gamma + 3)

  idx_g  <- 1:p_gamma
  idx_b  <- p_gamma + 1
  idx_c  <- p_gamma + 2
  idx_ph <- p_gamma + 3

  d_eta_gamma <- X_full
  d_eta_b     <- beta * g
  d_eta_c     <- beta * gc

  # --- gamma blocks ---
  H[idx_g, idx_g] <- crossprod(d_eta_gamma, u2 * d_eta_gamma)

  H[idx_g, idx_b] <- colSums(u2 * d_eta_gamma * d_eta_b)
  H[idx_b, idx_g] <- H[idx_g, idx_b]

  H[idx_g, idx_c] <- colSums(u2 * d_eta_gamma * d_eta_c)
  H[idx_c, idx_g] <- H[idx_g, idx_c]

  H[idx_g, idx_ph] <- colSums((q * r) * d_eta_gamma)
  H[idx_ph, idx_g] <- H[idx_g, idx_ph]

  # --- b blocks ---
  H[idx_b, idx_b] <- sum(u2 * d_eta_b^2 + u1 * beta * g)

  H[idx_b, idx_c] <- sum(u2 * d_eta_b * d_eta_c + u1 * beta * gc)
  H[idx_c, idx_b] <- H[idx_b, idx_c]

  H[idx_b, idx_ph] <- sum((q * r) * d_eta_b)
  H[idx_ph, idx_b] <- H[idx_b, idx_ph]

  # --- c blocks ---
  d2_eta_cc       <- beta * gcc
  H[idx_c, idx_c] <- sum(u2 * d_eta_c^2 + u1 * d2_eta_cc)

  H[idx_c, idx_ph] <- sum((q * r) * d_eta_c)
  H[idx_ph, idx_c] <- H[idx_c, idx_ph]

  # --- phi block ---
  H[idx_ph, idx_ph] <- sum(b_r * r^2 + a_r * r)

  0.5 * (H + t(H))
}

# -------------------------------------------------------------------
# Wald and sandwich covariance matrices
# -------------------------------------------------------------------

#' Wald and sandwich covariance matrices for the structured threshold model
#'
#' Computes two covariance matrix estimates for the internal parameter
#' vector \eqn{(\gamma, b, c, \phi)} of the structured threshold model:
#' a model-based Wald estimator and a robust sandwich estimator.
#'
#' @details
#' **Wald estimator**: The observed information matrix is
#' \eqn{I = -H}, where \eqn{H} is the Hessian from
#' \code{hessian_nb_threshold_lag_softmax}. The Wald covariance is
#' \eqn{V_\text{wald} = I^{-1}}, computed via the eigendecomposition
#' of \eqn{I}. Eigenvalues below \code{eig_tol} are replaced by
#' \code{eig_tol} before inversion to ensure numerical stability when
#' the information matrix is near-singular, which can occur in short
#' time series.
#'
#' **Sandwich estimator**: The sandwich covariance is
#' \deqn{V_\text{sand} = V_\text{wald} \cdot J \cdot V_\text{wald}}
#' where \eqn{J = S^\top S} is the outer product of the score matrix
#' from \code{score_matrix_nb_threshold_lag_softmax}. This estimator
#' is consistent under mild departures from the negative binomial
#' model, including residual overdispersion and moderate temporal
#' dependence in the mortality counts.
#'
#' Both estimators are on the internal parameter scale
#' \eqn{(\gamma, b, c, \phi)}. Standard errors for \eqn{\beta} and
#' \eqn{r} on the original scale are obtained by the delta method
#' in \code{se_from_vcov_softmax}.
#'
#' @inheritParams score_matrix_nb_threshold_lag_softmax
#' @param eig_tol Positive numeric scalar. Minimum eigenvalue threshold
#'   for regularizing the observed information matrix before inversion.
#'   Default is \code{1e-8}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{H}{The raw Hessian matrix of dimension
#'     \eqn{(p_\gamma + 3) \times (p_\gamma + 3)}.}
#'   \item{eigvals}{Numeric vector. Eigenvalues of the observed
#'     information matrix \eqn{-H} before regularization.}
#'   \item{eigvals_reg}{Numeric vector. Eigenvalues after
#'     regularization (floored at \code{eig_tol}).}
#'   \item{V_wald}{Numeric matrix. Wald (model-based) covariance
#'     estimate.}
#'   \item{V_sandwich}{Numeric matrix. Sandwich (robust) covariance
#'     estimate.}
#' }
#'
#' @noRd
vcov_wald_sandwich_lag_softmax <- function(gamma, beta, c, phi, eta,
                                           X_full, temp_lag_mat, y, off,
                                           k = 5,
                                           eig_tol = 1e-8) {
  H <- hessian_nb_threshold_lag_softmax(
    gamma        = gamma,
    beta         = beta,
    c            = c,
    phi          = phi,
    eta          = eta,
    X_full       = X_full,
    temp_lag_mat = temp_lag_mat,
    y            = y,
    off          = off,
    k            = k
  )

  info <- -H
  eig  <- eigen(info, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  vals_reg <- pmax(vals, eig_tol)
  V_wald   <- vecs %*% diag(1 / vals_reg) %*% t(vecs)

  S      <- score_matrix_nb_threshold_lag_softmax(
    gamma        = gamma,
    beta         = beta,
    c            = c,
    phi          = phi,
    eta          = eta,
    X_full       = X_full,
    temp_lag_mat = temp_lag_mat,
    y            = y,
    off          = off,
    k            = k
  )

  J      <- crossprod(S)
  V_sand <- V_wald %*% J %*% V_wald

  list(
    H           = H,
    eigvals     = vals,
    eigvals_reg = vals_reg,
    V_wald      = V_wald,
    V_sandwich  = V_sand
  )
}

# -------------------------------------------------------------------
# Compatibility stub
# -------------------------------------------------------------------

#' Delta method transformation of covariance matrix for lag weights (stub)
#'
#' Returns the covariance matrix unchanged along with an empty Jacobian.
#' Retained for interface compatibility with the softmax lag-weight
#' parameterization, in which a delta method step would transform the
#' covariance from the \eqn{\eta} (softmax) scale to the \eqn{w}
#' (weight) scale. In the single-lag model no such transformation is
#' needed.
#'
#' @param V Numeric matrix. Covariance matrix on the internal parameter
#'   scale.
#' @param p_gamma Integer. Number of spline coefficients.
#' @param p_eta Integer. Number of softmax parameters (zero in the
#'   single-lag model).
#' @param eta Ignored.
#'
#' @return A named list with components:
#' \describe{
#'   \item{J_w_eta}{A \eqn{1 \times 0} matrix. Empty Jacobian.}
#'   \item{V_transformed}{The input \code{V} unchanged.}
#' }
#'
#' @noRd
transform_vcov_eta_to_w_softmax <- function(V, p_gamma, p_eta, eta) {
  list(
    J_w_eta       = matrix(0, nrow = 1, ncol = 0),
    V_transformed = V
  )
}

# -------------------------------------------------------------------
# Standard errors
# -------------------------------------------------------------------

#' Extract standard errors from a covariance matrix
#'
#' Extracts standard errors for all parameters of the structured
#' threshold model from a covariance matrix on the internal parameter
#' scale \eqn{(\gamma, b, c, \phi)}.
#'
#' @details
#' The standard error for \eqn{\beta = \exp(b)} is obtained by the
#' delta method:
#' \deqn{\widehat{\mathrm{se}}(\hat\beta) =
#'   \hat\beta \cdot \widehat{\mathrm{se}}(\hat b)}
#'
#' The standard errors for \eqn{\gamma}, \eqn{c}, and \eqn{\phi}
#' are read directly from the square root of the diagonal of \code{V}.
#'
#' This function can be applied to either the Wald or sandwich
#' covariance matrix. The resulting standard errors can then be used
#' to construct Wald confidence intervals for the HOT:
#' \deqn{\hat c \pm z_{0.975} \cdot \widehat{\mathrm{se}}(\hat c)}
#'
#' @param V Numeric matrix of dimension \eqn{(p_\gamma + 3) \times
#'   (p_\gamma + 3)}. Covariance matrix on the internal parameter scale,
#'   as returned by \code{vcov_wald_sandwich_lag_softmax}.
#' @param p_gamma Integer. Number of columns in the full design matrix
#'   \code{X_full}, equal to the number of spline columns plus one
#'   for the intercept.
#' @param p_eta Integer. Number of softmax lag-weight parameters.
#'   Must be \code{0} in the single-lag model.
#' @param beta_hat Positive numeric scalar. Estimated heat effect
#'   \eqn{\hat\beta}, required for the delta method transformation.
#'
#' @return Named numeric vector of length \eqn{p_\gamma + 3} with
#'   entries:
#' \describe{
#'   \item{se_gamma1, ..., se_gammaP}{Standard errors for the spline
#'     coefficients.}
#'   \item{se_beta}{Delta method standard error for \eqn{\hat\beta}.}
#'   \item{se_c}{Standard error for the HOT \eqn{\hat c}.}
#'   \item{se_phi}{Standard error for \eqn{\hat\phi = \log\hat r}.}
#' }
#'
#' @noRd
se_from_vcov_softmax <- function(V, p_gamma, p_eta, beta_hat = NULL) {
  idx_gamma <- 1:p_gamma
  idx_b     <- p_gamma + 1
  idx_c     <- p_gamma + 2
  idx_phi   <- p_gamma + 3

  se_gamma <- sqrt(diag(V)[idx_gamma])
  se_b     <- sqrt(V[idx_b, idx_b])
  se_c     <- sqrt(V[idx_c, idx_c])
  se_phi   <- sqrt(V[idx_phi, idx_phi])

  if (is.null(beta_hat)) stop("beta_hat must be supplied")
  se_beta <- as.numeric(beta_hat * se_b)

  out <- c(se_gamma, se_beta, se_c, se_phi)
  names(out) <- c(
    paste0("se_gamma", seq_along(idx_gamma)),
    "se_beta", "se_c", "se_phi"
  )
  out
}

#' Standard errors from a covariance matrix (lag-weight interface)
#'
#' A thin wrapper around \code{se_from_vcov_softmax} that accepts
#' the lag count \code{L} instead of \code{p_eta}. Provided for
#' consistency with the calling convention used in the one-stop
#' wrapper \code{fit_hot}.
#'
#' @param V Numeric matrix. Covariance matrix on the internal
#'   parameter scale.
#' @param p_gamma Integer. Number of columns in the full design
#'   matrix.
#' @param L Non-negative integer. Lag used in the model. In the
#'   single-lag model this implies \code{p_eta = 0}.
#' @param beta_hat Positive numeric scalar. Estimated heat effect,
#'   passed through to \code{se_from_vcov_softmax}.
#'
#' @return Named numeric vector of standard errors as returned by
#'   \code{se_from_vcov_softmax}.
#'
#' @noRd
se_from_vcov_weights_softmax <- function(V, p_gamma, L, beta_hat = NULL) {
  se_from_vcov_softmax(V, p_gamma, 0, beta_hat = beta_hat)
}