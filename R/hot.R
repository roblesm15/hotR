# hot.R
# User-facing interface for the structured threshold model (STM).
#
# This file provides two functions:
#
#   fit_hot_internal()  internal wrapper that prepares data, constructs
#                       the spline design matrix, and calls the block
#                       coordinate descent engine in fit.R.
#
#   fit_hot()           the main exported function. Calls fit_hot_internal,
#                       computes Wald and sandwich covariance matrices,
#                       extracts standard errors, and assembles the full
#                       model object returned to the user.
#
# The target parameter of primary interest is the heat-onset risk
# threshold (HOT), returned as c_hat with Wald confidence interval
# ci_c_wald. Standard errors for all parameters are available on both
# the Wald and sandwich scales.
#
# Expected input format for city_data:
#   A data frame with columns:
#     date        : Date vector, one row per day
#     temperature : Numeric. Daily mean (or max) temperature in degrees C
#     deaths      : Non-negative integer. Daily all-cause mortality counts
#     population  : Positive numeric. Population at risk (used as offset)

# -------------------------------------------------------------------
# Internal fitting wrapper
# -------------------------------------------------------------------

#' Prepare data and fit the structured threshold model
#'
#' Internal wrapper that extracts variables from \code{city_data},
#' constructs the natural spline design matrix for the long-term
#' trend and seasonality component, builds the lagged temperature
#' exposure vector, and calls the block coordinate descent engine
#' \code{joint_nb_threshold_avg_lag_softmax}.
#'
#' @details
#' The number of spline degrees of freedom is computed as:
#' \deqn{df = \max\left(4,\, \text{round}(df\_per\_year \times
#'   \text{years})\right)}
#' where \emph{years} is the length of the series in years. The floor
#' of 4 ensures a minimum of flexibility even for very short series.
#' The spline is fit on the numeric representation of \code{date}
#' (i.e., days since the R epoch) and stored in the returned object
#' so that \code{predict} can be called on new data.
#'
#' @param city_data Data frame with columns \code{date} (Date),
#'   \code{temperature} (numeric), \code{deaths} (non-negative integer),
#'   and \code{population} (positive numeric).
#' @param df_per_year_spline Positive numeric scalar. Degrees of freedom
#'   per year allocated to the natural spline for long-term trend and
#'   seasonality. Typical values are 2--4. Default in \code{fit_hot}
#'   is \code{3}.
#' @param L Non-negative integer. Lag in days. The model links mortality
#'   at time \eqn{t} to temperature at time \eqn{t - L}.
#' @param k Numeric scalar. Smoothness parameter for the hinge
#'   transition. Default is \code{5}.
#' @param lower_q Numeric scalar in \eqn{(0, 1)}. Lower quantile of
#'   the lagged temperature exposure defining the left boundary of the
#'   search interval for the HOT. Default is \code{0.10}.
#' @param upper_q Numeric scalar in \eqn{(0, 1)}. Upper quantile of
#'   the lagged temperature exposure defining the right boundary of
#'   the search interval for the HOT. Default is \code{0.99}.
#' @param c0 Optional numeric scalar. Initial value for the threshold.
#'   If \code{NULL} (default), the median of the lagged temperature
#'   exposure is used.
#' @param eta0 Ignored. Included for interface compatibility.
#' @param tol Positive numeric scalar. Convergence tolerance on the
#'   absolute change in log-likelihood. Default is \code{1e-6}.
#' @param tol_c Positive numeric scalar. Convergence tolerance on the
#'   absolute change in the threshold. Default is \code{1e-5}.
#' @param tol_eta Ignored. Included for interface compatibility.
#' @param maxit Positive integer. Maximum number of block coordinate
#'   descent iterations. Default is \code{200}.
#' @param theta_update_every Positive integer. Number of iterations
#'   between dispersion updates. Default is \code{10}.
#' @param eta_maxit Ignored. Included for interface compatibility.
#' @param eta_reltol Ignored. Included for interface compatibility.
#' @param verbose Logical. If \code{TRUE}, prints iteration-level
#'   diagnostics. Default is \code{FALSE}.
#' @param return_data_mat Logical. If \code{TRUE}, the design matrix,
#'   lag matrix, response vector, offset vector, and spline object are
#'   retained in the returned list. If \code{FALSE} (default), these
#'   are set to \code{NULL} to reduce memory usage.
#'
#' @return A named list with components:
#' \describe{
#'   \item{c_hat}{Numeric scalar. Estimated HOT in degrees C.}
#'   \item{theta_hat}{Numeric scalar. Estimated negative binomial
#'     dispersion parameter \eqn{\hat r}.}
#'   \item{coef}{Numeric vector. Estimated coefficients
#'     \eqn{(\hat\gamma_1, \ldots, \hat\gamma_{p_\gamma}, \hat\beta)}.
#'     Note \eqn{\hat\beta} is on the original scale, not the log
#'     scale.}
#'   \item{b_hat}{Numeric scalar. \eqn{\hat b = \log\hat\beta}.}
#'   \item{lag_weights}{Numeric scalar 1. Retained for compatibility.}
#'   \item{z}{Numeric vector. Lagged temperature exposure used in
#'     fitting (if \code{return_data_mat = TRUE}, else \code{NULL}).}
#'   \item{df_loc}{Integer. Degrees of freedom used for the spline.}
#'   \item{X_spline}{Spline design matrix without intercept (if
#'     \code{return_data_mat = TRUE}, else \code{NULL}).}
#'   \item{temp_lag_mat}{Lagged temperature vector (if
#'     \code{return_data_mat = TRUE}, else \code{NULL}).}
#'   \item{ns_time_obj}{The fitted \code{splines::ns} object, used
#'     for prediction on new data (if \code{return_data_mat = TRUE},
#'     else \code{NULL}).}
#'   \item{y}{Response vector (if \code{return_data_mat = TRUE},
#'     else \code{NULL}).}
#'   \item{off}{Offset vector (if \code{return_data_mat = TRUE},
#'     else \code{NULL}).}
#'   \item{k}{Numeric scalar. Smoothness parameter used.}
#'   \item{L}{Integer. Lag used.}
#'   \item{converged}{Logical. Whether the algorithm converged.}
#'   \item{logLik}{Numeric scalar. Log-likelihood at convergence.}
#' }
#'
#' @noRd
fit_hot_internal <- function(
    city_data,
    df_per_year_spline,
    L,
    k = 5,
    lower_q = 0.10,
    upper_q = 0.99,
    c0 = NULL,
    eta0 = NULL,
    tol = 1e-6,
    tol_c = 1e-5,
    tol_eta = 1e-5,
    maxit = 200,
    theta_update_every = 10,
    eta_maxit = 10,
    eta_reltol = 1e-4,
    verbose = FALSE,
    return_data_mat = FALSE
) {
  if (L < 0) stop("L must be at least 0.")

  tnum <- as.numeric(city_data$date)
  temp <- city_data$temperature
  y    <- city_data$deaths
  off  <- log(city_data$population)

  no_years <- round(
    (max(city_data$date) - min(city_data$date) + 1) / 365.25
  )
  df_loc <- max(4, round(df_per_year_spline * no_years))

  ns_time_obj <- splines::ns(tnum, df = df_loc)
  X_spline    <- predict(ns_time_obj, newx = tnum)

  temp_lag_mat <- extract_single_lag(temp, L)

  fit <- joint_nb_threshold_avg_lag_softmax(
    y                  = y,
    off                = off,
    X_spline           = X_spline,
    temp_lag_mat       = temp_lag_mat,
    c0                 = c0,
    eta0               = eta0,
    k                  = k,
    lower_q            = lower_q,
    upper_q            = upper_q,
    tol                = tol,
    tol_c              = tol_c,
    tol_eta            = tol_eta,
    maxit              = maxit,
    theta_update_every = theta_update_every,
    eta_maxit          = eta_maxit,
    eta_reltol         = eta_reltol,
    verbose            = verbose
  )

  list(
    c_hat           = fit$c,
    theta_hat       = fit$r,
    coef            = c(fit$gamma, g = fit$beta),
    b_hat           = fit$b,
    lag_weight_par  = numeric(0),
    lag_weights     = 1,
    z               = if (return_data_mat) fit$z         else NULL,
    df_loc          = df_loc,
    X_spline        = if (return_data_mat) X_spline      else NULL,
    temp_lag_mat    = if (return_data_mat) temp_lag_mat  else NULL,
    ns_time_obj     = if (return_data_mat) ns_time_obj   else NULL,
    y               = if (return_data_mat) y             else NULL,
    off             = if (return_data_mat) off           else NULL,
    k               = k,
    L               = L,
    converged       = fit$converged,
    logLik          = fit$logLik,
    # retained internally for vcov step regardless of return_data_mat
    .y              = y,
    .off            = off,
    .X_spline       = X_spline,
    .temp_lag_mat   = temp_lag_mat
  )
}

# -------------------------------------------------------------------
# Main exported function
# -------------------------------------------------------------------

#' Fit the structured threshold model for heat-onset risk estimation
#'
#' Estimates the heat-onset risk threshold (HOT) from a daily mortality
#' time series using the structured threshold model (STM). The HOT is
#' the temperature above which heat-attributable mortality begins to
#' rise, estimated jointly with long-term trend, seasonality, and
#' overdispersion parameters.
#'
#' @details
#' ## Model
#' Mortality counts \eqn{y_t} are modeled as:
#' \deqn{y_t \sim \mathrm{NegBin}(\mu_t,\, r)}
#' \deqn{\log \mu_t = \mathbf{x}_t^\top \gamma +
#'   \beta \cdot g(T_{t-L};\, c, k) + \log(\mathrm{pop}_t)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\mathbf{x}_t} is a natural spline basis vector for
#'     time, capturing long-term trend and seasonality with
#'     \code{df_per_year_spline} degrees of freedom per year.
#'   \item \eqn{g(T_{t-L}; c, k)} is the smooth hinge function
#'     evaluated at the temperature lagged by \eqn{L} days. It is
#'     approximately zero for \eqn{T_{t-L} \leq c} and grows
#'     quadratically above \eqn{c}.
#'   \item \eqn{c} is the HOT, the primary parameter of interest.
#'   \item \eqn{\beta > 0} is the heat effect coefficient, constrained
#'     positive by the parameterization \eqn{\beta = \exp(b)}.
#'   \item \eqn{r} is the negative binomial dispersion (size) parameter.
#' }
#'
#' ## Estimation
#' Parameters are estimated by block coordinate descent, alternating
#' between a BFGS step for \eqn{(\gamma, b)}, a golden-section search
#' for \eqn{c}, and periodic profile ML updates for \eqn{r}. See
#' \code{joint_nb_threshold_avg_lag_softmax} for full algorithmic
#' details.
#'
#' ## Inference
#' Two covariance estimators are returned:
#' \itemize{
#'   \item **Wald**: based on the observed information matrix
#'     \eqn{-H(\hat\theta)}. Valid under correct negative binomial
#'     specification.
#'   \item **Sandwich**: robust to overdispersion misspecification and
#'     moderate temporal dependence in mortality counts.
#' }
#' Wald confidence intervals for the HOT are provided directly as
#' \code{ci_c_wald}. For short or noisy series the sandwich interval
#' is recommended.
#'
#' ## Data requirements
#' \code{city_data} must contain one row per day with no gaps in the
#' date sequence. Missing values in \code{temperature} or \code{deaths}
#' are not currently supported. The minimum recommended series length
#' is three years.
#'
#' @param city_data Data frame with columns:
#'   \describe{
#'     \item{date}{Date vector ordered chronologically with no gaps.}
#'     \item{temperature}{Numeric. Daily mean or maximum temperature
#'       in degrees C.}
#'     \item{deaths}{Non-negative integer. Daily all-cause or
#'       cause-specific mortality counts.}
#'     \item{population}{Positive numeric. Population at risk,
#'       used as a multiplicative offset on the log scale.}
#'   }
#' @param df_per_year_spline Positive numeric scalar. Degrees of freedom
#'   per year for the natural spline capturing long-term trend and
#'   seasonality. Smaller values (e.g., \code{2}) impose smoother
#'   baselines; larger values (e.g., \code{4}) allow more flexibility
#'   but increase the risk of overfitting in short series. Default
#'   is \code{3}.
#' @param L Non-negative integer. Lag in days between temperature
#'   exposure and mortality outcome. \code{L = 0} uses same-day
#'   temperature; \code{L = 1} uses the previous day, and so on.
#'   The choice of \code{L} should be guided by prior knowledge or
#'   cross-validation. Default is \code{0}.
#' @param k Numeric scalar. Smoothness parameter for the hinge
#'   transition around the HOT. Larger values make the transition
#'   sharper; the default \code{k = 5} provides a smooth but
#'   well-defined bend. In most applications this does not require
#'   tuning.
#' @param lower_q Numeric scalar in \eqn{(0, 1)}. Lower quantile of
#'   the lagged temperature distribution defining the left boundary
#'   of the search interval for the HOT. Temperatures below this
#'   quantile are excluded from consideration as threshold candidates.
#'   Default is \code{0.10}.
#' @param upper_q Numeric scalar in \eqn{(0, 1)}. Upper quantile
#'   defining the right boundary of the search interval. Default is
#'   \code{0.99}.
#' @param c0 Optional numeric scalar. Initial value for the HOT in
#'   degrees C. If \code{NULL} (default), the median of the lagged
#'   temperature series is used.
#' @param tol Positive numeric scalar. Convergence tolerance on the
#'   absolute change in log-likelihood between iterations. Default
#'   is \code{1e-6}.
#' @param tol_c Positive numeric scalar. Convergence tolerance on the
#'   absolute change in the HOT between iterations. Default is
#'   \code{1e-5}.
#' @param maxit Positive integer. Maximum number of block coordinate
#'   descent iterations. Default is \code{200}.
#' @param theta_update_every Positive integer. Number of iterations
#'   between profile ML updates of the dispersion parameter. Default
#'   is \code{10}.
#' @param verbose Logical. If \code{TRUE}, prints iteration diagnostics
#'   including log-likelihood, current HOT estimate, heat effect, and
#'   dispersion. Useful for monitoring convergence. Default is
#'   \code{FALSE}.
#' @param return_data_mat Logical. If \code{TRUE}, the design matrices,
#'   lag vector, response, and offset are included in the returned
#'   object. Useful for post-hoc diagnostics and visualization but
#'   increases memory usage. Default is \code{FALSE}.
#' @param eig_tol Positive numeric scalar. Minimum eigenvalue threshold
#'   for regularizing the observed information matrix before inversion.
#'   Rarely needs adjustment. Default is \code{1e-8}.
#'
#' @return An object of class \code{"hot_fit"}: a named list with
#'   components:
#' \describe{
#'   \item{c_hat}{Numeric scalar. Estimated HOT in degrees C.}
#'   \item{ci_c_wald}{Named numeric vector of length 2. Wald 95\%
#'     confidence interval for the HOT: \code{c(lower, upper)}.}
#'   \item{theta_hat}{Numeric scalar. Estimated negative binomial
#'     dispersion parameter \eqn{\hat r}.}
#'   \item{coef}{Numeric vector of estimated coefficients
#'     \eqn{(\hat\gamma_1, \ldots, \hat\gamma_{p_\gamma}, \hat\beta)}.}
#'   \item{b_hat}{Numeric scalar. \eqn{\hat b = \log\hat\beta}.}
#'   \item{se_wald}{Named numeric vector. Wald standard errors for
#'     all parameters on the internal scale, with \code{se_beta}
#'     obtained by the delta method and \code{se_c} the standard
#'     error of the HOT.}
#'   \item{se_sandwich}{Named numeric vector. Sandwich standard errors,
#'     structured identically to \code{se_wald}.}
#'   \item{V_wald}{Numeric matrix. Full Wald covariance matrix on the
#'     internal parameter scale \eqn{(\gamma, b, c, \phi)}.}
#'   \item{V_sandwich}{Numeric matrix. Full sandwich covariance matrix.}
#'   \item{eigvals}{Numeric vector. Eigenvalues of the observed
#'     information matrix before regularization. Values close to zero
#'     indicate near-singular information, which may occur in short
#'     series or when the threshold is weakly identified.}
#'   \item{eigvals_reg}{Numeric vector. Eigenvalues after regularization.}
#'   \item{df_loc}{Integer. Degrees of freedom used for the spline.}
#'   \item{L}{Integer. Lag used.}
#'   \item{k}{Numeric scalar. Smoothness parameter used.}
#'   \item{logLik}{Numeric scalar. Log-likelihood at convergence.}
#'   \item{converged}{Logical. Whether the algorithm converged within
#'     \code{maxit} iterations.}
#'   \item{X_full}{Numeric matrix. Full design matrix including intercept
#'     and spline columns. Always returned (used internally for vcov).}
#'   \item{X_spline, temp_lag_mat, ns_time_obj, y, off, z}{Retained
#'     if \code{return_data_mat = TRUE}, otherwise \code{NULL}.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Minimal example using a single city data frame
#' data(ahmedabad)  # example dataset included in hotR
#'
#' fit <- fit_hot(
#'   city_data          = ahmedabad,
#'   df_per_year_spline = 3,
#'   L                  = 1,
#'   verbose            = FALSE
#' )
#'
#' # Estimated heat-onset risk threshold
#' fit$c_hat
#'
#' # 95% Wald confidence interval for the HOT
#' fit$ci_c_wald
#'
#' # Standard errors (sandwich recommended for short series)
#' fit$se_sandwich["se_c"]
#' }
fit_hot <- function(
    city_data,
    df_per_year_spline = 3,
    L = 0,
    k = 5,
    lower_q = 0.10,
    upper_q = 0.99,
    c0 = NULL,
    tol = 1e-6,
    tol_c = 1e-5,
    maxit = 200,
    theta_update_every = 10,
    verbose = FALSE,
    return_data_mat = FALSE,
    eig_tol = 1e-8
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
  if (L < 0 || L != round(L)) {
    stop("L must be a non-negative integer.")
  }

  # --- fit ---
  fit <- fit_hot_internal(
    city_data          = city_data,
    df_per_year_spline = df_per_year_spline,
    L                  = L,
    k                  = k,
    lower_q            = lower_q,
    upper_q            = upper_q,
    c0                 = c0,
    tol                = tol,
    tol_c              = tol_c,
    maxit              = maxit,
    theta_update_every = theta_update_every,
    verbose            = verbose,
    return_data_mat    = return_data_mat
  )

  # --- covariance ---
  X_full  <- cbind(1, fit$.X_spline)
  p_gamma <- ncol(X_full)
  beta_hat <- as.numeric(fit$coef[p_gamma + 1])

  vc <- vcov_wald_sandwich_lag_softmax(
    gamma        = fit$coef[1:p_gamma],
    beta         = beta_hat,
    c            = fit$c_hat,
    phi          = log(fit$theta_hat),
    eta          = numeric(0),
    X_full       = X_full,
    temp_lag_mat = fit$.temp_lag_mat,
    y            = fit$.y,
    off          = fit$.off,
    k            = fit$k,
    eig_tol      = eig_tol
  )

  se_wald <- se_from_vcov_softmax(
    vc$V_wald, p_gamma, 0, beta_hat = beta_hat
  )
  se_sand <- se_from_vcov_softmax(
    vc$V_sandwich, p_gamma, 0, beta_hat = beta_hat
  )

  # --- assemble output ---
  fit$X_full          <- X_full
  fit$V_wald          <- vc$V_wald
  fit$V_sandwich      <- vc$V_sandwich
  fit$se_wald         <- se_wald
  fit$se_sandwich     <- se_sand
  fit$eigvals         <- vc$eigvals
  fit$eigvals_reg     <- vc$eigvals_reg
  fit$ci_c_wald       <- c(
    lower = fit$c_hat - qnorm(0.975) * se_wald["se_c"],
    upper = fit$c_hat + qnorm(0.975) * se_wald["se_c"]
  )

  # remove internal slots
  fit$.y            <- NULL
  fit$.off          <- NULL
  fit$.X_spline     <- NULL
  fit$.temp_lag_mat <- NULL
  fit$lag_weight_par <- NULL

  class(fit) <- "hot_fit"
  fit
}