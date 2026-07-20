test_that("smooth hinge derivatives agree with finite differences", {
  x <- c(20, 25, 30, 35)
  threshold <- 28
  k <- 5
  eps <- 1e-5

  first_numeric <- (
    hotR:::smooth_hinge(x, threshold + eps, k) -
      hotR:::smooth_hinge(x, threshold - eps, k)
  ) / (2 * eps)
  second_numeric <- (
    hotR:::dg_dc(x, threshold + eps, k) -
      hotR:::dg_dc(x, threshold - eps, k)
  ) / (2 * eps)

  expect_equal(hotR:::dg_dc(x, threshold, k), first_numeric,
               tolerance = 1e-6)
  expect_equal(hotR:::d2g_dc2(x, threshold, k), second_numeric,
               tolerance = 1e-5)
})

test_that("analytic likelihood scores agree with finite differences", {
  set.seed(2026)
  n <- 20
  X <- cbind(1, seq(-1, 1, length.out = n))
  temperature <- seq(22, 34, length.out = n)
  deaths <- stats::rnbinom(n, mu = 50, size = 25)
  offset <- rep(log(1e6), n)
  par <- c(-10, 0.02, log(0.01), 28, log(25))

  log_likelihood <- function(value) {
    gamma <- value[1:2]
    beta <- exp(value[3])
    eta <- drop(
      X %*% gamma +
        beta * hotR:::smooth_hinge(temperature, value[4], 5) +
        offset
    )
    hotR:::nb_loglik_full(eta, deaths, exp(value[5]))
  }

  eps <- 1e-5
  score_numeric <- vapply(seq_along(par), function(j) {
    upper <- lower <- par
    upper[j] <- upper[j] + eps
    lower[j] <- lower[j] - eps
    (log_likelihood(upper) - log_likelihood(lower)) / (2 * eps)
  }, numeric(1))

  score_analytic <- colSums(
    hotR:::score_matrix_nb_threshold_lag_softmax(
      gamma = par[1:2],
      beta = exp(par[3]),
      c = par[4],
      phi = par[5],
      eta = numeric(0),
      X_full = X,
      temp_lag_mat = temperature,
      y = deaths,
      off = offset,
      k = 5
    )
  )

  expect_equal(unname(score_analytic), score_numeric, tolerance = 1e-5)
})

test_that("lag matrices use the documented first-value padding", {
  temperature <- c(10, 11, 12, 13)

  expect_equal(
    hotR:::build_temp_lag_matrix(temperature, L = 2),
    cbind(
      c(10, 11, 12, 13),
      c(10, 10, 11, 12),
      c(10, 10, 10, 11)
    )
  )
  expect_equal(
    hotR:::extract_single_lag(temperature, L = 1),
    c(10, 10, 11, 12)
  )

  expect_error(
    hotR:::build_temp_lag_matrix(temperature, L = 4),
    "smaller than the number"
  )
  expect_error(
    hotR:::build_temp_lag_matrix(c(10, NA_real_), L = 0),
    "finite numeric vector"
  )
})
