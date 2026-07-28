# Everything CopRegBAYES() rejects, and why. Some of these messages come from
# endogCopula's shared model-frame builder rather than from this package.
# The cdf / ties / nboots refusals live in test-arguments.R.

fit_small <- function(formula, data, ...)
  CopRegBAYES(formula, data = data, iterations = 20, burnin = 5, thin = 1,
              verbose = FALSE, ...)

test_that("the MCMC settings are validated", {
  dat <- make_sim_data(n = 40)$data

  expect_error(
    CopRegBAYES(y ~ z | x, data = dat, iterations = 100, burnin = 100,
                verbose = FALSE),
    "'burnin' must be smaller than 'iterations'"
  )
  expect_error(
    CopRegBAYES(y ~ z | x, data = dat, iterations = 100, burnin = 10,
                thin = 0, verbose = FALSE),
    "'thin' must be at least 1"
  )
})

test_that("a variable missing from the data is reported", {
  dat <- make_sim_data(n = 40)$data

  expect_error(fit_small(y ~ z | missing_var, dat),
               "object 'missing_var' not found")
  expect_error(fit_small(y ~ missing_var | x, dat),
               "object 'missing_var' not found")
})

test_that("the formula shape is validated", {
  dat <- make_sim_data(n = 40)$data

  expect_error(fit_small(y | x ~ z | x, dat),
               "exactly one dependent variable")
  expect_error(fit_small(y ~ z | x | z | x, dat),
               "at most three parts")
  expect_error(fit_small(y ~ 1 | x, dat),
               "No endogenous regressor found")
})

test_that("a categorical variable cannot be endogenous", {
  dat <- make_sim_data(n = 40)$data
  dat$g <- factor(rep(c("a", "b"), length.out = nrow(dat)))

  # a factor expands into contrast columns, so it has no single copula term
  expect_error(fit_small(y ~ g | x, dat), "needs its own copula term")
  # behind the "|" the same variable is fine
  expect_s3_class(fit_small(y ~ z | x + g, dat), "copregbayes")
})

test_that("a non-numeric response is rejected", {
  dat <- make_sim_data(n = 40)$data
  dat$g <- factor(rep(c("a", "b"), length.out = nrow(dat)))

  expect_error(fit_small(g ~ z | x, dat), "dependent variable must be numeric")
})

test_that("a constant regressor cannot enter the copula", {
  dat <- make_sim_data(n = 40)$data
  dat$k <- 1

  expect_error(fit_small(y ~ k | x, dat),
               "Endogenous regressor 'k' is constant")
  expect_error(fit_small(y ~ z | k, dat),
               "single unique value cannot enter the copula")
})

test_that("a rank-deficient design is rejected rather than silently fitted", {
  dat <- make_sim_data(n = 40)$data
  dat$x2 <- dat$x

  expect_error(fit_small(y ~ z | x + x2, dat), "rank deficient")
})

test_that("too few complete observations are rejected", {
  dat <- make_sim_data(n = 40)$data

  expect_error(fit_small(y ~ z | x, dat[1:3, ]),
               "Not enough complete observations")

  dat$y <- NA_real_
  expect_error(suppressMessages(fit_small(y ~ z | x, dat)),
               "Not enough complete observations")
})

test_that("unknown prior and starting values are named back", {
  dat <- make_sim_data(n = 40)$data

  expect_error(fit_small(y ~ z | x, dat, prior.args = list(nope = 1)),
               "Unknown entries in 'prior.args': nope")
  expect_error(fit_small(y ~ z | x, dat, start = list(nope = 1)),
               "Unknown entries in 'start': nope")
})

test_that("prior.args and start override the defaults", {
  fit <- bayes_fit(prior.args = list(a = 0.01, b = 0.01),
                   start = list(coefficients = c(0, 0, 0), sigma2 = 2))
  expect_equal(fit$prior.args$a, 0.01)
  expect_equal(fit$prior.args$b, 0.01)
  expect_true(all(is.finite(fit$draws)))

  # an explicit NULL leaves the default in place rather than deleting a
  # hyperparameter the sampler needs
  expect_equal(bayes_fit(prior.args = list(a = NULL))$prior.args$a, 0.001)
})
