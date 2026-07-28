# cdf, ties and nboots do not apply to CopRegBAYES: the marginal CDFs are
# parameters drawn by the sampler, and inference is posterior rather than
# bootstrap, so supplying any of them is an error.

test_that("cdf, ties and nboots are rejected", {
  dat <- make_sim_data(n = 60, seed = 2)$data

  expect_error(
    CopRegBAYES(y ~ z | x, data = dat, cdf = "empirical",
                iterations = 50, burnin = 10, thin = 5, verbose = FALSE),
    "marginal CDFs"
  )
  expect_error(
    CopRegBAYES(y ~ z | x, data = dat, ties = "average",
                iterations = 50, burnin = 10, thin = 5, verbose = FALSE),
    "marginal CDFs"
  )
  expect_error(
    CopRegBAYES(y ~ z | x, data = dat, nboots = 100,
                iterations = 50, burnin = 10, thin = 5, verbose = FALSE),
    "posterior rather than bootstrap"
  )
})

test_that("an endogenous regressor with tied values still fits", {
  # marginal masses sit on unique values, so ties (a handful of distinct
  # levels repeated many times) are the normal case, not a problem
  set.seed(11)
  n <- 150
  z <- sample(1:5, n, replace = TRUE)
  x <- stats::rnorm(n)
  y <- 1 + 2 * z + 0.5 * x + stats::rnorm(n)
  dat <- data.frame(y = y, z = z, x = x)

  fit <- suppressWarnings(
    CopRegBAYES(y ~ z | x, data = dat, iterations = 80, burnin = 20,
                thin = 5, verbose = FALSE)
  )
  expect_s3_class(fit, "copregbayes")
  expect_equal(ncol(fit$lambda.draws[["z"]]), length(unique(z)))
})

test_that("a binary endogenous regressor still fits", {
  set.seed(12)
  n <- 150
  z <- rbinom(n, 1, 0.5)
  x <- stats::rnorm(n)
  y <- 1 + 2 * z + 0.5 * x + stats::rnorm(n)
  dat <- data.frame(y = y, z = z, x = x)

  fit <- suppressWarnings(
    CopRegBAYES(y ~ z | x, data = dat, iterations = 80, burnin = 20,
                thin = 5, verbose = FALSE)
  )
  expect_s3_class(fit, "copregbayes")
  expect_equal(ncol(fit$lambda.draws[["z"]]), 2L)
})
