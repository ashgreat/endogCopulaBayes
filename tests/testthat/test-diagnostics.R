# The sampler calls several symbols unqualified (dCopula, rmvnorm,
# rinvwishart, ...); until the NAMESPACE gains import directives, the
# packages must be attached.
suppressMessages({
  library(copula)
  library(mvtnorm)
  library(LaplacesDemon)
  library(invgamma)
})

seeded_fit <- function() {
  set.seed(321)
  n <- 80
  z <- rlnorm(n)
  x <- rnorm(n)
  e <- rnorm(n, sd = sqrt(2))
  dat <- data.frame(y = 2 - 4 * z + 6 * x + e, z = z, x = x)
  CopRegBayes(dat, iterations = 400, burnin = 100, thin = 2, seed = 42)
}

test_that("diagnostics() matches coda::effectiveSize and coda::geweke.diag", {
  skip_if_not_installed("coda")

  fit <- seeded_fit()
  d <- diagnostics(fit)

  draws <- fit$posterior[, fit$parameters, drop = FALSE]
  mc <- coda::as.mcmc(draws)

  ref_ess <- coda::effectiveSize(mc)
  ref_geweke <- coda::geweke.diag(mc, frac1 = 0.1, frac2 = 0.5)$z

  for (p in fit$parameters) {
    our_ess <- d[p, "ess"]
    their_ess <- unname(ref_ess[p])
    if (their_ess == 0) {
      expect_equal(our_ess, 0)
    } else {
      expect_lt(abs(our_ess - their_ess) / their_ess, 0.10)
    }

    our_z <- d[p, "geweke_z"]
    their_z <- unname(ref_geweke[p])
    if (is.na(their_z)) {
      # coda returns NaN (0/0) for a degenerate (constant) column, e.g.
      # rho_xe, which the sampler forces to zero; we report NA for the same
      # case rather than propagating a NaN.
      expect_true(is.na(our_z))
    } else {
      expect_lt(abs(our_z - their_z), 0.15)
    }
  }
})

test_that("diagnostics() returns the expected shape and names", {
  fit <- seeded_fit()
  d <- diagnostics(fit)

  expect_s3_class(d, "data.frame")
  expect_equal(colnames(d), c("ess", "geweke_z", "n_draws"))
  expect_equal(rownames(d), fit$parameters)
  expect_equal(nrow(d), length(fit$parameters))
  expect_true(all(d$n_draws == nrow(fit$posterior)))
  expect_true(all(is.finite(d$ess)))
  expect_true(all(d$ess >= 0))
})

test_that("diagnostics() respects formula-interface variable labels", {
  set.seed(9)
  n <- 40
  dat <- data.frame(sales = rnorm(n), price = rlnorm(n), income = rnorm(n))
  fit <- CopRegBayes(sales ~ price | income, data = dat,
                     iterations = 60, burnin = 10, thin = 2, seed = 3)

  d <- diagnostics(fit)
  expect_true("price" %in% rownames(d))
  expect_true("income" %in% rownames(d))
  expect_false("beta_z" %in% rownames(d))
  expect_false("beta_x" %in% rownames(d))
})

test_that("plot.endog_copula_bayes() runs without error", {
  fit <- seeded_fit()

  pdf(NULL)
  on.exit(dev.off())

  expect_no_error(plot(fit))
  expect_no_error(plot(fit, params = c("beta_z", "beta_x")))
  expect_no_error(plot(fit, which = "trace"))
  expect_no_error(plot(fit, which = "density"))
  expect_error(plot(fit, params = "not_a_parameter"), "Unknown parameter")
})

test_that("degenerate constant-column diagnostics do not error", {
  fit <- seeded_fit()

  # Force a constant column to exercise the degenerate-spectrum guard.
  fit$posterior[, "rho_xe"] <- 0

  d <- expect_no_error(diagnostics(fit))
  expect_equal(d["rho_xe", "ess"], 0)
  expect_true(is.na(d["rho_xe", "geweke_z"]))

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(fit, params = "rho_xe"))
})
