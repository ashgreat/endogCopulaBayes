# validity() carries the identification checks and the convergence
# diagnostics in one object, plus the plot methods for looking at the chain
# itself.

test_that("validity() reports identification and convergence together", {
  fit <- bayes_fit(iterations = 600, burnin = 100)
  v <- validity(fit)

  expect_s3_class(v, "copregbayes.validity")
  expect_equal(v$endo.names, "z")
  expect_equal(v$ndraws, fit$ndraws)

  # [1] nonnormality of the endogenous regressors, which is what identifies
  # the correction; z is chi-squared, so it is skewed
  expect_equal(rownames(v$nonnormality), "z")
  expect_true(all(c("skewness", "ex.kurtosis") %in% colnames(v$nonnormality)))
  expect_gt(v$nonnormality[["z", "skewness"]], 0.5)

  # [2] posterior of the endogeneity correlations
  expect_equal(rownames(v$endogeneity), fit$rho.names)
  expect_equal(colnames(v$endogeneity),
               c("P. Mean", "2.5%", "97.5%", "P(rho > 0)"))
  expect_true(all(v$endogeneity[, "P(rho > 0)"] >= 0 &
                    v$endogeneity[, "P(rho > 0)"] <= 1))

  # [3] convergence, one row per retained parameter
  expect_s3_class(v$convergence, "data.frame")
  expect_equal(rownames(v$convergence), colnames(fit$draws))
  expect_equal(colnames(v$convergence), c("Geweke", "ESS", "AC(1)"))
  expect_true(all(is.finite(as.matrix(v$convergence))))
  expect_true(all(v$convergence$ESS > 0))
  expect_true(all(abs(v$convergence$`AC(1)`) <= 1))

  # Gelman-Rubin is skipped unless asked for, since it costs extra chains
  expect_null(v$gelman.rubin)
  expect_output(print(v), "Nonnormality of the endogenous regressors")
  expect_output(print(v), "Convergence of the chain")
  expect_output(print(v), "Gelman-Rubin was not computed")
})

test_that("validity(chains = ) runs extra chains from dispersed starts", {
  dat <- make_sim_data(n = 80)$data
  set.seed(6)
  fit <- CopRegBAYES(y ~ z | x, data = dat, iterations = 400, burnin = 100,
                     thin = 2, verbose = FALSE)

  v <- validity(fit, chains = 2, verbose = FALSE)

  expect_equal(names(v$gelman.rubin), colnames(fit$draws))
  expect_true(all(is.finite(v$gelman.rubin)))
  expect_true(all(v$gelman.rubin > 0))
  expect_output(print(v), "Gelman-Rubin from dispersed starts")

  expect_error(validity(fit, chains = 1), "at least two chains")
})

test_that("plot() draws every type without error", {
  fit <- bayes_fit()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  for (ty in c("trace", "acf", "pacf", "density")) {
    expect_no_error(plot(fit, type = ty))
  }
  # the default is the endogenous regressors, the ones the correction is about
  expect_no_error(plot(fit, which = c("z", "sigma2"), ask = FALSE))
  expect_no_error(plot(fit, which = 1))

  # type = "cdf" plots the estimated marginal distribution instead, which is
  # matched against the margins rather than against the parameter names
  expect_no_error(plot(fit, type = "cdf"))
  expect_no_error(plot(fit, which = c("z", "x"), type = "cdf", ask = FALSE))

  expect_error(plot(fit, which = "nope"), "Not a parameter of this model")
  expect_error(plot(fit, which = "nope", type = "cdf"),
               "No margin is estimated for")
  expect_error(plot(fit, type = "nope"), "'arg' should be one of")
})

# The convergence statistics are computed from the chain itself rather than
# taken from coda, so they are checked against series whose answers are known.

test_that("the effective sample size tracks the autocorrelation", {
  n <- 4000
  set.seed(1)
  iid <- stats::rnorm(n)
  # an AR(1) with phi = 0.8 has theoretical ESS/n = (1 - phi) / (1 + phi) =
  # 1/9; the initial positive sequence estimator computed here sits slightly
  # above that, since it sums a finite, truncated run of sample
  # autocorrelations rather than the infinite theoretical series
  ar1 <- as.numeric(stats::filter(stats::rnorm(n, sd = sqrt(1 - 0.8^2)),
                                  0.8, method = "recursive"))

  expect_gt(.bayes_ess(iid), 0.5 * n)
  expect_lt(abs(.bayes_ess(ar1) / n - 1 / 9), 0.04)
  expect_lt(.bayes_ess(ar1), .bayes_ess(iid))

  # truncating at the first negative autocorrelation bounds the sum below by
  # zero, so no series can be reported as carrying more effective draws than
  # it has draws -- including one that anticorrelates from lag one onwards
  set.seed(2)
  alternating <- rep(c(-1, 1), n / 2) + stats::rnorm(n, sd = 0.1)
  expect_lte(.bayes_ess(iid), n)
  expect_lte(.bayes_ess(ar1), n)
  expect_lte(.bayes_ess(alternating), n)

  # a chain that never moves has non-finite autocorrelations, and is
  # reported as not available rather than as NaN
  expect_true(is.na(.bayes_ess(rep(1, 500))))
})

test_that("Geweke's statistic flags a chain that has not settled", {
  set.seed(1)
  settled <- stats::rnorm(2000)
  expect_lt(abs(.bayes_geweke(settled)), 3)

  # the first tenth of the chain sits three standard deviations away
  drifting <- settled
  drifting[1:200] <- drifting[1:200] + 3
  expect_gt(abs(.bayes_geweke(drifting)), 5)

  # a chain stuck at a single value has no spectrum to divide by
  expect_true(is.na(.bayes_geweke(rep(1, 500))))
})
