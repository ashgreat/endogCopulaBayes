test_that("CopRegBAYES fits a well-formed posterior and corrects the endogeneity bias", {
  sim  <- make_sim_data()
  dat  <- sim$data
  true <- sim$true

  iterations <- 3000
  burnin     <- 500
  thin       <- 10

  fit <- CopRegBAYES(y ~ z | x, data = dat, iterations = iterations,
                      burnin = burnin, thin = thin, verbose = FALSE)

  expect_s3_class(fit, "copregbayes")

  # coef() is the posterior mean, finite and named
  cf <- coef(fit)
  expect_identical(cf, fit$coefficients)
  expect_true(is.numeric(cf))
  expect_true(all(is.finite(cf)))
  expect_named(cf)

  # the corrected coefficient on the endogenous regressor beats OLS
  ols <- stats::coef(stats::lm(y ~ z + x, data = dat))
  expect_lt(abs(cf[["z"]] - true[["z"]]), abs(ols[["z"]] - true[["z"]]))

  # draws matrices have the dimensions implied by iterations/burnin/thin
  ndraw <- length(seq.int(burnin + 1L, iterations, by = thin))
  expect_equal(fit$ndraws, ndraw)
  expect_equal(dim(fit$coefficient.draws), c(ndraw, length(cf)))
  expect_length(fit$sigma2.draws, ndraw)
  expect_equal(nrow(fit$Sigma.draws), ndraw)
  expect_equal(nrow(fit$draws), ndraw)

  # summary()
  s <- summary(fit)
  expect_s3_class(s, "summary.copregbayes")
  expect_equal(dim(s$coefficients), c(length(cf), 5L))
  expect_equal(rownames(s$coefficients), names(cf))

  # confint()
  ci <- confint(fit)
  expect_equal(nrow(ci), ncol(fit$draws))
  expect_equal(rownames(ci), colnames(fit$draws))
  expect_equal(colnames(ci), c("2.5%", "97.5%"))

  # vcov()
  vc <- vcov(fit)
  expect_equal(dim(vc), c(length(cf), length(cf)))

  # residuals() / fitted()
  res    <- residuals(fit)
  fitted <- fitted(fit)
  expect_length(res, nrow(dat))
  expect_length(fitted, nrow(dat))
  expect_equal(unname(res), unname(dat$y - fitted))

  # predict()
  pred <- predict(fit, newdata = dat[1:5, ])
  expect_length(pred, 5L)
  expect_equal(unname(pred), unname(fitted[1:5]))
  expect_equal(unname(predict(fit)), unname(fitted))

  # nobs()
  expect_equal(nobs(fit), nrow(dat))

  # validity()
  v <- validity(fit)
  expect_s3_class(v, "copregbayes.validity")
  expect_equal(nrow(v$nonnormality), length(fit$endo.names))
  expect_equal(nrow(v$endogeneity), length(fit$rho.names))
  expect_equal(nrow(v$convergence), ncol(fit$draws))
})
