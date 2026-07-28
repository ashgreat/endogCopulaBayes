# The two-part formula y ~ endogenous | exogenous is the whole interface:
# position around the "|" decides what gets a copula term, and every
# non-intercept column of the design matrix enters the copula.

test_that("the fit records what sits on each side of the '|'", {
  fit <- bayes_fit()

  expect_equal(fit$endo.names, "z")
  expect_equal(fit$exog.names, "x")
  expect_equal(fit$rho.names, "rho(z*, xi*)")
  expect_equal(colnames(fit$draws),
               c("(Intercept)", "z", "x", "sigma2",
                 "rho(z*, x*)", "rho(z*, xi*)"))
  expect_named(fit$lambda.draws, c("z", "x"))
  expect_equal(deparse(formula(fit)), "y ~ z | x")
})

test_that("exogeneity is imposed exactly, not approximately", {
  fit <- bayes_fit()

  # rho(x*, xi*) is held at zero by the block-recursive parameterisation
  # rather than merely shrunk towards it, so this is == and not a tolerance.
  expect_equal(colnames(fit$Sigma.draws)[!fit$Sigma.free], "rho(x*, xi*)")
  expect_true(all(fit$Sigma.draws[, !fit$Sigma.free] == 0))

  # the free entries are the ones carried into $draws
  expect_true(all(colnames(fit$Sigma.draws)[fit$Sigma.free] %in%
                    colnames(fit$draws)))
  expect_false(any(colnames(fit$Sigma.draws)[!fit$Sigma.free] %in%
                     colnames(fit$draws)))
})

test_that("several endogenous and exogenous regressors are all carried through", {
  dat <- make_sim_data()$data
  set.seed(3)
  dat$z2 <- exp(stats::rnorm(nrow(dat)))
  dat$w  <- stats::rnorm(nrow(dat))

  fit <- bayes_fit(y ~ z + z2 | x + w, data = dat)

  expect_equal(names(coef(fit)), c("(Intercept)", "z", "z2", "x", "w"))
  expect_equal(fit$endo.names, c("z", "z2"))
  expect_equal(fit$exog.names, c("x", "w"))
  expect_equal(fit$rho.names, c("rho(z*, xi*)", "rho(z2*, xi*)"))
  expect_named(fit$lambda.draws, c("z", "z2", "x", "w"))

  # one correlation per pair of the four regressors and the error, with both
  # exogenous-with-error entries restricted to zero
  expect_equal(ncol(fit$Sigma.draws), choose(5, 2))
  expect_equal(colnames(fit$Sigma.draws)[!fit$Sigma.free],
               c("rho(x*, xi*)", "rho(w*, xi*)"))
  expect_true(all(fit$Sigma.draws[, !fit$Sigma.free] == 0))
})

test_that("a factor behind the '|' expands into contrast columns that enter the copula", {
  dat <- make_sim_data()$data
  set.seed(4)
  dat$g <- factor(sample(c("a", "b", "c"), nrow(dat), replace = TRUE))

  fit <- bayes_fit(y ~ z | x + g, data = dat)

  expect_equal(names(coef(fit)), c("(Intercept)", "z", "x", "gb", "gc"))
  expect_equal(fit$exog.names, c("x", "gb", "gc"))
  # each dummy is its own margin, with two unique values
  expect_equal(vapply(fit$lambda.draws[c("gb", "gc")], ncol, integer(1)),
               c(gb = 2L, gc = 2L))
  expect_equal(fit$xlevels, list(g = c("a", "b", "c")))
})

test_that("a transformed regressor gets its own copula term", {
  dat <- make_sim_data()$data
  fit <- bayes_fit(y ~ I(z^2) | x, data = dat)

  expect_equal(fit$endo.names, "I(z^2)")
  expect_equal(fit$rho.names, "rho(I(z^2)*, xi*)")
  expect_true("rho(I(z^2)*, xi*)" %in% colnames(fit$draws))
  expect_equal(fit$margins[["I(z^2)"]]$uv, sort(unique(dat$z^2)))
})

test_that("a one-part formula fits with no exogenous regressors", {
  fit <- bayes_fit(y ~ z)

  expect_equal(colnames(fit$draws),
               c("(Intercept)", "z", "sigma2", "rho(z*, xi*)"))
  expect_equal(fit$exog.names, character(0))
  expect_true(all(fit$Sigma.free))
  expect_true(all(is.finite(fit$draws)))
})

test_that("intercept removal is honoured", {
  fit <- bayes_fit(y ~ z | x - 1)

  expect_equal(names(coef(fit)), c("z", "x"))
  expect_false("(Intercept)" %in% colnames(fit$draws))
  expect_equal(fit$exog.names, "x")
})

test_that("subset selects the rows used", {
  dat <- make_sim_data()$data
  fit <- bayes_fit(data = dat, subset = seq_len(60))

  expect_equal(nobs(fit), 60)
  expect_equal(nrow(fit$X), 60)
})

test_that("incomplete rows are dropped, with a message", {
  dat <- make_sim_data()$data
  dat$x[1:3] <- NA

  set.seed(1)
  expect_message(
    fit <- CopRegBAYES(y ~ z | x, data = dat, iterations = 50, burnin = 10,
                       thin = 1, verbose = FALSE),
    "removed because of missing values"
  )
  expect_equal(nobs(fit), nrow(dat) - 3L)
  expect_true(all(is.finite(fit$draws)))
})
