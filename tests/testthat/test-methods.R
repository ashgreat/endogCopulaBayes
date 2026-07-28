# The extractor and summary methods. test-bayes.R checks that they all run on
# a real fit; this file checks that what they return is what the posterior
# draws say.

test_that("the point summaries are the summaries of the draws", {
  fit <- bayes_fit()

  expect_equal(coef(fit), colMeans(fit$coefficient.draws))
  expect_equal(unname(fit$std.error),
               unname(apply(fit$coefficient.draws, 2, stats::sd)))
  expect_equal(vcov(fit), stats::var(fit$coefficient.draws))
  expect_equal(fit$sigma2, mean(fit$sigma2.draws))
  expect_equal(fitted(fit), as.vector(fit$X %*% coef(fit)))
  expect_equal(residuals(fit), fit$y - fitted(fit))

  expect_named(fit$acceptance, c("coefficients", "sigma2"))
  expect_true(all(fit$acceptance >= 0 & fit$acceptance <= 1))
  expect_true(all(fit$sigma2.draws > 0))
})

test_that("the estimated margins are proper distributions over the data", {
  dat <- make_sim_data()$data
  fit <- bayes_fit(data = dat)

  expect_named(fit$lambda.draws, c("z", "x"))
  for (v in c("z", "x")) {
    L <- fit$lambda.draws[[v]]
    expect_equal(dim(L), c(fit$ndraws, length(unique(dat[[v]]))))
    expect_true(all(L >= 0))
    expect_equal(rowSums(L), rep(1, fit$ndraws))
    expect_equal(fit$margins[[v]]$uv, sort(unique(dat[[v]])))
  }
})

test_that("confint() returns posterior quantiles of the draws", {
  fit <- bayes_fit()
  ci <- confint(fit)

  expect_equal(ci["z", ],
               stats::quantile(fit$draws[, "z"], c(0.025, 0.975)))
  expect_true(all(ci[, 1] <= ci[, 2]))

  # a narrower level gives a narrower interval
  ci90 <- confint(fit, level = 0.9)
  expect_equal(colnames(ci90), c("5%", "95%"))
  expect_gte(ci90["z", 1], ci["z", 1])
  expect_lte(ci90["z", 2], ci["z", 2])

  # parm selects, by name or by position
  expect_equal(rownames(confint(fit, parm = "z")), "z")
  expect_equal(rownames(confint(fit, parm = 1:2)), c("(Intercept)", "z"))
})

test_that("predict() handles factors in newdata through the stored levels", {
  dat <- make_sim_data()$data
  set.seed(8)
  dat$g <- factor(sample(c("a", "b", "c"), nrow(dat), replace = TRUE))
  fit <- bayes_fit(y ~ z | x + g, data = dat)

  # a subset holding only one level still gets the full contrast coding
  keep <- which(dat$g == "a")[1:5]
  expect_equal(predict(fit, newdata = dat[keep, ]), unname(fitted(fit)[keep]))
})

test_that("summary() tabulates the draws at the requested level", {
  fit <- bayes_fit()
  s <- summary(fit)

  expect_equal(colnames(s$coefficients),
               c("P. Mean", "P. Median", "Sd", "2.5%", "97.5%"))
  expect_equal(s$coefficients[, "P. Mean"], coef(fit))
  expect_equal(rownames(s$sigma2), "sigma2")
  expect_equal(unname(s$sigma2[, "P. Mean"]), fit$sigma2)

  # only the endogeneity correlations are tabulated; the correlations among
  # the regressors are counted, and the restricted ones are named
  expect_equal(rownames(s$rho), fit$rho.names)
  expect_equal(s$n.other, 1L)
  expect_equal(s$restricted, "rho(x*, xi*)")
  expect_equal(s$level, 0.95)
  expect_equal(s$n, fit$n)
  expect_equal(s$ndraws, fit$ndraws)

  s90 <- summary(fit, level = 0.9)
  expect_equal(colnames(s90$coefficients),
               c("P. Mean", "P. Median", "Sd", "5%", "95%"))
})

test_that("print() shows the call, the posterior means and the draw count", {
  dat <- make_sim_data(n = 60)$data
  set.seed(2)
  fit <- CopRegBAYES(y ~ z | x, data = dat, iterations = 200, burnin = 50,
                     thin = 2, verbose = FALSE)

  expect_output(print(fit), "Bayesian copula correction")
  expect_output(print(fit), "CopRegBAYES\\(formula = y ~ z \\| x")
  expect_output(print(fit), "Posterior means")
  expect_output(print(fit), "\\(Intercept\\)")
  expect_output(print(fit), "75 draws from 200 iterations")

  s <- summary(fit)
  expect_output(print(s), "Regression coefficients")
  expect_output(print(s), "Structural error variance")
  expect_output(print(s), "rho\\(z\\*, xi\\*\\)")
  expect_output(print(s), "95% credible")
  expect_output(print(s), "horseshoe")
})

test_that("the same seed reproduces the chain and a different one does not", {
  dat <- make_sim_data(n = 60)$data

  a <- bayes_fit(data = dat, iterations = 100, burnin = 20, seed = 7)
  b <- bayes_fit(data = dat, iterations = 100, burnin = 20, seed = 7)
  expect_identical(a$draws, b$draws)

  d <- bayes_fit(data = dat, iterations = 100, burnin = 20, seed = 8)
  expect_false(identical(a$draws, d$draws))
})

test_that("horseshoe = FALSE fits with the normal prior instead", {
  fit <- bayes_fit(horseshoe = FALSE)

  expect_false(fit$horseshoe)
  expect_true(all(is.finite(fit$draws)))
  expect_output(print(summary(fit)), "normal with inverse Gamma variance")
})

test_that("verbose = TRUE announces the run", {
  dat <- make_sim_data(n = 50)$data
  set.seed(1)
  expect_message(
    capture.output(
      CopRegBAYES(y ~ z | x, data = dat, iterations = 50, burnin = 10,
                  thin = 1, verbose = TRUE)),
    "Sampling 50 iterations"
  )
})
