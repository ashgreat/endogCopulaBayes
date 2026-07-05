# The sampler calls several symbols unqualified (dCopula, rmvnorm,
# rinvwishart, ...); until the NAMESPACE gains import directives, the
# packages must be attached.
suppressMessages({
  library(copula)
  library(mvtnorm)
  library(LaplacesDemon)
  library(invgamma)
})

test_that("formula interface matches the data-frame interface exactly", {
  set.seed(321)
  n <- 40
  price <- rlnorm(n)
  income <- rnorm(n)
  e <- rnorm(n, sd = sqrt(2))
  sales <- 2 - 4 * price + 6 * income + e

  dat <- data.frame(sales = sales, price = price, income = income)

  fit_formula <- CopRegBayes(sales ~ price | income, data = dat,
                             iterations = 25, burnin = 5, thin = 1, seed = 7)

  dat_renamed <- data.frame(y = sales, z = price, x = income)
  fit_df <- CopRegBayes(dat_renamed, iterations = 25, burnin = 5, thin = 1,
                        seed = 7)

  expect_identical(fit_formula$chain, fit_df$chain)
  expect_identical(fit_formula$posterior, fit_df$posterior)
})

test_that("formula method records original variable names", {
  set.seed(4)
  n <- 30
  dat <- data.frame(sales = rnorm(n), price = rlnorm(n), income = rnorm(n))

  fit <- CopRegBayes(sales ~ price | income, data = dat,
                     iterations = 10, burnin = 2, thin = 1, seed = 1)

  expect_type(fit$variables, "list")
  expect_equal(fit$variables$response, "sales")
  expect_equal(fit$variables$endogenous, "price")
  expect_equal(fit$variables$exogenous, "income")
})

test_that("legacy data-frame interface leaves variables NULL and labels unchanged", {
  set.seed(5)
  n <- 25
  dat <- data.frame(y = rnorm(n), z = rlnorm(n), x = rnorm(n))

  fit <- CopRegBayes(dat, iterations = 10, burnin = 2, thin = 1, seed = 1)

  expect_null(fit$variables)
  s <- summary(fit)
  expect_true(all(c("beta_z", "beta_x") %in% rownames(s$statistics)))
})

test_that("summary/print label beta_z / beta_x with the real variable names", {
  set.seed(6)
  n <- 30
  dat <- data.frame(sales = rnorm(n), price = rlnorm(n), income = rnorm(n))

  fit <- CopRegBayes(sales ~ price | income, data = dat,
                     iterations = 10, burnin = 2, thin = 1, seed = 1)

  s <- summary(fit)
  expect_true("price" %in% rownames(s$statistics))
  expect_true("income" %in% rownames(s$statistics))
  expect_false("beta_z" %in% rownames(s$statistics))
  expect_false("beta_x" %in% rownames(s$statistics))

  expect_output(print(fit), "price")
  expect_output(print(fit), "income")
  expect_output(print(s), "price")
  expect_output(print(s), "income")
})

test_that("formula interface errors when there is no '|' separator", {
  dat <- data.frame(y = rnorm(10), a = rnorm(10), b = rnorm(10))
  expect_error(
    CopRegBayes(y ~ a + b, data = dat),
    "must separate endogenous and exogenous regressors"
  )
})

test_that("formula interface errors with two endogenous regressors", {
  dat <- data.frame(y = rnorm(10), a = rnorm(10), b = rnorm(10), c = rnorm(10))
  expect_error(
    CopRegBayes(y ~ a + b | c, data = dat),
    "exactly one endogenous and exactly one exogenous"
  )
})

test_that("formula interface errors with zero exogenous regressors", {
  dat <- data.frame(y = rnorm(10), a = rnorm(10))
  expect_error(
    CopRegBayes(y ~ a | 1, data = dat),
    "exactly one endogenous and exactly one exogenous"
  )
})

test_that("formula interface errors on transformed regressors instead of silently dropping them", {
  dat <- data.frame(y = rnorm(10), a = rlnorm(10), b = rnorm(10))
  expect_error(
    CopRegBayes(y ~ log(a) | b, data = dat),
    "single.*untransformed column name"
  )
  expect_error(
    CopRegBayes(y ~ I(a^2) | b, data = dat),
    "single.*untransformed column name"
  )
  expect_error(
    CopRegBayes(y ~ a | log(b), data = dat),
    "single.*untransformed column name"
  )
})

test_that("formula interface errors on intercept-removal syntax", {
  dat <- data.frame(y = rnorm(10), a = rnorm(10), b = rnorm(10))
  expect_error(
    CopRegBayes(y ~ a - 1 | b, data = dat),
    "intercept removal"
  )
  expect_error(
    CopRegBayes(y ~ 0 + a | b, data = dat),
    "intercept removal"
  )
  expect_error(
    CopRegBayes(y ~ a | b - 1, data = dat),
    "intercept removal"
  )
  expect_error(
    CopRegBayes(y ~ a | 0 + b, data = dat),
    "intercept removal"
  )
})

test_that("formula interface errors on a transformed response", {
  dat <- data.frame(y = rlnorm(10), a = rnorm(10), b = rnorm(10))
  expect_error(
    CopRegBayes(log(y) ~ a | b, data = dat),
    "single untransformed column name"
  )
  expect_error(
    CopRegBayes(cbind(y, a) ~ a | b, data = dat),
    "single untransformed column name"
  )
})

test_that("formula interface errors when a variable is missing from data", {
  dat <- data.frame(y = rnorm(10), a = rnorm(10))
  expect_error(
    CopRegBayes(y ~ a | missing_var, data = dat),
    "Missing variables in data: missing_var"
  )
})

test_that("formula interface errors on factor variables", {
  n <- 20
  dat <- data.frame(y = rnorm(n), a = rnorm(n),
                    b = factor(sample(c("lo", "hi"), n, replace = TRUE)))
  expect_error(
    CopRegBayes(y ~ a | b, data = dat),
    "Factor variables are not supported: b"
  )
})
