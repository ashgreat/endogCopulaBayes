test_that("CopRegBayes errors when required columns are missing", {
  dat <- data.frame(y = rnorm(10), z = rnorm(10))
  expect_error(CopRegBayes(dat), "Missing variables in data: x")

  dat2 <- data.frame(a = rnorm(10))
  expect_error(CopRegBayes(dat2), "Missing variables in data: y, z, x")
})

test_that("CopRegBayes errors when no complete observations remain", {
  dat <- data.frame(y = c(NA_real_, NA_real_),
                    z = c(1, NA_real_),
                    x = c(NA_real_, 2))
  expect_error(CopRegBayes(dat), "No complete observations")
})

test_that("CopRegBayes errors on a startvalue of the wrong length", {
  set.seed(1)
  n <- 20
  dat <- data.frame(z = rlnorm(n), x = rnorm(n))
  dat$y <- 1 + 2 * dat$z - dat$x + rnorm(n)
  expect_error(
    CopRegBayes(dat, iterations = 5, startvalue = rep(0.1, 5)),
    "Length of 'startvalue' does not match expected dimension"
  )
})
