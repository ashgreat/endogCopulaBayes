# The sampler calls several symbols unqualified (dCopula, rmvnorm,
# rinvwishart, ...); until the NAMESPACE gains import directives, the
# packages must be attached.
suppressMessages({
  library(copula)
  library(mvtnorm)
  library(LaplacesDemon)
  library(invgamma)
})

test_that("seeded smoke run produces a well-formed fit", {
  set.seed(321)
  n <- 80
  z <- rlnorm(n)
  x <- rnorm(n)
  e <- rnorm(n, sd = sqrt(2))
  dat <- data.frame(y = 2 - 4 * z + 6 * x + e, z = z, x = x)

  iterations <- 150
  burnin <- 50
  thin <- 2
  fit <- CopRegBayes(dat, iterations = iterations, burnin = burnin,
                     thin = thin, seed = 42)

  expect_s3_class(fit, "endog_copula_bayes")
  expect_equal(dim(fit$chain), c(iterations + 1, (7 + 2 * n) + 3))
  expect_equal(nrow(fit$posterior), length(seq(1, iterations + 1 - burnin, by = thin)))
  expect_equal(ncol(fit$posterior), (7 + 2 * n) + 3)
  expect_true(all(is.finite(fit$chain)))
  expect_true(all(is.finite(fit$posterior)))

  expect_equal(colnames(fit$chain)[1:7],
               c("beta_0", "beta_z", "beta_x", "sigma2",
                 "rho_zx", "rho_ze", "rho_xe"))
  expect_equal(colnames(fit$chain)[(2 * n + 8):(2 * n + 10)],
               c("hyper_a", "hyper_b1", "hyper_b2"))
  expect_equal(fit$n, n)
  expect_equal(fit$seed, 42)

  # supplying the same seed reproduces the chain
  fit2 <- CopRegBayes(dat, iterations = 25, burnin = 5, thin = 1, seed = 7)
  fit3 <- CopRegBayes(dat, iterations = 25, burnin = 5, thin = 1, seed = 7)
  expect_identical(fit2$chain, fit3$chain)

  # summary and print methods work
  s <- summary(fit)
  expect_s3_class(s, "summary.endog_copula_bayes")
  expect_equal(rownames(s$statistics), fit$parameters)
  expect_equal(colnames(s$statistics),
               c("mean", "sd", "median", "2.5%", "97.5%"))
  expect_true(all(is.finite(s$statistics)))
  expect_equal(s$statistics[, "mean"],
               colMeans(fit$posterior[, fit$parameters]))

  expect_output(print(fit), "Bayesian Gaussian copula regression")
  expect_output(print(fit), "beta_z")
  expect_output(print(s), "Posterior summary")
  expect_output(print(s), "97.5%")
})
