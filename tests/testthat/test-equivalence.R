# Exact-chain equivalence against the reference implementation in
# Copula-based-endogeneity-corrections-main/CopRegBAYES.R. The reference
# sampler reads the sample size N from its enclosing (global) scope, so it is
# assigned into the environment holding the parsed reference functions.
suppressMessages({
  library(copula)
  library(mvtnorm)
  library(LaplacesDemon)
  library(invgamma)
})

test_that("metropolis_Gibbs_MCMC1 reproduces the reference chain exactly", {
  skip_if_no_reference()
  ref <- load_reference_functions("CopRegBAYES.R")

  # simulate a small data set following the reference DGP
  set.seed(99)
  n <- 60
  sig2 <- 5
  sigma_mat1 <- matrix(c(1, .7, 0,
                         .7, 1, .3,
                         0, .3, 1), ncol = 3, nrow = 3, byrow = TRUE)
  eps <- mvtnorm::rmvnorm(n = n, mean = c(0, 0, 0), sigma = sigma_mat1,
                          method = "chol")
  dat <- data.frame(const = rep(1, n))
  dat$z <- (qlnorm(p = pnorm(eps[, 2])) - exp(.5)) / sqrt((exp(1) - 1) * exp(1))
  dat$x <- qnorm(pnorm(eps[, 3]))
  dat$e <- qnorm(pnorm(eps[, 1]), sd = sqrt(sig2))
  dat$y <- as.vector(cbind(dat$const, dat$z, dat$x) %*% c(2, -4, 6)) + dat$e

  # identical starting values for both samplers
  mod1 <- lm(y ~ z + x, dat)
  startvalue <- c(mod1$coefficients, var(mod1$residuals), rep(0, 3),
                  c(LaplacesDemon::rdirichlet(n = 1, rep(1, n))),
                  c(LaplacesDemon::rdirichlet(n = 1, rep(1, n))),
                  rep(1000, 3))

  iterations <- 100

  # reference sampler uses the global N; provide it in its environment
  assign("N", n, envir = ref)
  set.seed(42)
  ref_chain <- ref$metropolis_Gibbs_MCMC1(startvalue = startvalue,
                                          iterations = iterations,
                                          data = dat)

  set.seed(42)
  pkg_chain <- metropolis_Gibbs_MCMC1(startvalue = startvalue,
                                      iterations = iterations,
                                      data = dat)

  expect_equal(dim(pkg_chain), dim(ref_chain))
  expect_equal(pkg_chain, ref_chain, tolerance = 1e-10)
})
