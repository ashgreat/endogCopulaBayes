# Small simulated data set with a known true coefficient vector. The
# endogenous regressor z is non-normal (chi-squared) and correlated with the
# structural error: z is built by rank-matching a draw from rchisq() to a
# standard normal eps[, 1] that is itself correlated with the structural
# error eps[, 2], so z keeps an exact chi-squared marginal while carrying the
# same rank correlation with the error that the Gaussian copula model is
# designed to pick up.
make_sim_data <- function(n = 180, seed = 1, rho = 0.6,
                           beta0 = 1, beta_z = 2, beta_x = 0.5) {
  set.seed(seed)
  x   <- stats::rnorm(n)
  eps <- matrix(stats::rnorm(2 * n), n, 2)
  eps[, 2] <- rho * eps[, 1] + sqrt(1 - rho^2) * eps[, 2]
  z <- sort(stats::rchisq(n, df = 3))[rank(eps[, 1])] + 0.3 * x
  e <- eps[, 2]
  y <- beta0 + beta_z * z + beta_x * x + e
  list(
    data = data.frame(y = y, z = z, x = x),
    true = c(`(Intercept)` = beta0, z = beta_z, x = beta_x)
  )
}
