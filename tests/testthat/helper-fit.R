# Short seeded run, for the many tests that only need a well-formed fit.
# CopRegBAYES() has no 'seed' argument -- it draws from the global stream --
# so the seed is set here rather than passed through. Tests that depend on the
# recorded call (print, and the extra chains validity() runs for Gelman-Rubin)
# call CopRegBAYES() directly instead.
bayes_fit <- function(formula = y ~ z | x, data = make_sim_data()$data,
                      iterations = 300, burnin = 50, thin = 2, seed = 42,
                      ...) {
  force(data)
  set.seed(seed)
  CopRegBAYES(formula, data = data, iterations = iterations, burnin = burnin,
              thin = thin, verbose = FALSE, ...)
}
