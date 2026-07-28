# Changelog

## endogCopulaBayes 0.2.0

### Breaking changes

- `CopRegBayes()` is replaced by
  [`CopRegBAYES()`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md),
  ported from the rewritten reference implementation. It takes the same
  two part formula used across the toolbox,
  `y ~ endogenous | exogenous`, and returns an object of class
  `"copregbayes"` with [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`validity()`](https://ashgreat.github.io/endogCopula/reference/validity.html)
  methods.
- `diagnostics()` is replaced by
  [`validity()`](https://ashgreat.github.io/endogCopula/reference/validity.html),
  which also reports chain convergence and can compute the Gelman-Rubin
  statistic on request.

### Changed results

- The 0.1.0 sampler had been ported from a superseded revision of the
  reference code, and it differed from the current reference
  implementation in all four sampler blocks. The coefficient step used a
  plain symmetric random walk instead of the Hastings corrected IWLS
  proposal. The variance step used an inverse gamma independence
  proposal instead of the Laplace proposal. Exogeneity was imposed by
  zeroing entries of an unrestricted inverse Wishart draw and rejecting,
  instead of the block recursive decomposition that now holds the zeros
  exactly. The marginal transform used the cumulative right endpoint
  instead of the cell midpoint. On the same data and seed the two
  samplers produce different posteriors, and the 0.2.0 sampler is the
  one that reproduces the published reference to 1e-12.

### Bug fixes

- A regressor with tied values used to raise a hard error, because the
  0.1.0 code added a vector of length N to a table of shorter length.
  The marginal probability masses now sit on the unique observed values,
  so ties are the normal case. A count, a dummy, a Likert item or a
  rounded price all fit without error.

### Dependencies

- The package now imports endogCopula for the two shared core helpers.
  The sampler itself is base R only, so copula, LaplacesDemon, invgamma
  and mvtnorm are dropped.
- The test suite is replaced, since the old tests targeted the removed
  API.
