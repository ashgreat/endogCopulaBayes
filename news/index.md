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

### Differences from the reference implementation

- The package deliberately differs from the reference implementation in
  three places. All three sit inside
  [`validity()`](https://ashgreat.github.io/endogCopula/reference/validity.html).
  None of them changes any estimate: the sampler itself, and so every
  coefficient, draw and posterior summary, is unchanged.
- The effective sample size function `.bayes_ess()` truncates its
  autocorrelation sum at the first negative lag unconditionally. The
  reference implementation only truncates when that lag is not the first
  one, which let it report an effective sample size larger than the
  number of retained draws, up to 390 times the number of draws, on 45%
  of chains in a simulation check. The package instead always truncates
  at the first negative lag, so the reported value never exceeds the
  number of retained draws.
- [`validity.copregbayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/validity.copregbayes.md)
  captures the calling frame with `where <- parent.frame()` in its own
  body, before the [`lapply()`](https://rdrr.io/r/base/lapply.html) that
  runs the extra chains for the Gelman-Rubin statistic, and evaluates
  the original call there. The reference implementation calls
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html) inside the
  [`lapply()`](https://rdrr.io/r/base/lapply.html) callback instead,
  which resolves to the [`lapply()`](https://rdrr.io/r/base/lapply.html)
  call’s own frame rather than the frame
  [`validity()`](https://ashgreat.github.io/endogCopula/reference/validity.html)
  was called from. That made `validity(fit, chains = 2)` fail with
  “object ‘d’ not found” whenever it was called from inside a function
  whose data was not in the global environment. The package instead
  finds the data reliably regardless of where
  [`validity()`](https://ashgreat.github.io/endogCopula/reference/validity.html)
  is called from.
- `.bayes_ess()` returns `NA_real_` for a chain that never moves, since
  [`stats::acf()`](https://rdrr.io/r/stats/acf.html) of a constant
  series is `NaN` and the truncation above never fires for it. The
  reference implementation returns `NaN` in that case. The package
  instead reports `NA`, matching what `.bayes_geweke()` already returns
  for the same input, so the convergence table reads as “not available”
  rather than showing `NaN`.

### Dependencies

- The package now imports endogCopula for the two shared core helpers.
  The sampler itself is base R only, so copula, LaplacesDemon, invgamma
  and mvtnorm are dropped.
- The test suite is replaced, since the old tests targeted the removed
  API.
