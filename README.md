# endogCopulaBayes

<!-- badges: start -->
[![R-CMD-check](https://github.com/ashgreat/endogCopulaBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ashgreat/endogCopulaBayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`endogCopulaBayes` implements the Bayesian Gaussian copula endogeneity
correction of Haschka (2022b; published 2026 in the *Oxford Bulletin of
Economics and Statistics* as "Bayesian Inference for Joint Estimation Models
Using Copulas to Handle Endogenous Regressors"). The estimator handles an
endogenous regressor without instruments by modelling the joint distribution
of the regressors and the regression error with a Gaussian copula and
sampling all unknowns — regression coefficients, residual variance, copula
correlations, and the (Dirichlet) probability masses describing the marginal
distributions of the regressors — with a Metropolis-within-Gibbs MCMC
sampler.

The code is a faithful port of the published replication scripts
(`CopRegBAYES.R`) from the copula-based endogeneity corrections research
code; the packaged sampler reproduces the reference chain draw for draw.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ashgreat/endogCopulaBayes")
```

## Data format

`CopRegBayes()` expects a data frame with exactly these columns (extra
columns are ignored):

| Column | Role |
|--------|------|
| `y`    | Continuous dependent variable |
| `z`    | Continuous endogenous regressor (must be non-normally distributed for identification) |
| `x`    | Continuous exogenous regressor |

Notes:

- The intercept column (`const`) is added internally — do **not** include
  one in your data.
- Rows with missing values in `y`, `z`, or `x` are dropped before sampling.
- The model estimated is `y = beta_0 + beta_z * z + beta_x * x + e`, with a
  Gaussian copula linking `z`, `x`, and `e`; the correlation between `x` and
  `e` is restricted to zero (exogeneity of `x`).

## Worked example

```r
library(endogCopulaBayes)

# Simulate data with an endogenous, log-normally distributed regressor z
set.seed(1)
n <- 200
sigma <- matrix(c(1, .7, 0,
                  .7, 1, .3,
                  0, .3, 1), nrow = 3, byrow = TRUE)
eps <- mvtnorm::rmvnorm(n, sigma = sigma)
z <- qlnorm(pnorm(eps[, 2]))          # endogenous (correlated with the error)
x <- qnorm(pnorm(eps[, 3]))           # exogenous
e <- qnorm(pnorm(eps[, 1]), sd = sqrt(2))
dat <- data.frame(y = 2 - 4 * z + 6 * x + e, z = z, x = x)

# Fit (use many more iterations in real applications)
fit <- CopRegBayes(dat, iterations = 10000, burnin = 2000, thin = 10,
                   seed = 42)

fit             # posterior means of the named parameters
summary(fit)    # mean, sd, median, and 95% credible intervals
```

The returned object (class `endog_copula_bayes`) carries the full MCMC
`chain` and the burned-in, thinned `posterior` matrix. Columns 1–7 are the
named model parameters (`beta_0`, `beta_z`, `beta_x`, `sigma2`, `rho_zx`,
`rho_ze`, `rho_xe`), followed by `2 * N` unnamed Dirichlet probability
masses for the marginals of `z` and `x`, and three hyperprior variances
(`hyper_a`, `hyper_b1`, `hyper_b2`). See `?CopRegBayes` for details.

## Related packages

- [endogCopula](https://github.com/ashgreat/endogCopula) — cross-sectional
  copula endogeneity corrections (Park & Gupta 2012; Yang et al. 2025; Hu et
  al. 2025; Breitung et al. 2024; Haschka 2024; Liengaard et al. 2025).
- [endogCopulaPanel](https://github.com/ashgreat/endogCopulaPanel) — the
  fixed-effects panel maximum-likelihood estimator of Haschka (2022).

## References

- Haschka, R. E. (2026). Bayesian Inference for Joint Estimation Models
  Using Copulas to Handle Endogenous Regressors. *Oxford Bulletin of
  Economics and Statistics*. (Working paper version: Haschka 2022b,
  <https://ssrn.com/abstract=4235194>.)
- Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
  estimation using copulas. *Marketing Science* 31(4), 567–586.

## License

MIT + file LICENSE. See [LICENSE.md](LICENSE.md) for the full text.
