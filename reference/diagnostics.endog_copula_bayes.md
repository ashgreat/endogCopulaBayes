# MCMC Convergence Diagnostics for a Bayesian Copula Regression Fit

Computes the effective sample size (ESS) and Geweke (1992) convergence
z-score for each of the named model parameters (the regression
coefficients, residual variance, copula correlations, and hyperprior
variances) using the thinned posterior sample stored in
`object$posterior`. The `2 * N` unnamed Dirichlet probability masses are
not included.

## Usage

``` r
# S3 method for class 'endog_copula_bayes'
diagnostics(object, ...)
```

## Arguments

- object:

  An object of class `endog_copula_bayes`.

- ...:

  Ignored.

## Value

A `data.frame` with one row per named parameter (row names give the
parameter names, respecting variable labels for a formula fit) and
columns:

- `ess`:

  Effective sample size.

- `geweke_z`:

  Geweke (1992) convergence z-score comparing the first 10% and last 50%
  of the thinned posterior sample.

- `n_draws`:

  Number of thinned posterior draws used.

## Details

Both statistics are implemented in base R using the same
spectral-density -at-zero estimator as
[`coda::spectrum0.ar()`](https://rdrr.io/pkg/coda/man/spectrum0.ar.html):
an AR(p) model is fit to each parameter's draws with the order selected
by AIC (via [`stats::ar()`](https://rdrr.io/r/stats/ar.html)), and the
resulting `var.pred / (1 - sum(ar))^2` is used as the long-run variance
in both the ESS ratio and the Geweke z-score denominator. A column whose
posterior draws are (numerically) constant has an undefined spectral
density; its ESS is reported as `0` and its Geweke z-score as `NA`,
rather than raising an error.

## See also

[`CopRegBayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md),
[`plot.endog_copula_bayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/plot.endog_copula_bayes.md)

## Examples

``` r
set.seed(1)
n <- 60
z <- rlnorm(n)
x <- rnorm(n)
dat <- data.frame(y = 2 - 4 * z + 6 * x + rnorm(n, sd = sqrt(2)),
                  z = z, x = x)
fit <- CopRegBayes(dat, iterations = 100, burnin = 20, thin = 2, seed = 42)
diagnostics(fit)
#>                ess   geweke_z n_draws
#> beta_0    1.505765 -2.4577231      41
#> beta_z    1.788386  2.9249946      41
#> beta_x    4.991725 -2.1334245      41
#> sigma2   41.000000 -0.6542289      41
#> rho_zx   63.540290  0.1953849      41
#> rho_ze    2.698485 -3.8641608      41
#> rho_xe    0.000000         NA      41
#> hyper_a  41.000000 -0.9103161      41
#> hyper_b1 41.000000 -0.9979605      41
#> hyper_b2 41.000000  3.4375755      41
```
