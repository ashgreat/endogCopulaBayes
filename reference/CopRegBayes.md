# Bayesian Gaussian Copula Sampler

Metropolis-within-Gibbs sampler for the copula-based endogeneity
correction of Haschka. The chain columns are laid out as in the
reference implementation: columns 1–3 hold the regression coefficients
(`beta_0`, `beta_z`, `beta_x`), column 4 the residual variance
(`sigma2`), columns 5–7 the copula correlations (`rho_zx` between
endogenous and exogenous regressor, `rho_ze` between endogenous
regressor and error, `rho_xe` between exogenous regressor and error, the
latter forced to zero), columns `8:(N + 7)` the Dirichlet probability
masses for the distribution of `z`, columns `(N + 8):(2 * N + 7)` the
Dirichlet probability masses for the distribution of `x`, and the final
three columns the hyperprior variances of the normal priors on the
regression coefficients (`hyper_a`, `hyper_b1`, `hyper_b2`). The
Dirichlet mass columns are left unnamed.

## Usage

``` r
CopRegBayes(x, ...)

# S3 method for class 'formula'
CopRegBayes(
  formula,
  data,
  iterations = 10000,
  burnin = 2000,
  thin = 10,
  startvalue = NULL,
  seed = NULL,
  ...
)

# S3 method for class 'data.frame'
CopRegBayes(
  x,
  iterations = 10000,
  burnin = 2000,
  thin = 10,
  startvalue = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  Either a `formula` of the form `y ~ endog | exog`, or a `data.frame`
  with columns `y`, `z`, and `x` (legacy interface).

- ...:

  Passed on to methods (currently unused).

- formula:

  A two-part formula of the form `y ~ endog | exog`, with exactly one
  endogenous regressor before the bar and exactly one exogenous
  regressor after it (the Haschka (2026) sampler supports exactly one of
  each).

- data:

  A `data.frame` containing the variables referenced in `formula`
  (formula method) or the columns `y`, `z`, and `x` (data-frame method).

- iterations:

  Total number of MCMC iterations.

- burnin:

  Number of initial iterations to discard.

- thin:

  Thinning interval applied after burn-in.

- startvalue:

  Optional numeric vector of starting values; if `NULL`, a default based
  on OLS estimates and Dirichlet draws is used.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before any random
  number is drawn; if `NULL` (default) the current RNG state is used.

## Value

An object of class `endog_copula_bayes`: a list with components

- `chain`:

  Numeric matrix of dimension `(iterations + 1) x (2 * N + 10)` holding
  the full MCMC chain, one draw per row. The columns follow the layout
  described above: the seven named model parameters, then the `2 * N`
  unnamed Dirichlet probability masses (`N` for `z`, `N` for `x`), then
  the three hyperprior variances.

- `posterior`:

  Numeric matrix with the same columns as `chain`, after discarding the
  first `burnin` rows and keeping every `thin`-th remaining row.

- `parameters`:

  Character vector with the names of the ten tracked (named) parameters
  summarised by
  [`summary.endog_copula_bayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/summary.endog_copula_bayes.md).

- `n`:

  Number of complete observations used.

- `iterations`, `burnin`, `thin`, `seed`:

  The sampler settings as supplied.

- `variables`:

  `NULL` for the data-frame interface; for the formula interface, a list
  with the original `response`, `endogenous`, and `exogenous` variable
  names.

## Details

`CopRegBayes()` is an S3 generic. Use the formula interface
(`y ~ endog | exog`) or the legacy data-frame interface (a data frame
with columns `y`, `z`, and `x`); both dispatch to the same underlying
sampler so a seeded fit is identical either way. The formula interface
additionally records the original variable names (`$variables`) and uses
them to label the `beta_z` / `beta_x` rows in
[`print.endog_copula_bayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/print.endog_copula_bayes.md)
and
[`summary.endog_copula_bayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/summary.endog_copula_bayes.md).

## See also

[`summary.endog_copula_bayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/summary.endog_copula_bayes.md)
and
[`print.endog_copula_bayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/print.endog_copula_bayes.md)
for posterior summaries.

## Examples

``` r
set.seed(1)
n <- 60
z <- rlnorm(n)
x <- rnorm(n)
dat <- data.frame(y = 2 - 4 * z + 6 * x + rnorm(n, sd = sqrt(2)),
                  z = z, x = x)
fit <- CopRegBayes(dat, iterations = 100, burnin = 20, thin = 2, seed = 42)
fit
#> Bayesian Gaussian copula regression (Haschka)
#> Observations: 60
#> Iterations: 100 (burn-in 20, thinning 2, 41 posterior draws)
#> 
#> Posterior means:
#>    beta_0    beta_z    beta_x    sigma2    rho_zx    rho_ze    rho_xe   hyper_a 
#>    2.1708   -4.1770    5.6646    2.2923   -0.0621    0.1038    0.0000  230.4490 
#>  hyper_b1  hyper_b2 
#> 1744.2163  547.1008 
summary(fit)
#> Bayesian Gaussian copula regression (Haschka)
#> Observations: 60
#> Iterations: 100 (burn-in 20, thinning 2, 41 posterior draws)
#> 
#> Posterior summary:
#>               mean        sd  median    2.5%     97.5%
#> beta_0      2.1708    0.4746  1.9885  1.7091    3.0271
#> beta_z     -4.1770    0.2729 -4.0512 -4.6916   -3.8495
#> beta_x      5.6646    0.1503  5.6206  5.4646    5.9248
#> sigma2      2.2923    0.3857  2.1896  1.5984    3.0272
#> rho_zx     -0.0621    0.1324 -0.0587 -0.2963    0.1994
#> rho_ze      0.1038    0.2357  0.0090 -0.2172    0.5247
#> rho_xe      0.0000    0.0000  0.0000  0.0000    0.0000
#> hyper_a   230.4490  586.0057  8.4011  0.8192 2364.4946
#> hyper_b1 1744.2163 9555.6623 56.5276  4.2783 4077.6518
#> hyper_b2  547.1008 1882.5218 79.5511 11.3015 1344.0010

# \donttest{
# Equivalent formula interface (variable names are used as labels)
fit_f <- CopRegBayes(y ~ z | x, data = dat, iterations = 100, burnin = 20,
                     thin = 2, seed = 42)
summary(fit_f)
#> Bayesian Gaussian copula regression (Haschka)
#> Observations: 60
#> Iterations: 100 (burn-in 20, thinning 2, 41 posterior draws)
#> 
#> Posterior summary:
#>               mean        sd  median    2.5%     97.5%
#> beta_0      2.1708    0.4746  1.9885  1.7091    3.0271
#> z          -4.1770    0.2729 -4.0512 -4.6916   -3.8495
#> x           5.6646    0.1503  5.6206  5.4646    5.9248
#> sigma2      2.2923    0.3857  2.1896  1.5984    3.0272
#> rho_zx     -0.0621    0.1324 -0.0587 -0.2963    0.1994
#> rho_ze      0.1038    0.2357  0.0090 -0.2172    0.5247
#> rho_xe      0.0000    0.0000  0.0000  0.0000    0.0000
#> hyper_a   230.4490  586.0057  8.4011  0.8192 2364.4946
#> hyper_b1 1744.2163 9555.6623 56.5276  4.2783 4077.6518
#> hyper_b2  547.1008 1882.5218 79.5511 11.3015 1344.0010
# }
```
