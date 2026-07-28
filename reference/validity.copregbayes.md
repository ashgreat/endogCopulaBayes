# Identification and convergence checks for a copregbayes fit

[`validity()`](https://ashgreat.github.io/endogCopula/reference/validity.html)
for a
[`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md)
fit reports two things at once. First, the identification requirement of
the model: the copula correction is identified off the nonnormality of
the endogenous regressors, since under normality \\E(e\|z)\\ is linear
and the copula cannot separate regressor variation from error variation;
the posterior of the endogeneity correlations (the correlation between
an endogenous regressor's normal score and that of the structural error)
is reported alongside, since the Bayesian version is expected to degrade
gracefully – concentrating near zero rather than producing confident
nonsense – when identification is weak. Second, convergence of the MCMC
chain itself: Geweke's statistic, the effective sample size and the
lag-one autocorrelation, all computed from the chain that is already
there, plus – on request, since it needs several chains from dispersed
starting values – the Gelman-Rubin statistic, which costs one further
run of the sampler per additional chain.

[`print()`](https://rdrr.io/r/base/print.html) on the result formats all
of that.

## Usage

``` r
# S3 method for class 'copregbayes'
validity(object, chains = FALSE, power = 0.8, verbose = interactive(), ...)

# S3 method for class 'copregbayes.validity'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An object of class `"copregbayes"`, as returned by
  [`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md).

- chains:

  Whether to also compute the Gelman-Rubin statistic from dispersed
  starting values. `FALSE` (the default) skips it. `TRUE` runs 4
  additional chains; a number runs that many additional chains instead
  (at least 2). Each additional chain re-runs the sampler with the same
  call but a fresh, dispersed starting point drawn from
  \\N\\(coefficients, posterior sd) for the regression coefficients, an
  LKJ(1) draw for the copula correlation matrix, and a fresh Dir(1) draw
  per regressor for the probability masses.

- power:

  Statistical power used for the sample-size-dependent nonnormality
  thresholds (Becker, Proksch & Ringle 2022). Defaults to `0.8`.

- verbose:

  If `TRUE`, print a message before each additional chain required by
  `chains`. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- ...:

  Not used; present for S3 method consistency.

- x:

  An object of class `"copregbayes.validity"`, as returned by
  `validity.copregbayes`.

- digits:

  Number of significant digits to print.

## Value

For `validity.copregbayes`, an object of class `"copregbayes.validity"`,
a list with `nonnormality` (a table of skewness, excess kurtosis,
Anderson-Darling and Cramer-von Mises statistics and a
Kolmogorov-Smirnov p-value for each endogenous regressor, with pass/fail
flags against two thresholds), `thresholds` (the
sample-size-and-power-dependent thresholds used), `endogeneity`
(posterior mean, 95% interval and \\P(\rho \> 0)\\ for each endogeneity
correlation), `convergence` (a data frame of Geweke, effective sample
size and lag-one autocorrelation for every retained parameter),
`gelman.rubin` (a named vector of Gelman-Rubin statistics, or `NULL` if
`chains = FALSE`), `acceptance`, `ndraws` and `endo.names`.

For `print.copregbayes.validity`, `x` is returned invisibly; called for
its side effect of printing.

## References

Haschka, R. E. (2025). Bayesian inference for joint estimation models
using copulas to handle endogenous regressors. *Oxford Bulletin of
Economics and Statistics*.
[doi:10.1111/obes.70023](https://doi.org/10.1111/obes.70023)

## Examples

``` r
set.seed(1)
n <- 60
x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
                    iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
validity(fit)
#> 
#> Checks for the Bayesian copula correction
#> 
#> [1] Nonnormality of the endogenous regressors
#>     Identification comes from the nonlinearity of E(e|z), which is
#>     present only when z is nonnormal; under normality the copula cannot
#>     tell regressor variation from error variation.
#>   skewness ex.kurtosis    AD  CvM  KS p Yang ok Becker ok
#> z   -0.691       0.741 0.451 0.06 0.812   FALSE     FALSE
#> 
#> [2] Posterior of the endogeneity correlations
#>              P. Mean 2.5%    97.5%   P(rho > 0)
#> rho(z*, xi*) -0.4564 -0.7113  0.1041  0.1000   
#>     An interval covering zero says the data carry no evidence of
#>     endogeneity, which is the honest reading; it is not a test decision.
#> 
#> [3] Convergence of the chain
#>     Geweke compares the first tenth with the last half of the draws and
#>     is a z statistic, so |Geweke| > 2 is a warning sign. ESS is the
#>     effective number of independent draws behind 30 retained ones.
#>                Geweke     ESS    AC(1)
#> (Intercept)    4.5229  18.240  0.29615
#> z            -11.1126   5.897  0.76881
#> x             13.0992   5.391  0.66270
#> sigma2        -7.5163   5.991  0.58845
#> rho(z*, x*)   -0.9601 798.730 -0.06026
#> rho(z*, xi*)   7.1434   6.304  0.76326
#>     Acceptance rates: coefficients 0.835, sigma2 0.865
#> 
#> [4] Gelman-Rubin was not computed
#>     It needs several chains from dispersed starts, so it costs one
#>     further run of the sampler per chain. Pass chains = TRUE, or a
#>     number, to compute it.
#> 
#> Source: Haschka (2025), Oxford Bulletin of Economics and Statistics
```
