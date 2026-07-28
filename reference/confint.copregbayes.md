# Credible intervals for a copregbayes fit

Posterior credible intervals, taken directly as quantiles of the
retained draws in `object$draws` – not a normal approximation to them,
and no asymptotic argument is involved anywhere.

## Usage

``` r
# S3 method for class 'copregbayes'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `"copregbayes"`, as returned by
  [`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md).

- parm:

  Optional character or integer vector selecting which parameters
  (columns of `object$draws`) to return intervals for. If missing, all
  parameters are returned.

- level:

  Credible level, e.g. `0.95` for a 95% interval. Defaults to `0.95`.

- ...:

  Not used; present for S3 method consistency.

## Value

A numeric matrix with one row per parameter and two columns, the lower
and upper quantiles of the posterior draws at the requested credible
level.

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
confint(fit, level = 0.9)
#>                       5%        95%
#> (Intercept)   0.58081263 1.08908201
#> z             0.80466858 2.35495730
#> x            -0.03377693 1.14262702
#> sigma2        1.00784023 3.39397825
#> rho(z*, x*)   0.39537919 0.62547322
#> rho(z*, xi*) -0.68733513 0.05093047
```
