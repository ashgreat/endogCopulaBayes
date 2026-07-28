# Extract components of a copregbayes fit

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the posterior mean
of the regression coefficients,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) their posterior covariance
matrix, [`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
observations used, [`formula()`](https://rdrr.io/r/stats/formula.html)
the model formula,
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) the fitted
values at the posterior mean coefficients, and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) the
corresponding residuals.

## Usage

``` r
# S3 method for class 'copregbayes'
coef(object, ...)

# S3 method for class 'copregbayes'
vcov(object, ...)

# S3 method for class 'copregbayes'
nobs(object, ...)

# S3 method for class 'copregbayes'
formula(x, ...)

# S3 method for class 'copregbayes'
fitted(object, ...)

# S3 method for class 'copregbayes'
residuals(object, ...)
```

## Arguments

- object:

  An object of class `"copregbayes"`, as returned by
  [`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md).

- ...:

  Not used; present for S3 method consistency.

- x:

  An object of class `"copregbayes"` (the argument name used by the
  [`formula()`](https://rdrr.io/r/stats/formula.html) method).

## Value

For `coef.copregbayes`, a named numeric vector: the posterior mean of
each regression coefficient.

For `vcov.copregbayes`, a numeric matrix: the posterior covariance of
the regression coefficients.

For `nobs.copregbayes`, an integer: the number of observations.

For `formula.copregbayes`, the two-part model formula.

For `fitted.copregbayes`, a numeric vector of fitted values.

For `residuals.copregbayes`, a numeric vector of residuals (response
minus fitted values, at the posterior mean coefficients).

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
coef(fit); vcov(fit); nobs(fit); formula(fit)
#> (Intercept)           z           x 
#>   0.8558353   1.7525212   0.3936036 
#>             (Intercept)           z           x
#> (Intercept)  0.03493359 -0.03511199  0.01822975
#> z           -0.03511199  0.25725760 -0.17790201
#> x            0.01822975 -0.17790201  0.15221237
#> [1] 60
#> y ~ z | x
#> <environment: 0x55cbec1b4538>
head(fitted(fit)); head(residuals(fit))
#> [1]  3.7202734  1.1811879  0.2712549  4.3285814  0.2603980 -0.5741297
#>          1          2          3          4          5          6 
#> -2.0775207  1.4898976 -0.4673522 -0.2895742  0.5551536  0.8346516 
```
