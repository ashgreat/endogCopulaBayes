# Predict from a copregbayes fit

Predictions from the structural regression only: the copula terms never
enter [`predict()`](https://rdrr.io/r/stats/predict.html) or
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) because they
are endogeneity controls, not part of the causal model being predicted.

## Usage

``` r
# S3 method for class 'copregbayes'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  An object of class `"copregbayes"`, as returned by
  [`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md).

- newdata:

  An optional `data.frame` of new predictor values. If `NULL` (the
  default), the fitted values of the original data are returned.

- ...:

  Not used; present for S3 method consistency.

## Value

A numeric vector of predictions, one per row of `newdata` (or of the
original data if `newdata` is `NULL`), computed at the posterior mean
coefficients.

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
dat <- data.frame(y = y, z = z, x = x)
fit <- CopRegBAYES(y ~ z | x, data = dat,
                    iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
predict(fit, newdata = dat[1:5, ])
#> [1] 3.7202734 1.1811879 0.2712549 4.3285814 0.2603980
```
