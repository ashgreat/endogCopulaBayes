# Plot a Bayesian Copula Regression Fit

Base-graphics traceplots and posterior density plots for the named
parameters of a fitted `endog_copula_bayes` object (the `2 * N` unnamed
Dirichlet probability masses are excluded).

## Usage

``` r
# S3 method for class 'endog_copula_bayes'
plot(x, params = NULL, which = c("trace", "density"), ...)
```

## Arguments

- x:

  An object of class `endog_copula_bayes`.

- params:

  Character vector of parameter names to plot (after variable-label
  substitution for a formula fit); `NULL` (default) plots all named
  parameters in `x$parameters`.

- which:

  Character vector selecting which diagnostic plots to draw: `"trace"`
  for traceplots, `"density"` for posterior density plots, or both (the
  default).

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`CopRegBayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md),
[`diagnostics()`](https://ashgreat.github.io/endogCopulaBayes/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
n <- 60
z <- rlnorm(n)
x <- rnorm(n)
dat <- data.frame(y = 2 - 4 * z + 6 * x + rnorm(n, sd = sqrt(2)),
                  z = z, x = x)
fit <- CopRegBayes(dat, iterations = 100, burnin = 20, thin = 2, seed = 42)
plot(fit, params = c("beta_z", "beta_x"))
```
