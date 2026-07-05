# Summarise a Bayesian Copula Regression Fit

Posterior summaries (mean, standard deviation, median, and 95% credible
interval) of the named model parameters, computed from the thinned
posterior sample.

## Usage

``` r
# S3 method for class 'endog_copula_bayes'
summary(object, ...)
```

## Arguments

- object:

  An object of class `endog_copula_bayes`.

- ...:

  Ignored.

## Value

An object of class `summary.endog_copula_bayes`: a list with the matrix
of summary `statistics` (one row per named parameter, columns `mean`,
`sd`, `median`, `2.5%`, and `97.5%`), the number of observations `n`,
the number of retained posterior draws `n_draws`, and the `iterations`,
`burnin`, and `thin` settings.

## See also

[`CopRegBayes()`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md)
