# Print a copregbayes fit

Prints the call and the posterior mean of the regression coefficients
for a
[`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md)
fit.

## Usage

``` r
# S3 method for class 'copregbayes'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"copregbayes"`, as returned by
  [`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md).

- digits:

  Number of significant digits to print.

- ...:

  Not used; present for S3 method consistency.

## Value

`x`, invisibly. Called for its side effect of printing.

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
print(fit)
#> 
#> Bayesian copula correction (Haschka 2025)
#> 
#> Call:
#> CopRegBAYES(formula = y ~ z | x, data = data.frame(y, z, x), 
#>     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#> 
#> Posterior means:
#> (Intercept)            z            x  
#>      0.8558       1.7525       0.3936  
#> 
#> 30 draws from 200 iterations.
```
