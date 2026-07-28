# Summarise a copregbayes fit

[`summary()`](https://rdrr.io/r/base/summary.html) tabulates the
posterior mean, median, standard deviation and credible interval of the
regression coefficients, the residual variance and the endogeneity
correlations; [`print()`](https://rdrr.io/r/base/print.html) on the
result formats that table.

## Usage

``` r
# S3 method for class 'copregbayes'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.copregbayes'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An object of class `"copregbayes"`, as returned by
  [`CopRegBAYES`](https://ashgreat.github.io/endogCopulaBayes/reference/CopRegBayes.md).

- level:

  Credible level for the reported interval, e.g. `0.95` for a 95%
  interval. Defaults to `0.95`.

- ...:

  Not used; present for S3 method consistency.

- x:

  An object of class `"summary.copregbayes"`, as returned by
  `summary.copregbayes`.

- digits:

  Number of significant digits to print.

## Value

For `summary.copregbayes`, an object of class `"summary.copregbayes"`, a
list with `coefficients` (the posterior mean/median/sd/credible-interval
table for the regression coefficients), `sigma2` (the same for the
residual variance), `rho` (the same for the endogeneity correlations
named in `object$rho.names`; correlations among the regressors
themselves are nuisance parameters and are not shown here, though they
remain in `object$Sigma.draws`), `n.other` (how many of those nuisance
correlations there are), `restricted` (names of the correlations held at
exactly zero by the exogeneity restriction), `acceptance`, `horseshoe`,
`n`, `ndraws`, `level`, `iterations`, `burnin`, `thin`, `call` and
`method`.

For `print.summary.copregbayes`, `x` is returned invisibly; called for
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
summary(fit)
#> 
#> Bayesian copula correction (Haschka 2025)
#> 
#> Call:
#> CopRegBAYES(formula = y ~ z | x, data = data.frame(y, z, x), 
#>     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#> 
#> Regression coefficients:
#>             P. Mean  P. Median Sd       2.5%     97.5%   
#> (Intercept)  0.85584  0.82703   0.18691  0.53041  1.12715
#> z            1.75252  1.85036   0.50721  0.67135  2.37521
#> x            0.39360  0.27318   0.39014 -0.04615  1.18508
#> 
#> Structural error variance:
#>        P. Mean P. Median Sd     2.5%   97.5% 
#> sigma2 1.9613  1.8970    0.7394 0.9596 3.7874
#> 
#> Endogeneity: rho(P*, xi*) is the correlation between the normal score 
#>   of an endogenous regressor and that of the structural error.
#>              P. Mean P. Median Sd      2.5%    97.5%  
#> rho(z*, xi*) -0.4564 -0.5433    0.2409 -0.7113  0.1041
#>   1 further correlations among the regressors are estimated jointly and
#>   sit in $Sigma.draws; 1 more are held at zero, which is what makes
#>   the exogenous regressors exogenous.
#> 
#> The two quantile columns are 95% credible intervals: no asymptotic
#>   argument is involved, and the marginal CDFs are estimated jointly rather
#>   than plugged in.
#> 
#> 60 observations. 30 draws kept from 200 iterations (burn-in 50, thinning 5).
#> Prior on the coefficients: horseshoe. Acceptance coefficients 0.835, sigma2 0.865.
```
