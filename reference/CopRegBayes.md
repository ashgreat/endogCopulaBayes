# Bayesian Gaussian copula endogeneity correction

Joint Bayesian estimation of a linear regression together with a
Gaussian copula correction for endogenous regressors (Haschka 2025).
Unlike the other estimators in the toolbox, the marginal CDFs of the
regressors are not plugged in from a first stage: each regressor gets
probability masses on its uniquely observed values (Section 3.1 of the
paper), and those masses are drawn by MCMC together with the regression
coefficients, the residual variance and the full copula correlation
matrix. There is therefore no first-stage plug-in uncertainty to
bootstrap back in; the posterior carries all of it.

There is no `cdf` argument and no `ties` argument: the CDF of each
regressor is a parameter of the model rather than a choice the user
makes, and ties are the normal case rather than a problem to work around
– the probability masses sit on unique values, so a binary regressor
simply has two of them. There is likewise no `nboots` argument:
inference here is posterior rather than bootstrap, governed instead by
`iterations`, `burnin` and `thin`. Supplying any of `cdf`, `ties` or
`nboots` is an error.

Exogeneity of the regressors listed after `|` is imposed by holding the
corresponding entries of the copula covariance matrix at exactly zero.
This is done by a block-recursive reparameterisation of the covariance
matrix (endogenous regressors regressed on the exogenous ones and on the
structural error, all independent by construction) rather than by an
inverse Wishart draw, so the zeros hold exactly rather than only
approximately.

## Usage

``` r
CopRegBAYES(
  formula,
  data,
  iterations = 102000,
  burnin = 2000,
  thin = 100,
  horseshoe = TRUE,
  prior.args = list(),
  start = NULL,
  subset = NULL,
  contrasts = NULL,
  verbose = interactive(),
  cdf,
  ties,
  nboots
)
```

## Arguments

- formula:

  Two-part formula, `y ~ endogenous | exogenous`, as everywhere else in
  the toolbox: factors, interactions, transformations and `-1` all work,
  and the position of a term around `|` decides whether it is treated as
  endogenous (and so gets a copula term) or exogenous. Every
  non-intercept column of the design matrix enters the copula, so a
  four-level factor on the exogenous side contributes three columns with
  two unique values each.

- data:

  A `data.frame` holding the variables in `formula`.

- iterations:

  Total number of MCMC iterations. Defaults to 102,000, the setting used
  in the paper's simulations.

- burnin:

  Number of iterations discarded from the start of the chain. Defaults
  to 2,000.

- thin:

  Keep every `thin`-th iteration after burn-in. Defaults to 100, which
  together with the defaults above leaves 1,000 draws. Draws are thinned
  as they are produced rather than after the fact, so only the retained
  iterates are ever held in memory.

- horseshoe:

  If `TRUE` (the default), a horseshoe prior is placed on every
  regression coefficient, including the intercept. If `FALSE`, each
  coefficient instead gets an independent \\N(0, \phi^2)\\ prior with
  \\\phi^2\\ drawn by a conjugate Gibbs step from an inverse Gamma full
  conditional.

- prior.args:

  A list overriding any of the hyperparameters, all of which default to
  the paper's own settings: `a`, `b` for the inverse Gamma prior on the
  residual variance; `Psi1`, `nu1` for the exogenous-regressor
  covariance block; `Psi2`, `nu2` for the copula error-covariance block;
  `a_e`, `b_e` for the structural-error variance in that
  reparameterisation; `M0`, `V0` for the
  regression-on-exogenous-and-error prior; and `a_phi`, `b_phi` for the
  normal-prior variance used when `horseshoe = FALSE`.

- start:

  A list overriding any of the starting values (any of `coefficients`,
  `sigma2`, `Sigma`, `lambda`). The defaults follow Web Appendix D of
  the paper: OLS for the coefficients and the residual variance, the
  identity matrix for `Sigma`, and one draw from a Dir(1) distribution
  per regressor for the probability masses.

- subset:

  An optional logical or index vector selecting the rows of `data` to
  use, as in [`lm()`](https://rdrr.io/r/stats/lm.html).

- contrasts:

  An optional list passed to
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) for the
  coding of factors.

- verbose:

  If `TRUE`, print a starting message and show a text progress bar while
  sampling. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- cdf:

  Not a valid argument for this estimator: the marginal CDF is a
  parameter drawn by the sampler rather than a plug-in choice. Supplying
  it raises an error.

- ties:

  Not a valid argument for this estimator, for the same reason as `cdf`.
  Supplying it raises an error.

- nboots:

  Not a valid argument for this estimator: inference is posterior rather
  than bootstrap. Use `iterations`, `burnin` and `thin` instead.
  Supplying it raises an error.

## Value

An object of class `"copregbayes"`, a list including

- `coefficients`:

  Posterior mean of the regression coefficients (what
  [`coef()`](https://rdrr.io/r/stats/coef.html) returns).

- `posterior.median`:

  Posterior median of the coefficients.

- `std.error`:

  Posterior standard deviation of the coefficients.

- `sigma2`:

  Posterior mean of the residual variance.

- `vcov`:

  Posterior covariance matrix of the coefficients (what
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns).

- `draws`:

  Matrix of thinned post-burn-in draws (one row per retained iteration)
  of the coefficients, `sigma2`, and the free (non-zero-restricted)
  entries of the copula correlation matrix; the basis for
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html).

- `coefficient.draws`, `sigma2.draws`:

  The coefficient and residual-variance draws separately.

- `Sigma.draws`:

  Matrix of the thinned draws of every lower-triangular entry of the
  copula correlation matrix, free and zero-restricted alike, named
  `"rho(a*, b*)"`.

- `Sigma.free`:

  Logical vector marking which columns of `Sigma.draws` are free; the
  rest are held at exactly zero by the exogeneity restriction.

- `lambda.draws`:

  A list, one matrix per endogenous or exogenous regressor entering the
  copula, of the thinned draws of its Dirichlet probability masses (one
  column per unique observed value).

- `margins`:

  The unique values, grouping and counts underlying each regressor's
  masses (see
  [`plot.copregbayes`](https://ashgreat.github.io/endogCopulaBayes/reference/plot.copregbayes.md),
  `type = "cdf"`).

- `acceptance`:

  Metropolis-Hastings acceptance rates for the regression-coefficient
  and log-residual-variance updates.

- `fitted.values`, `residuals`:

  Fitted values and residuals at the posterior mean coefficients.

- `endo.names`, `exog.names`:

  Names of the regressors entering the copula on the endogenous and
  exogenous side.

- `rho.names`:

  Names, in `draws` and `Sigma.draws`, of the correlations between each
  endogenous regressor's normal score and that of the structural error –
  the endogeneity correlations.

- `horseshoe`, `prior.args`, `iterations`, `burnin`, `thin`, `ndraws`:

  The settings the model was fit with, and the number of retained draws.

- `method`, `call`:

  A description string and the matched call.

Also carries `y`, `X`, `n`, `formula`, `terms`, `mf`, `xlevels`,
`contrasts` and `na.action` for use by
[`predict()`](https://rdrr.io/r/stats/predict.html) and related methods.

## References

Haschka, R. E. (2025). Bayesian inference for joint estimation models
using copulas to handle endogenous regressors. *Oxford Bulletin of
Economics and Statistics*.
[doi:10.1111/obes.70023](https://doi.org/10.1111/obes.70023)

## Examples

``` r
# \donttest{
set.seed(1)
n   <- 150
x   <- rnorm(n)
rho <- 0.6
eps <- matrix(rnorm(2 * n), n, 2)
eps[, 2] <- rho * eps[, 1] + sqrt(1 - rho^2) * eps[, 2]
z   <- 1 + x + eps[, 1]              # endogenous, correlated with e below
e   <- eps[, 2]
y   <- 1 + 2 * z + 0.5 * x + e
dat <- data.frame(y = y, z = z, x = x)

fit <- CopRegBAYES(y ~ z | x, data = dat,
                    iterations = 3000, burnin = 500, thin = 10)
print(fit)
#> 
#> Bayesian copula correction (Haschka 2025)
#> 
#> Call:
#> CopRegBAYES(formula = y ~ z | x, data = dat, iterations = 3000, 
#>     burnin = 500, thin = 10)
#> 
#> Posterior means:
#> (Intercept)            z            x  
#>     0.43217      2.55402     -0.04718  
#> 
#> 250 draws from 3,000 iterations.
summary(fit)
#> 
#> Bayesian copula correction (Haschka 2025)
#> 
#> Call:
#> CopRegBAYES(formula = y ~ z | x, data = dat, iterations = 3000, 
#>     burnin = 500, thin = 10)
#> 
#> Regression coefficients:
#>             P. Mean  P. Median Sd       2.5%     97.5%   
#> (Intercept)  0.43217  0.45108   0.26719 -0.07975  0.92392
#> z            2.55402  2.53253   0.24978  2.06653  3.05664
#> x           -0.04718 -0.02867   0.23309 -0.56803  0.41392
#> 
#> Structural error variance:
#>        P. Mean P. Median Sd     2.5%   97.5% 
#> sigma2 0.7518  0.7418    0.1363 0.5715 1.0912
#> 
#> Endogeneity: rho(P*, xi*) is the correlation between the normal score 
#>   of an endogenous regressor and that of the structural error.
#>              P. Mean  P. Median Sd       2.5%     97.5%   
#> rho(z*, xi*)  0.04049  0.05296   0.22706 -0.42343  0.42077
#>   1 further correlations among the regressors are estimated jointly and
#>   sit in $Sigma.draws; 1 more are held at zero, which is what makes
#>   the exogenous regressors exogenous.
#> 
#> The two quantile columns are 95% credible intervals: no asymptotic
#>   argument is involved, and the marginal CDFs are estimated jointly rather
#>   than plugged in.
#> 
#> 150 observations. 250 draws kept from 3,000 iterations (burn-in 500, thinning 10).
#> Prior on the coefficients: horseshoe. Acceptance coefficients 0.806, sigma2 0.939.
# }
```
