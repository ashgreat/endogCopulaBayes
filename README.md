# endogCopulaBayes

<!-- badges: start -->
[![R-CMD-check](https://github.com/ashgreat/endogCopulaBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ashgreat/endogCopulaBayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Package website: <https://ashgreat.github.io/endogCopulaBayes/>

endogCopulaBayes implements the Bayesian Gaussian copula endogeneity
correction of Haschka (2025, Oxford Bulletin of Economics and Statistics).
The estimator handles an endogenous regressor without instruments. It
models the joint distribution of the regressors and the regression error
with a Gaussian copula, and draws every unknown quantity with a
Metropolis-within-Gibbs MCMC sampler: the regression coefficients, the
residual variance, the copula correlation matrix, and the Dirichlet
probability masses describing the marginal distribution of each regressor.

The package builds on endogCopula, which supplies the shared model parsing
and diagnostic helpers used across the whole toolbox.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ashgreat/endogCopulaBayes")
```

## Usage

The exported function is `CopRegBAYES()`. The model is specified with the
same two part formula used across the toolbox, `y ~ endogenous | exogenous`,
with the endogenous regressor before the bar and the exogenous regressor
after it.

There is no `cdf` argument and no `ties` argument. The CDF of each
regressor is drawn by the sampler rather than plugged in, so ties are the
normal case and a binary regressor simply gets two probability masses.
There is also no `nboots` argument. Inference is posterior, controlled by
`iterations`, `burnin` and `thin` instead. Supplying `cdf`, `ties` or
`nboots` raises an error.

Exogeneity of the regressor after the bar is imposed by holding the
matching entries of the copula correlation matrix at exactly zero.

```r
library(endogCopulaBayes)

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
fit
```

This is the real output of the example above, run with `iterations = 3000`
to keep it fast. Use many more iterations in practice. The default is
102,000.

```
Bayesian copula correction (Haschka 2025)

Call:
CopRegBAYES(formula = y ~ z | x, data = dat, iterations = 3000, 
    burnin = 500, thin = 10)

Posterior means:
(Intercept)            z            x  
    0.43217      2.55402     -0.04718  

250 draws from 3,000 iterations.
```

`summary(fit)` adds posterior standard deviations, medians, and 95 percent
credible intervals for the coefficients, the residual variance, and the
endogeneity correlation:

```
Regression coefficients:
            P. Mean  P. Median Sd       2.5%     97.5%   
(Intercept)  0.43217  0.45108   0.26719 -0.07975  0.92392
z            2.55402  2.53253   0.24978  2.06653  3.05664
x           -0.04718 -0.02867   0.23309 -0.56803  0.41392

Structural error variance:
       P. Mean P. Median Sd     2.5%   97.5% 
sigma2 0.7518  0.7418    0.1363 0.5715 1.0912

Endogeneity: rho(P*, xi*) is the correlation between the normal score 
  of an endogenous regressor and that of the structural error.
             P. Mean  P. Median Sd       2.5%     97.5%   
rho(z*, xi*)  0.04049  0.05296   0.22706 -0.42343  0.42077
```

## Checking convergence and identification

`validity(fit)`, an S3 method re-exported from endogCopula, reports three
things: non-normality of each endogenous regressor, since identification
needs that non-normality; the posterior of the endogeneity correlations, so
an interval that covers zero says the data carry no evidence of
endogeneity; and convergence of the chain, through Geweke's statistic, the
effective sample size, the lag-one autocorrelation, and the
Metropolis-Hastings acceptance rates. The effective sample size is bounded
above by the number of retained draws, and is `NA` for a chain that does
not move. The Gelman-Rubin statistic across several chains from dispersed
starting values is available on request by passing `chains = TRUE`.

The package deliberately differs from the reference implementation it
reproduces in three small places, all inside `validity()` and none of them
affecting any estimate; see "Differences from the reference implementation"
in [NEWS.md](NEWS.md) for the list.

## Related packages

- endogCopula (<https://github.com/ashgreat/endogCopula>): cross-sectional
  Gaussian copula corrections with a shared formula interface.
- endogCopulaPanel (<https://github.com/ashgreat/endogCopulaPanel>): the
  fixed effects panel maximum likelihood estimator of Haschka (2022).

## References

- Haschka, R. E. (2025). Bayesian inference for joint estimation models
  using copulas to handle endogenous regressors. Oxford Bulletin of
  Economics and Statistics. doi:10.1111/obes.70023
- Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
  estimation using copulas. Marketing Science, 31(4), 567-586.

## License

MIT + file LICENSE. See [LICENSE.md](LICENSE.md) for the full text.
