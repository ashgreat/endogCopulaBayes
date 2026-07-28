# cran-comments

## Package purpose

endogCopulaBayes implements the Bayesian Gaussian copula endogeneity
correction of Haschka (2025, Oxford Bulletin of Economics and Statistics,
<doi:10.1111/obes.70023>) via a Metropolis-within-Gibbs sampler. The
sampler jointly estimates the regression coefficients, the residual
variance, the copula correlation matrix, and the Dirichlet probability
masses describing the marginal distributions of the regressors, returning
the full MCMC chain alongside a thinned posterior sample with print,
summary, plot, and validity methods.

This is a new submission.

## Test environments

* local: macOS 15 (Darwin 25.5.0), R 4.5.0
* GitHub Actions: ubuntu-latest, R release
* GitHub Actions: macos-latest, R release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.
