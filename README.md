# endogCopulaBayes

Bayesian Gaussian copula regression sampler migrated from the original
EndogCopula scripts. This package exposes the MCMC estimator as a standalone
module for future refinement.

## Installation (development)
```r
# install.packages("remotes")
remotes::install_github("ashgreat/endogCopulaBayes")
```

## Usage
```r
library(endogCopulaBayes)
# CopRegBayes(data, iterations = 5000, ...)
```

API details will stabilise as the refactor progresses.
