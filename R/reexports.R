## =============================================================================
##  Re-exports and generic imports
## =============================================================================
##  validity() is defined as an S3 generic in endogCopula, and copreg-bayes.R
##  registers a validity.copregbayes method for it (importFrom endogCopula
##  validity). That import makes the generic usable from within this package's
##  own code, but does not put it in this package's exported namespace, so a
##  user who has only library(endogCopulaBayes) attached cannot call validity()
##  on a copregbayes fit. Re-exporting the generic here fixes that.

#' @importFrom endogCopula validity
#' @export
endogCopula::validity

## nobs.copregbayes (in copreg-bayes.R) registers an S3 method for the stats
## generic nobs(). Loading a NAMESPACE resolves S3method(nobs, *) generics
## against a built-in table of well-known base/stats/utils/graphics generics
## when the package does not import them explicitly (this is how
## coef.copregbayes, predict.copregbayes, vcov.copregbayes and friends resolve
## without an explicit import). nobs is not on that table, so without this
## importFrom, loadNamespace() fails with "object 'nobs' not found whilst
## loading namespace" -- reproducible even though stats is already listed in
## Imports, since that table lookup is separate from DESCRIPTION.
#' @importFrom stats nobs
NULL
