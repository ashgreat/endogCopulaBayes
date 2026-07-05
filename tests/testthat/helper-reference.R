# Locate the reference research code (ground truth) used for equivalence
# tests. Set ENDOG_REF_DIR to point at the
# "Copula-based-endogeneity-corrections-main" folder, otherwise the fallback
# assumes the package lives next to it inside the EndogCopula project.
endog_ref_dir <- function() {
  dir <- Sys.getenv("ENDOG_REF_DIR", unset = "")
  if (nzchar(dir)) {
    return(dir)
  }
  testthat::test_path("..", "..", "..", "Copula-based-endogeneity-corrections-main")
}

skip_if_no_reference <- function() {
  if (!dir.exists(endog_ref_dir())) {
    testthat::skip("Reference implementation not available (set ENDOG_REF_DIR)")
  }
}

# Parse a reference .R file and evaluate ONLY top-level assignments whose
# right-hand side is a function definition, into a fresh environment whose
# parent is the global environment. This keeps the top-level demo code in the
# reference scripts (simulations, plots, long MCMC runs) from executing.
load_reference_functions <- function(filename) {
  path <- file.path(endog_ref_dir(), filename)
  exprs <- parse(path)
  env <- new.env(parent = globalenv())
  for (expr in exprs) {
    if (is.call(expr) &&
        (identical(expr[[1]], as.name("<-")) || identical(expr[[1]], as.name("="))) &&
        length(expr) == 3L &&
        is.call(expr[[3]]) &&
        identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, envir = env)
    }
  }
  env
}
