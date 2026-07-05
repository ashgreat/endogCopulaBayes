#' MCMC Convergence Diagnostics
#'
#' Internal helpers implementing the spectral-density-at-zero variance
#' estimator underlying both the effective sample size (ESS) and the
#' Geweke (1992) convergence diagnostic, mirroring `coda::spectrum0.ar()`
#' in base R so the package does not need a hard dependency on coda.
#' @noRd

# Spectral density at frequency zero via an AR(p) fit with AIC order
# selection, following coda::spectrum0.ar(): regress the column on a linear
# trend and treat it as degenerate (no variability to estimate a spectrum
# for) whenever the residuals are (numerically) constant.
spectrum0_ar <- function(x) {
  n <- length(x)
  z <- seq_len(n)
  fit <- stats::lm(x ~ z)
  if (isTRUE(all.equal(stats::sd(stats::residuals(fit)), 0))) {
    return(list(spec = 0, order = 0L))
  }
  ar_fit <- stats::ar(x, aic = TRUE)
  spec <- ar_fit$var.pred / (1 - sum(ar_fit$ar))^2
  list(spec = spec, order = ar_fit$order)
}

# Effective sample size for a single numeric vector, following
# coda::effectiveSize(): n * var(x) / spectrum0(x), with a zero spectral
# density (degenerate / constant column) mapped to an ESS of zero.
ess_vector <- function(x) {
  n <- length(x)
  if (n < 2 || !is.finite(stats::var(x)) || stats::var(x) == 0) {
    return(0)
  }
  spec <- spectrum0_ar(x)$spec
  if (!is.finite(spec) || spec == 0) {
    return(0)
  }
  n * stats::var(x) / spec
}

# Geweke (1992) z-score for a single numeric vector, following
# coda::geweke.diag() with the default frac1 = 0.1 (start window) and
# frac2 = 0.5 (end window), using integer draw indices 1:n in place of
# coda's mcmc time index (equivalent when thin = 1, as here, since the
# diagnostic already operates on the thinned posterior sample).
geweke_z_vector <- function(x, frac1 = 0.1, frac2 = 0.5) {
  n <- length(x)
  if (n < 4) {
    return(NA_real_)
  }
  start_idx <- 1L
  end_idx <- n
  first_end <- ceiling(start_idx + frac1 * (end_idx - start_idx))
  second_start <- floor(end_idx - frac2 * (end_idx - start_idx))

  first <- x[start_idx:first_end]
  second <- x[second_start:end_idx]

  if (length(first) < 2 || length(second) < 2) {
    return(NA_real_)
  }

  spec1 <- spectrum0_ar(first)$spec / length(first)
  spec2 <- spectrum0_ar(second)$spec / length(second)

  denom <- sqrt(spec1 + spec2)
  if (!is.finite(denom) || denom == 0) {
    return(NA_real_)
  }

  (mean(first) - mean(second)) / denom
}

#' MCMC Convergence Diagnostics
#'
#' Generic function computing per-parameter MCMC convergence diagnostics.
#'
#' @param object A fitted model object.
#' @param ... Passed on to methods.
#' @return A `data.frame` of per-parameter diagnostics.
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' MCMC Convergence Diagnostics for a Bayesian Copula Regression Fit
#'
#' Computes the effective sample size (ESS) and Geweke (1992) convergence
#' z-score for each of the named model parameters (the regression
#' coefficients, residual variance, copula correlations, and hyperprior
#' variances) using the thinned posterior sample stored in `object$posterior`.
#' The `2 * N` unnamed Dirichlet probability masses are not included.
#'
#' Both statistics are implemented in base R using the same spectral-density
#' -at-zero estimator as `coda::spectrum0.ar()`: an AR(p) model is fit to
#' each parameter's draws with the order selected by AIC (via [stats::ar()]),
#' and the resulting `var.pred / (1 - sum(ar))^2` is used as the long-run
#' variance in both the ESS ratio and the Geweke z-score denominator. A
#' column whose posterior draws are (numerically) constant has an
#' undefined spectral density; its ESS is reported as `0` and its Geweke
#' z-score as `NA`, rather than raising an error.
#'
#' @param object An object of class `endog_copula_bayes`.
#' @param ... Ignored.
#' @return A `data.frame` with one row per named parameter (row names give
#'   the parameter names, respecting variable labels for a formula fit) and
#'   columns:
#'   \describe{
#'     \item{`ess`}{Effective sample size.}
#'     \item{`geweke_z`}{Geweke (1992) convergence z-score comparing the
#'       first 10% and last 50% of the thinned posterior sample.}
#'     \item{`n_draws`}{Number of thinned posterior draws used.}
#'   }
#' @seealso [CopRegBayes()], [plot.endog_copula_bayes()]
#' @examples
#' set.seed(1)
#' n <- 60
#' z <- rlnorm(n)
#' x <- rnorm(n)
#' dat <- data.frame(y = 2 - 4 * z + 6 * x + rnorm(n, sd = sqrt(2)),
#'                   z = z, x = x)
#' fit <- CopRegBayes(dat, iterations = 100, burnin = 20, thin = 2, seed = 42)
#' diagnostics(fit)
#' @importFrom stats ar lm residuals sd var
#' @export
diagnostics.endog_copula_bayes <- function(object, ...) {
  draws <- object$posterior[, object$parameters, drop = FALSE]
  n_draws <- nrow(draws)

  ess <- apply(draws, 2, ess_vector)
  geweke_z <- apply(draws, 2, geweke_z_vector)

  out <- data.frame(
    ess = as.numeric(ess),
    geweke_z = as.numeric(geweke_z),
    n_draws = rep(n_draws, ncol(draws)),
    row.names = label_parameters(colnames(draws), object$variables)
  )
  out
}

#' Plot a Bayesian Copula Regression Fit
#'
#' Base-graphics traceplots and posterior density plots for the named
#' parameters of a fitted `endog_copula_bayes` object (the `2 * N` unnamed
#' Dirichlet probability masses are excluded).
#'
#' @param x An object of class `endog_copula_bayes`.
#' @param params Character vector of parameter names to plot (after
#'   variable-label substitution for a formula fit); `NULL` (default) plots
#'   all named parameters in `x$parameters`.
#' @param which Character vector selecting which diagnostic plots to draw:
#'   `"trace"` for traceplots, `"density"` for posterior density plots, or
#'   both (the default).
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @seealso [CopRegBayes()], [diagnostics()]
#' @examples
#' set.seed(1)
#' n <- 60
#' z <- rlnorm(n)
#' x <- rnorm(n)
#' dat <- data.frame(y = 2 - 4 * z + 6 * x + rnorm(n, sd = sqrt(2)),
#'                   z = z, x = x)
#' fit <- CopRegBayes(dat, iterations = 100, burnin = 20, thin = 2, seed = 42)
#' plot(fit, params = c("beta_z", "beta_x"))
#' @importFrom graphics par plot
#' @importFrom stats density
#' @export
plot.endog_copula_bayes <- function(x, params = NULL,
                                    which = c("trace", "density"), ...) {
  which <- match.arg(which, several.ok = TRUE)
  all_labels <- label_parameters(x$parameters, x$variables)
  if (is.null(params)) {
    params <- all_labels
  }
  missing_params <- setdiff(params, all_labels)
  if (length(missing_params) > 0) {
    stop(sprintf("Unknown parameter(s): %s", paste(missing_params, collapse = ", ")),
         call. = FALSE)
  }

  draws <- x$posterior[, x$parameters, drop = FALSE]
  colnames(draws) <- all_labels

  n_plot_cols <- length(which)
  n_params <- length(params)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # One row per parameter, one column per requested plot type ("trace" and/or
  # "density"), capped at four rows per screen so large parameter sets still
  # page reasonably via the graphics device's own paging.
  n_row <- min(n_params, 4L)
  graphics::par(mfrow = c(n_row, n_plot_cols), mar = c(4, 4, 2, 1))

  for (p in params) {
    v <- draws[, p]
    if ("trace" %in% which) {
      graphics::plot(v, type = "l", main = paste("Trace:", p),
                     xlab = "Iteration", ylab = p)
    }
    if ("density" %in% which) {
      if (isTRUE(all.equal(stats::sd(v), 0))) {
        graphics::plot(v, rep(0, length(v)), type = "n",
                       main = paste("Density:", p), xlab = p, ylab = "Density")
      } else {
        d <- stats::density(v)
        graphics::plot(d, main = paste("Density:", p), xlab = p)
      }
    }
  }

  invisible(x)
}
