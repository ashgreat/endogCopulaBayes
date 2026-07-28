## =============================================================================
##  Copreg_bayes.R
##
##  Bayesian joint estimation of the Gaussian copula endogeneity correction,
##  Haschka (2025), Oxford Bulletin of Economics and Statistics.
##
##  Requires copreg-core.R to be sourced first.  Base R only.
## =============================================================================
##  Ported into the endogCopulaBayes package: the ".copreg_model" and
##  ".validity_nonnormality" helpers below now come from endogCopula (Imports)
##  and are called as endogCopula::.copreg_model() / endogCopula::.validity_nonnormality(),
##  so the "source copreg-core.R first" guard that used to sit here has been
##  removed -- the package DESCRIPTION's Imports: endogCopula plays that role.


## -----------------------------------------------------------------------------
##  What this estimator does, and why it is its own class
## -----------------------------------------------------------------------------
##  Every other estimator in the toolbox computes marginal CDFs first, builds a
##  copula term from them, runs a regression and bootstraps the whole thing.
##  This one does none of that.  The CDFs of the regressors are not plugged in
##  but are unknown parameters: each regressor gets probability masses on its
##  uniquely observed values (Section 3.1), and those masses are drawn along
##  with the regression coefficients, the error variance and the full copula
##  correlation matrix in one MCMC run.  There is no first stage and nothing is
##  estimated a priori, so there are no plug-in estimates whose uncertainty
##  would have to be bootstrapped back in; the posterior carries all of it.
##
##  Consequences for the interface:
##
##    * cdf and ties do not exist.  The CDF is a parameter, and ties are the
##      normal case -- the masses sit on unique values, so a binary regressor
##      simply has two of them (the application in the paper has 121 weekly
##      dummies of exactly that kind).
##    * nboots is replaced by iterations, burnin and thin.
##    * rho IS reported here, unlike in 2sCOPE-np: the copula correlation
##      matrix is a parameter of the model, so its posterior comes for free.
##    * confint() returns credible intervals, not confidence intervals, and no
##      asymptotic argument is involved anywhere.
##
##  Ordering convention throughout, taken from the paper: z (endogenous), then
##  x (exogenous), then e.  The intercept is a regressor but not part of the
##  copula -- it has one unique value, so it carries no distributional
##  information (Web Appendix E does the same).


## =============================================================================
##  1.  Base R replacements for invgamma, mvtnorm, copula and LaplacesDemon
## =============================================================================
##  The reference implementation pulls in four packages for draws and densities
##  that are a few lines each.  Writing them out keeps the dependency list at
##  minimal, and the copula density in particular gets much faster: the
##  Gaussian copula log density is given in closed form in Equation (3), so
##  there is no reason to route N observations through a generic dCopula().

.bayes_rinvgamma <- function(n, shape, rate) 1 / stats::rgamma(n, shape, rate = rate)

## Wishart by the Bartlett decomposition.  V = LL' lower Cholesky, A lower
## triangular with chi distributed diagonal, then LA(LA)' ~ W(nu, V).
.bayes_rwish <- function(nu, V) {
  p <- ncol(V)
  L <- t(chol(V))
  A <- matrix(0, p, p)
  diag(A) <- sqrt(stats::rchisq(p, nu - seq_len(p) + 1))
  if (p > 1L) A[lower.tri(A)] <- stats::rnorm(p * (p - 1L) / 2L)
  tcrossprod(L %*% A)
}

## X ~ IW(nu, S)  <=>  X^{-1} ~ W(nu, S^{-1}).  Scale parameterisation, i.e.
## E(X) = S / (nu - p - 1), which is what the paper's Equation (8) uses.
.bayes_riwish <- function(nu, S) {
  Si <- chol2inv(chol(S))
  chol2inv(chol(.bayes_rwish(nu, Si)))
}

## Multivariate normal draw and log density, both from one Cholesky factor.
.bayes_rmvn <- function(mu, R) as.vector(mu + crossprod(R, stats::rnorm(ncol(R))))
.bayes_dmvn <- function(x, mu, R) {
  z <- backsolve(R, x - mu, transpose = TRUE)
  -0.5 * length(x) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z * z)
}


## =============================================================================
##  2.  Margins: probability masses on the unique values of a regressor
## =============================================================================
##  Section 3.1.  For a generic regressor w the density is treated as a
##  nonparametric function putting mass lambda_j on the j-th smallest unique
##  value, and the margin entering the copula is the cumulative sum.
##
##  Two details that the reference implementations disagree on, resolved here
##  in favour of the version that survives ties:
##
##    * cells are the unique VALUES, not the observations.  With m unique
##      values out of N observations the mass vector has length m, and the
##      counts v enter the Dirichlet full conditional.  Indexing by rank()
##      instead breaks as soon as two observations share a value, which for
##      the dummies of the empirical application is every observation.
##    * the margin is evaluated at the midpoint of the mass, cumsum(l) - l/2,
##      not at its upper edge.  For a continuous regressor with all masses of
##      order 1/N the difference is negligible; for a binary one it is the
##      difference between a sensible normal score and Phi^{-1}(1) = Inf.

.bayes_margin <- function(v) {
  uv  <- sort(unique(v))
  grp <- match(v, uv)
  m   <- length(uv)
  list(uv = uv, grp = grp, m = m, cnt = tabulate(grp, nbins = m),
       f = factor(grp, levels = seq_len(m)))
}

## normal scores of the margin implied by the current masses
.bayes_xi <- function(lambda, mg) {
  Fmid <- (cumsum(lambda) - lambda / 2) * (mg$m / (mg$m + 0.01))
  stats::qnorm(Fmid[mg$grp])
}

## Gibbs draw of lambda, Web Appendix B following Ng et al. (2011): a copula
## correlated uniform per observation, pushed through the Gamma(1,1) quantile
## function, aggregated within value cells, plus one Gamma(1,1) prior
## pseudo-count per cell for the Dir(1) prior.  Normalising gives a draw from
## Dir(1 + v).
.bayes_draw_lambda <- function(u, mg) {
  g <- stats::qgamma(u, shape = 1, rate = 1)
  G <- as.numeric(tapply(g, mg$f, sum))
  G[is.na(G)] <- 0
  G <- G + stats::rgamma(mg$m, shape = 1, rate = 1)
  G / sum(G)
}


## =============================================================================
##  3.  Log posterior, score and working weight
## =============================================================================
##  Only the parts that depend on (alpha, beta, delta) and sigma^2 are needed:
##  the two Metropolis-Hastings steps condition on Sigma and lambda, so the
##  categorical likelihood terms and the Dirichlet priors of Equation (4) are
##  constant and cancel in the acceptance ratio.  Dropping them is not an
##  approximation, it is the same ratio computed with fewer operations.
##
##  Writing A = Sigma^{-1} - I and splitting xi = (xi_zx, xi_e), the copula
##  exponent is
##
##      xi' A xi  =  xi_zx' A_zx,zx xi_zx  +  2 xi_e (A_zx,e' xi_zx)  +  A_ee xi_e^2
##
##  whose first term is constant within the step.  Precomputing q = A_zx,e' xi_zx
##  once turns each evaluation into O(N) instead of O(N d^2).

.bayes_lpost <- function(cf, s2, y, X, q, aee, pr, horseshoe, phi2) {
  if (s2 <= 0) return(-Inf)
  e   <- y - as.vector(X %*% cf)
  xe  <- e / sqrt(s2)
  q2  <- sum(xe * xe)
  lp  <- -0.5 * (2 * sum(q * xe) + aee * q2) -
    0.5 * length(e) * log(s2) - 0.5 * q2
  ## inverse Gamma prior on sigma^2, Equation (7)
  lp <- lp - (pr$a + 1) * log(s2) - pr$b / s2
  ## prior on the regression coefficients, intercept included
  lp + if (horseshoe) {
    ## The horseshoe of Equation (5) with phi, tau and iota marginalised out.
    ## That marginal has no closed form; Carvalho, Polson & Scott (2010,
    ## Thm. 1) bracket it by
    ##
    ##   (K/2) log(1 + 4/gamma^2)  <  p(gamma)  <  K log(1 + 2/gamma^2),
    ##
    ## and the upper bound is the usual working form.  Both ends are right:
    ## p behaves like 1/gamma^2 far out, which is the Cauchy tail Equation (5)
    ## puts there, and diverges at zero, which is what does the shrinking.
    ##
    ## Reading Equation (6) literally, as p(gamma) proportional to 1/|gamma|,
    ## does not work.  That density is not integrable at zero, so the spike
    ## there is infinitely deep rather than merely tall.  The proposal moves
    ## the whole coefficient vector at once, so as soon as one coefficient
    ## lands near zero the prior term explodes, the move is accepted, and
    ## every later move away from it is rejected: acceptance collapses to a
    ## fraction of a percent and the trace becomes a flat line.  The log-log
    ## divergence below shrinks just as hard and is integrable.
    sum(log(log1p(2 / pmax(cf * cf, 1e-300))))
  } else {
    sum(stats::dnorm(cf, 0, sqrt(phi2), log = TRUE))
  }
}

## First and second derivative of the log full conditional of tau = log sigma^2.
## Web Appendix C, plus the Jacobian of the transformation: the target in tau
## is p(sigma^2) * sigma^2, so the prior contributes -(a+1) + b/sigma^2 + 1,
## which is the -a + b/sigma^2 below.  The acceptance ratio must carry the
## matching (tau_new - tau_cur); the reference implementation omits it there
## while using this score, which biases sigma^2.
.bayes_tau_derivs <- function(cf, s2, y, X, q, aee, pr) {
  e   <- y - as.vector(X %*% cf)
  xe  <- e / sqrt(s2)
  qxe <- sum(q * xe)
  q2  <- sum(xe * xe)
  list(score = 0.5 * qxe + 0.5 * (aee + 1) * q2 - 0.5 * length(e) -
         pr$a + pr$b / s2,
       f2 = -0.25 * qxe - 0.5 * (aee + 1) * q2 - pr$b / s2)
}


## =============================================================================
##  4.  Block-recursive Gibbs step for the copula covariance matrix
## =============================================================================
##  Exogeneity of x means Cov(xi_x, xi_e) = 0.  Web Appendix A of the paper
##  draws W from a plain inverse Wishart full conditional, which does not
##  respect those zeros: the inverse Wishart family has no way to hold entries
##  of a covariance matrix at zero.
##
##  Instead W is reparameterised hierarchically,
##
##      xi_z = B_x xi_x + beta_e xi_e + eps,
##      xi_x ~ MN(0, Sigma_xx),  xi_e ~ N(0, sigma_e^2),  eps ~ MN(0, Omega),
##
##  with the three components independent.  The implied W has the zero blocks
##  by construction, and the map from (B_x, beta_e, Omega, Sigma_xx, sigma_e^2)
##  to W is bijective given the restrictions, so nothing is lost.  With
##  C = [B_x | beta_e] the four conjugate full conditionals are
##
##      Sigma_xx | .  ~  IW(nu1 + N, Psi1 + S_xx)
##      sigma_e^2 | . ~  IG(a + N/2, b + S_ee/2)
##      C | Omega, .  ~  MN(Mbar, Omega, Vbar)
##      Omega | C, .  ~  IW(nu2 + N, Psi2 + S_epseps)
##
##  drawn in that order, C using the Omega of the previous cycle.  Sigma is
##  then the correlation matrix of the reconstructed W.

.bayes_draw_W <- function(xi_z, xi_x, xi_e, Omega, pr) {
  
  N <- length(xi_e)
  K <- ncol(xi_z)
  L <- if (is.null(xi_x)) 0L else ncol(xi_x)
  
  Xt <- if (L > 0L) cbind(xi_x, xi_e) else matrix(xi_e, N, 1L)
  St <- crossprod(Xt)
  Sz <- crossprod(xi_z, Xt)
  
  Sxx <- if (L > 0L) .bayes_riwish(pr$nu1 + N, pr$Psi1 + crossprod(xi_x)) else NULL
  se2 <- .bayes_rinvgamma(1, pr$a_e + N / 2, pr$b_e + sum(xi_e^2) / 2)
  
  Vb <- chol2inv(chol(pr$V0inv + St))
  Vb <- (Vb + t(Vb)) / 2
  Mb <- (pr$M0 %*% pr$V0inv + Sz) %*% Vb
  C  <- Mb + t(chol(Omega)) %*% matrix(stats::rnorm(K * (L + 1L)), K, L + 1L) %*%
    chol(Vb)
  
  R     <- xi_z - Xt %*% t(C)
  Omega <- .bayes_riwish(pr$nu2 + N, pr$Psi2 + crossprod(R))
  
  Bx <- C[, seq_len(L), drop = FALSE]
  be <- C[, L + 1L, drop = FALSE]
  
  W11 <- Omega + se2 * tcrossprod(be)
  if (L > 0L) W11 <- W11 + Bx %*% Sxx %*% t(Bx)
  W13 <- se2 * be
  W <- if (L > 0L) {
    W12 <- Bx %*% Sxx
    rbind(cbind(W11, W12, W13),
          cbind(t(W12), Sxx, matrix(0, L, 1L)),
          cbind(t(W13), matrix(0, 1L, L), se2))
  } else {
    rbind(cbind(W11, W13), cbind(t(W13), se2))
  }
  
  s <- sqrt(diag(W))
  list(Sigma = W / outer(s, s), Omega = Omega)
}


## =============================================================================
##  5.  The sampler
## =============================================================================
##  Algorithm 1 of Web Appendix D: one Metropolis-Hastings step for the
##  regression coefficients using the iteratively weighted least squares
##  proposal of Gamerman (1997), one for log sigma^2 using a Newton proposal,
##  then Gibbs steps for W and for the masses.
##
##  Draws are thinned as they are produced rather than afterwards.  The
##  reference implementation keeps every iterate and thins at the end, which
##  for the empirical application of the paper -- 1,002,000 iterations and some
##  2,600 mass columns -- would need around 21 GB.  Keeping only the retained
##  iterates brings that down to a few MB.

.bayes_sampler <- function(y, X, cop, mgs, iterations, burnin, thin,
                           horseshoe, pr, start, verbose) {
  
  N  <- length(y)
  p  <- ncol(X)
  K  <- length(cop$endo)
  L  <- length(cop$exog)
  d  <- K + L + 1L
  ci <- seq_len(d - 1L)                       # copula columns other than e
  
  ndraw <- length(seq.int(burnin + 1L, iterations, by = thin))
  out <- list(
    coefficients = matrix(NA_real_, ndraw, p, dimnames = list(NULL, colnames(X))),
    sigma2       = numeric(ndraw),
    Sigma        = matrix(NA_real_, ndraw, d * (d - 1L) / 2L),
    lambda       = lapply(mgs, function(m) matrix(NA_real_, ndraw, m$m)),
    accept       = c(coefficients = 0, sigma2 = 0))
  names(out$lambda) <- names(mgs)
  
  cf   <- start$coefficients
  s2   <- start$sigma2
  Sig  <- start$Sigma
  lam  <- start$lambda
  phi2 <- start$phi2
  Om   <- diag(K)
  
  ## The same progress bar the bootstrap uses.  It cannot be redrawn on every
  ## iteration here: at 102,000 of them that is 102,000 writes to the console
  ## for a bar some 190 characters wide, which costs more than it shows.
  ## Redrawing at most a thousand times is beyond the resolution of the bar
  ## and free.
  keep <- 0L
  nxt  <- burnin + 1L
  step <- max(1L, iterations %/% 1000L)
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  for (it in seq_len(iterations)) {
    
    ## --- normal scores of the margins under the current masses --------------
    Xi <- matrix(0, N, d)
    for (j in seq_len(d - 1L)) Xi[, j] <- .bayes_xi(lam[[j]], mgs[[j]])
    
    Sinv <- chol2inv(chol(Sig))
    A    <- Sinv - diag(d)
    q    <- as.vector(Xi[, ci, drop = FALSE] %*% A[ci, d])
    aee  <- A[d, d]
    
    ## --- Metropolis-Hastings for the regression coefficients ---------------
    ## IWLS proposal.  The working weight is the exact observed information
    ## per observation, (Sigma^{-1})_{dd} / sigma^2, so the proposal covariance
    ## does not depend on the coefficients and reusing it in the reverse
    ## density is exact rather than approximate.
    e  <- y - as.vector(X %*% cf)
    sg <- sqrt(s2)
    nu <- (q + aee * e / sg) / sg + e / s2
    Rp <- chol(chol2inv(chol((Sinv[d, d] / s2) * crossprod(X))))
    
    mu  <- cf + as.vector(crossprod(Rp, Rp %*% crossprod(X, nu)))
    cfp <- .bayes_rmvn(mu, Rp)
    
    ep  <- y - as.vector(X %*% cfp)
    nup <- (q + aee * ep / sg) / sg + ep / s2
    mup <- cfp + as.vector(crossprod(Rp, Rp %*% crossprod(X, nup)))
    
    la <- .bayes_lpost(cfp, s2, y, X, q, aee, pr, horseshoe, phi2) -
      .bayes_lpost(cf, s2, y, X, q, aee, pr, horseshoe, phi2) +
      .bayes_dmvn(cf, mup, Rp) - .bayes_dmvn(cfp, mu, Rp)
    if (is.finite(la) && log(stats::runif(1)) < la) {
      cf <- cfp
      out$accept[["coefficients"]] <- out$accept[["coefficients"]] + 1
    }
    
    ## --- Gibbs for the coefficient variances, normal prior only ------------
    if (!horseshoe)
      phi2 <- .bayes_rinvgamma(p, pr$a_phi + 0.5, pr$b_phi + cf^2 / 2)
    
    ## --- Metropolis-Hastings for log sigma^2 -------------------------------
    tau <- log(s2)
    dv  <- .bayes_tau_derivs(cf, s2, y, X, q, aee, pr)
    Pc  <- if (dv$f2 < -1e-12) -1 / dv$f2 else 1
    mc  <- tau + Pc * dv$score
    tp  <- stats::rnorm(1, mc, sqrt(Pc))
    s2p <- exp(tp)
    
    dvp <- .bayes_tau_derivs(cf, s2p, y, X, q, aee, pr)
    Pp  <- if (dvp$f2 < -1e-12) -1 / dvp$f2 else 1
    mp  <- tp + Pp * dvp$score
    
    la <- .bayes_lpost(cf, s2p, y, X, q, aee, pr, horseshoe, phi2) -
      .bayes_lpost(cf, s2, y, X, q, aee, pr, horseshoe, phi2) +
      stats::dnorm(tau, mp, sqrt(Pp), log = TRUE) -
      stats::dnorm(tp, mc, sqrt(Pc), log = TRUE) +
      (tp - tau)                                  # Jacobian of tau = log sigma^2
    if (is.finite(la) && log(stats::runif(1)) < la) {
      s2 <- s2p
      out$accept[["sigma2"]] <- out$accept[["sigma2"]] + 1
    }
    
    ## --- Gibbs for W, hence Sigma ------------------------------------------
    ## xi_e is e / sigma directly.  Routing it through qnorm(pnorm(e, sd = s))
    ## is the same number in exact arithmetic but returns +-Inf once |e/sigma|
    ## exceeds about 8.3, where pnorm saturates at 0 or 1.
    Xi[, d] <- (y - as.vector(X %*% cf)) / sqrt(s2)
    dw  <- .bayes_draw_W(Xi[, seq_len(K), drop = FALSE],
                         if (L > 0L) Xi[, K + seq_len(L), drop = FALSE] else NULL,
                         Xi[, d], Om, pr)
    Sig <- dw$Sigma
    Om  <- dw$Omega
    
    ## --- Gibbs for the masses ----------------------------------------------
    Rs <- chol(Sig)
    U  <- stats::pnorm(matrix(stats::rnorm(N * d), N, d) %*% Rs)
    for (j in seq_len(d - 1L)) lam[[j]] <- .bayes_draw_lambda(U[, j], mgs[[j]])
    
    ## --- store -------------------------------------------------------------
    if (it == nxt) {
      keep <- keep + 1L
      out$coefficients[keep, ] <- cf
      out$sigma2[keep]         <- s2
      out$Sigma[keep, ]        <- Sig[lower.tri(Sig)]
      for (j in seq_len(d - 1L)) out$lambda[[j]][keep, ] <- lam[[j]]
      nxt <- nxt + thin
    }
    
    if (verbose && (it %% step == 0L || it == iterations))
      utils::setTxtProgressBar(pb, it)
  }
  
  out$accept <- out$accept / iterations
  out
}


## =============================================================================
##  6.  CopRegBAYES
## =============================================================================
##
##  formula        y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##                 as everywhere else in the toolbox: factors, interactions,
##                 transformations and "-1" all work, and the position around
##                 the "|" decides what is treated as endogenous.  Every
##                 non-intercept column of the design matrix enters the copula,
##                 so a four level factor contributes three columns with two
##                 unique values each, exactly as the weekly dummies of the
##                 paper's application do.
##
##  iterations,    102,000 iterations, the first 2,000 discarded and every
##  burnin, thin   100th of the rest retained, which is the setting of the
##                 paper's simulations and leaves 1,000 draws.  Draws are
##                 thinned while sampling, so only the retained ones are ever
##                 held in memory.
##
##  horseshoe      TRUE uses the horseshoe prior of Equation (5) in the
##                 marginalised form of Equation (6), p(gamma) proportional to
##                 1/|gamma|, on every coefficient including the intercept.
##                 FALSE puts an independent N(0, phi^2) on each coefficient
##                 with phi^2 ~ IG, drawn by a conjugate Gibbs step.
##
##  prior.args     list overriding the hyperparameters, all defaulting to the
##                 paper: a and b for sigma^2 (0.001 each, Section 3.2.2),
##                 Psi1, nu1, Psi2, nu2, a_e, b_e, M0 and V0 for the blocks of
##                 W, and a_phi, b_phi for the normal prior.
##
##  start          list with any of coefficients, sigma2, Sigma, lambda.
##                 Defaults follow Web Appendix D: OLS for the coefficients and
##                 the error variance, the identity for Sigma, and one draw
##                 from Dir(1) per regressor for the masses.
##
##  Use validity() for the identification checks and the convergence
##  diagnostics, including Gelman-Rubin.

#' @title Bayesian Gaussian copula endogeneity correction
#'
#' @description
#' Joint Bayesian estimation of a linear regression together with a Gaussian
#' copula correction for endogenous regressors (Haschka 2025). Unlike the
#' other estimators in the toolbox, the marginal CDFs of the regressors are
#' not plugged in from a first stage: each regressor gets probability masses
#' on its uniquely observed values (Section 3.1 of the paper), and those
#' masses are drawn by MCMC together with the regression coefficients, the
#' residual variance and the full copula correlation matrix. There is
#' therefore no first-stage plug-in uncertainty to bootstrap back in; the
#' posterior carries all of it.
#'
#' There is no \code{cdf} argument and no \code{ties} argument: the CDF of
#' each regressor is a parameter of the model rather than a choice the user
#' makes, and ties are the normal case rather than a problem to work around
#' -- the probability masses sit on unique values, so a binary regressor
#' simply has two of them. There is likewise no \code{nboots} argument:
#' inference here is posterior rather than bootstrap, governed instead by
#' \code{iterations}, \code{burnin} and \code{thin}. Supplying any of
#' \code{cdf}, \code{ties} or \code{nboots} is an error.
#'
#' Exogeneity of the regressors listed after \code{|} is imposed by holding
#' the corresponding entries of the copula covariance matrix at exactly zero.
#' This is done by a block-recursive reparameterisation of the covariance
#' matrix (endogenous regressors regressed on the exogenous ones and on the
#' structural error, all independent by construction) rather than by an
#' inverse Wishart draw, so the zeros hold exactly rather than only
#' approximately.
#'
#' @param formula Two-part formula, \code{y ~ endogenous | exogenous}, as
#'   everywhere else in the toolbox: factors, interactions, transformations
#'   and \code{-1} all work, and the position of a term around \code{|}
#'   decides whether it is treated as endogenous (and so gets a copula term)
#'   or exogenous. Every non-intercept column of the design matrix enters the
#'   copula, so a four-level factor on the exogenous side contributes three
#'   columns with two unique values each.
#' @param data A \code{data.frame} holding the variables in \code{formula}.
#' @param iterations Total number of MCMC iterations. Defaults to 102,000,
#'   the setting used in the paper's simulations.
#' @param burnin Number of iterations discarded from the start of the chain.
#'   Defaults to 2,000.
#' @param thin Keep every \code{thin}-th iteration after burn-in. Defaults to
#'   100, which together with the defaults above leaves 1,000 draws. Draws
#'   are thinned as they are produced rather than after the fact, so only the
#'   retained iterates are ever held in memory.
#' @param horseshoe If \code{TRUE} (the default), a horseshoe prior is placed
#'   on every regression coefficient, including the intercept. If
#'   \code{FALSE}, each coefficient instead gets an independent
#'   \eqn{N(0, \phi^2)} prior with \eqn{\phi^2} drawn by a conjugate Gibbs
#'   step from an inverse Gamma full conditional.
#' @param prior.args A list overriding any of the hyperparameters, all of
#'   which default to the paper's own settings: \code{a}, \code{b} for the
#'   inverse Gamma prior on the residual variance; \code{Psi1}, \code{nu1}
#'   for the exogenous-regressor covariance block; \code{Psi2}, \code{nu2}
#'   for the copula error-covariance block; \code{a_e}, \code{b_e} for the
#'   structural-error variance in that reparameterisation; \code{M0},
#'   \code{V0} for the regression-on-exogenous-and-error prior; and
#'   \code{a_phi}, \code{b_phi} for the normal-prior variance used when
#'   \code{horseshoe = FALSE}.
#' @param start A list overriding any of the starting values (any of
#'   \code{coefficients}, \code{sigma2}, \code{Sigma}, \code{lambda}). The
#'   defaults follow Web Appendix D of the paper: OLS for the coefficients
#'   and the residual variance, the identity matrix for \code{Sigma}, and one
#'   draw from a Dir(1) distribution per regressor for the probability
#'   masses.
#' @param subset An optional logical or index vector selecting the rows of
#'   \code{data} to use, as in \code{lm()}.
#' @param contrasts An optional list passed to \code{model.matrix()} for the
#'   coding of factors.
#' @param verbose If \code{TRUE}, print a starting message and show a text
#'   progress bar while sampling. Defaults to \code{interactive()}.
#' @param cdf Not a valid argument for this estimator: the marginal CDF is a
#'   parameter drawn by the sampler rather than a plug-in choice. Supplying
#'   it raises an error.
#' @param ties Not a valid argument for this estimator, for the same reason
#'   as \code{cdf}. Supplying it raises an error.
#' @param nboots Not a valid argument for this estimator: inference is
#'   posterior rather than bootstrap. Use \code{iterations}, \code{burnin}
#'   and \code{thin} instead. Supplying it raises an error.
#'
#' @return An object of class \code{"copregbayes"}, a list including
#'   \describe{
#'     \item{\code{coefficients}}{Posterior mean of the regression
#'       coefficients (what \code{coef()} returns).}
#'     \item{\code{posterior.median}}{Posterior median of the coefficients.}
#'     \item{\code{std.error}}{Posterior standard deviation of the
#'       coefficients.}
#'     \item{\code{sigma2}}{Posterior mean of the residual variance.}
#'     \item{\code{vcov}}{Posterior covariance matrix of the coefficients
#'       (what \code{vcov()} returns).}
#'     \item{\code{draws}}{Matrix of thinned post-burn-in draws (one row per
#'       retained iteration) of the coefficients, \code{sigma2}, and the free
#'       (non-zero-restricted) entries of the copula correlation matrix; the
#'       basis for \code{confint()}, \code{plot()} and \code{summary()}.}
#'     \item{\code{coefficient.draws}, \code{sigma2.draws}}{The coefficient
#'       and residual-variance draws separately.}
#'     \item{\code{Sigma.draws}}{Matrix of the thinned draws of every
#'       lower-triangular entry of the copula correlation matrix, free and
#'       zero-restricted alike, named \code{"rho(a*, b*)"}.}
#'     \item{\code{Sigma.free}}{Logical vector marking which columns of
#'       \code{Sigma.draws} are free; the rest are held at exactly zero by
#'       the exogeneity restriction.}
#'     \item{\code{lambda.draws}}{A list, one matrix per endogenous or
#'       exogenous regressor entering the copula, of the thinned draws of its
#'       Dirichlet probability masses (one column per unique observed
#'       value).}
#'     \item{\code{margins}}{The unique values, grouping and counts underlying
#'       each regressor's masses (see \code{\link{plot.copregbayes}},
#'       \code{type = "cdf"}).}
#'     \item{\code{acceptance}}{Metropolis-Hastings acceptance rates for the
#'       regression-coefficient and log-residual-variance updates.}
#'     \item{\code{fitted.values}, \code{residuals}}{Fitted values and
#'       residuals at the posterior mean coefficients.}
#'     \item{\code{endo.names}, \code{exog.names}}{Names of the regressors
#'       entering the copula on the endogenous and exogenous side.}
#'     \item{\code{rho.names}}{Names, in \code{draws} and \code{Sigma.draws},
#'       of the correlations between each endogenous regressor's normal score
#'       and that of the structural error -- the endogeneity correlations.}
#'     \item{\code{horseshoe}, \code{prior.args}, \code{iterations},
#'       \code{burnin}, \code{thin}, \code{ndraws}}{The settings the model was
#'       fit with, and the number of retained draws.}
#'     \item{\code{method}, \code{call}}{A description string and the matched
#'       call.}
#'   }
#'   Also carries \code{y}, \code{X}, \code{n}, \code{formula}, \code{terms},
#'   \code{mf}, \code{xlevels}, \code{contrasts} and \code{na.action} for use
#'   by \code{predict()} and related methods.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n   <- 150
#' x   <- rnorm(n)
#' rho <- 0.6
#' eps <- matrix(rnorm(2 * n), n, 2)
#' eps[, 2] <- rho * eps[, 1] + sqrt(1 - rho^2) * eps[, 2]
#' z   <- 1 + x + eps[, 1]              # endogenous, correlated with e below
#' e   <- eps[, 2]
#' y   <- 1 + 2 * z + 0.5 * x + e
#' dat <- data.frame(y = y, z = z, x = x)
#'
#' fit <- CopRegBAYES(y ~ z | x, data = dat,
#'                     iterations = 3000, burnin = 500, thin = 10)
#' print(fit)
#' summary(fit)
#' }
#'
#' @export
CopRegBAYES <- function(formula, data,
                        iterations = 102000,
                        burnin = 2000,
                        thin = 100,
                        horseshoe = TRUE,
                        prior.args = list(),
                        start = NULL,
                        subset = NULL,
                        contrasts = NULL,
                        verbose = interactive(),
                        cdf, ties, nboots) {
  
  cl <- match.call()
  
  if (!missing(cdf) || !missing(ties))
    stop("The Bayesian approach estimates the marginal CDFs rather than ",
         "plugging them in: the\n  probability masses defining them are ",
         "parameters, drawn along with everything else.\n  The 'cdf' and ",
         "'ties' arguments therefore do not apply here.", call. = FALSE)
  if (!missing(nboots))
    stop("Inference here is posterior rather than bootstrap. Use 'iterations',",
         " 'burnin' and\n  'thin' instead of 'nboots'.", call. = FALSE)
  
  iterations <- as.integer(iterations)
  burnin     <- as.integer(burnin)
  thin       <- as.integer(thin)
  if (burnin >= iterations)
    stop("'burnin' must be smaller than 'iterations'.", call. = FALSE)
  if (thin < 1L) stop("'thin' must be at least 1.", call. = FALSE)
  
  info <- endogCopula::.copreg_model(formula, data, subset, contrasts)
  X <- info$X; y <- info$y
  N <- length(y); p <- ncol(X)
  
  ## --- which columns enter the copula, in the order z, x, e -----------------
  endo <- info$endo_cols
  exog <- setdiff(info$exo_cols, if (info$has_intercept) 1L else integer(0))
  K <- length(endo); L <- length(exog)
  if (K < 1L) stop("No endogenous regressor found.", call. = FALSE)
  
  const <- vapply(c(endo, exog), function(j) length(unique(X[, j])) < 2L,
                  logical(1))
  if (any(const))
    stop("Regressor(s) with a single unique value cannot enter the copula: ",
         paste(colnames(X)[c(endo, exog)][const], collapse = ", "),
         ".\n  A constant carries no distributional information.", call. = FALSE)
  
  ## P for the endogenous regressors, W for the exogenous ones and xi for the
  ## structural error, as everywhere else in the toolbox; the paper writes
  ## z, x and e for the same three things.  The stars mark normal scores:
  ## these are correlations of Phi^{-1}(F(.)), which is what the frequentist
  ## estimators report as rho(P*, xi*) as well.
  cop <- list(endo = endo, exog = exog,
              names = c(colnames(X)[endo], colnames(X)[exog], "xi"))
  mgs <- lapply(c(endo, exog), function(j) .bayes_margin(X[, j]))
  names(mgs) <- colnames(X)[c(endo, exog)]
  d <- K + L + 1L
  
  ## --- hyperparameters, defaults as in the paper ---------------------------
  pr <- list(a = 0.001, b = 0.001,                       # sigma^2, Equation (7)
             a_phi = 0.001, b_phi = 0.001,               # normal prior variance
             nu1 = L + 2, Psi1 = diag(max(L, 1L)),       # Sigma_xx
             nu2 = K + 2, Psi2 = diag(K),                # Omega
             a_e = 0.001, b_e = 0.001,                   # sigma_e^2
             M0 = matrix(0, K, L + 1L), V0 = diag(L + 1L))
  ## An explicit NULL means "leave the default", not "delete it": assigning
  ## NULL into a list drops the element, which would leave the sampler without
  ## a hyperparameter it needs.
  prior.args <- prior.args[!vapply(prior.args, is.null, logical(1))]
  if (length(prior.args) > 0L) {
    bad <- setdiff(names(prior.args), names(pr))
    if (length(bad) > 0L)
      stop("Unknown entries in 'prior.args': ", paste(bad, collapse = ", "),
           ".\n  Available: ", paste(names(pr), collapse = ", "), ".",
           call. = FALSE)
    pr[names(prior.args)] <- prior.args
  }
  pr$V0inv <- chol2inv(chol(pr$V0))
  
  ## --- starting values, Web Appendix D --------------------------------------
  ols <- stats::lm.fit(X, y)
  st <- list(coefficients = unname(ols$coefficients),
             sigma2 = sum(ols$residuals^2) / max(1L, N - p),
             Sigma = diag(d),
             lambda = lapply(mgs, function(m)
             { g <- stats::rgamma(m$m, 1); g / sum(g) }),
             phi2 = rep(1000, p))
  start <- start[!vapply(start, is.null, logical(1))]
  if (length(start) > 0L) {
    bad <- setdiff(names(start), names(st))
    if (length(bad) > 0L)
      stop("Unknown entries in 'start': ", paste(bad, collapse = ", "), ".",
           call. = FALSE)
    st[names(start)] <- start
  }
  if (anyNA(st$coefficients))
    stop("The design matrix is rank deficient, so OLS gives no starting ",
         "values.\n  Drop the collinear column(s) or supply 'start'.",
         call. = FALSE)
  
  if (verbose)
    message("Sampling ", format(iterations, big.mark = ","), " iterations, ",
            "keeping ", length(seq.int(burnin + 1L, iterations, by = thin)),
            " draws ...")
  
  fit <- .bayes_sampler(y, X, cop, mgs, iterations, burnin, thin,
                        horseshoe, pr, st, verbose)
  
  ## --- labels for the free elements of Sigma --------------------------------
  ij <- which(lower.tri(diag(d)), arr.ind = TRUE)
  colnames(fit$Sigma) <- paste0("rho(", cop$names[ij[, 2L]], "*, ",
                                cop$names[ij[, 1L]], "*)")
  free <- !(ij[, 2L] > K & ij[, 2L] < d & ij[, 1L] == d)
  
  draws <- cbind(fit$coefficients, sigma2 = fit$sigma2,
                 fit$Sigma[, free, drop = FALSE])
  
  cf <- colMeans(fit$coefficients)
  names(cf) <- colnames(X)
  fitted <- as.vector(X %*% cf)
  
  structure(list(
    coefficients = cf,
    posterior.median = apply(fit$coefficients, 2, stats::median),
    std.error = apply(fit$coefficients, 2, stats::sd),
    sigma2 = mean(fit$sigma2),
    vcov = stats::var(fit$coefficients),
    draws = draws,
    coefficient.draws = fit$coefficients,
    sigma2.draws = fit$sigma2,
    Sigma.draws = fit$Sigma,
    Sigma.free = free,
    lambda.draws = fit$lambda,
    margins = mgs,
    copula.names = cop$names,
    acceptance = fit$accept,
    fitted.values = fitted,
    residuals = y - fitted,
    y = y, X = X, n = N,
    formula = formula, terms = info$terms, mf = info$mf,
    xlevels = info$xlevels, contrasts = info$contrasts,
    na.action = info$na.action,
    endo.names = colnames(X)[endo], exog.names = colnames(X)[exog],
    rho.names = paste0("rho(", colnames(X)[endo], "*, xi*)"),
    horseshoe = horseshoe, prior.args = pr,
    iterations = iterations, burnin = burnin, thin = thin,
    ndraws = nrow(fit$coefficients),
    method = "Bayesian copula correction (Haschka 2025)",
    call = cl), class = "copregbayes")
}


## =============================================================================
##  7.  Methods
## =============================================================================

#' @title Print a copregbayes fit
#'
#' @description
#' Prints the call and the posterior mean of the regression coefficients for
#' a \code{\link{CopRegBAYES}} fit.
#'
#' @param x An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param digits Number of significant digits to print.
#' @param ... Not used; present for S3 method consistency.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' print(fit)
#'
#' @export
print.copregbayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nPosterior means:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n", x$ndraws, " draws from ", format(x$iterations, big.mark = ","),
      " iterations.\n", sep = "")
  invisible(x)
}


#' @title Summarise a copregbayes fit
#'
#' @description
#' \code{summary()} tabulates the posterior mean, median, standard deviation
#' and credible interval of the regression coefficients, the residual
#' variance and the endogeneity correlations; \code{print()} on the result
#' formats that table.
#'
#' @param object An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param level Credible level for the reported interval, e.g. \code{0.95}
#'   for a 95% interval. Defaults to \code{0.95}.
#' @param x An object of class \code{"summary.copregbayes"}, as returned by
#'   \code{summary.copregbayes}.
#' @param digits Number of significant digits to print.
#' @param ... Not used; present for S3 method consistency.
#'
#' @return For \code{summary.copregbayes}, an object of class
#'   \code{"summary.copregbayes"}, a list with \code{coefficients} (the
#'   posterior mean/median/sd/credible-interval table for the regression
#'   coefficients), \code{sigma2} (the same for the residual variance),
#'   \code{rho} (the same for the endogeneity correlations named in
#'   \code{object$rho.names}; correlations among the regressors themselves
#'   are nuisance parameters and are not shown here, though they remain in
#'   \code{object$Sigma.draws}), \code{n.other} (how many of those nuisance
#'   correlations there are), \code{restricted} (names of the correlations
#'   held at exactly zero by the exogeneity restriction), \code{acceptance},
#'   \code{horseshoe}, \code{n}, \code{ndraws}, \code{level},
#'   \code{iterations}, \code{burnin}, \code{thin}, \code{call} and
#'   \code{method}.
#'
#'   For \code{print.summary.copregbayes}, \code{x} is returned invisibly;
#'   called for its side effect of printing.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' summary(fit)
#'
#' @rdname summary.copregbayes
#' @export
summary.copregbayes <- function(object, level = 0.95, ...) {
  a <- (1 - level) / 2
  tab <- t(apply(object$draws, 2, function(v)
    c(`P. Mean` = mean(v), `P. Median` = stats::median(v), `Sd` = stats::sd(v),
      stats::quantile(v, c(a, 1 - a)))))
  ## Only the correlations between an endogenous regressor and the structural
  ## error are shown.  Those are the endogeneity, and the whole point of the
  ## model; the correlations among the regressors themselves are nuisance
  ## parameters that the joint estimation happens to produce as well, and with
  ## twelve regressors there are sixty-six of them.  They stay on the object.
  p   <- length(object$coefficients)
  rho <- tab[rownames(tab) %in% object$rho.names, , drop = FALSE]
  structure(list(call = object$call, method = object$method,
                 coefficients = tab[seq_len(p), , drop = FALSE],
                 sigma2 = tab[p + 1L, , drop = FALSE],
                 rho = rho,
                 n.other = nrow(tab) - p - 1L - nrow(rho),
                 restricted = colnames(object$Sigma.draws)[!object$Sigma.free],
                 acceptance = object$acceptance, horseshoe = object$horseshoe,
                 n = object$n, ndraws = object$ndraws, level = level,
                 iterations = object$iterations, burnin = object$burnin,
                 thin = object$thin),
            class = "summary.copregbayes")
}


#' @rdname summary.copregbayes
#' @export
print.summary.copregbayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nRegression coefficients:\n")
  print(format(x$coefficients, digits = digits), quote = FALSE)
  cat("\nStructural error variance:\n")
  print(format(x$sigma2, digits = digits), quote = FALSE)
  cat("\nEndogeneity: rho(P*, xi*) is the correlation between the normal score",
      "\n  of an endogenous regressor and that of the structural error.\n")
  print(format(x$rho, digits = digits), quote = FALSE)
  if (x$n.other > 0L)
    cat("  ", x$n.other, " further correlations among the regressors are ",
        "estimated jointly and\n  sit in $Sigma.draws; ", length(x$restricted),
        " more are held at zero, which is what makes\n  the exogenous ",
        "regressors exogenous.\n", sep = "")
  cat("\nThe two quantile columns are ", format(100 * x$level), "% credible ",
      "intervals: no asymptotic\n  argument is involved, and the marginal ",
      "CDFs are estimated jointly rather\n  than plugged in.\n", sep = "")
  cat("\n", x$n, " observations. ", x$ndraws, " draws kept from ",
      format(x$iterations, big.mark = ","), " iterations (burn-in ",
      format(x$burnin, big.mark = ","), ", thinning ", x$thin, ").\n", sep = "")
  cat("Prior on the coefficients: ",
      if (x$horseshoe) "horseshoe" else "normal with inverse Gamma variance",
      ". Acceptance ", paste0(names(x$acceptance), " ",
                              formatC(x$acceptance, format = "f", digits = 3),
                              collapse = ", "), ".\n", sep = "")
  invisible(x)
}


#' @title Extract components of a copregbayes fit
#'
#' @description
#' \code{coef()} returns the posterior mean of the regression coefficients,
#' \code{vcov()} their posterior covariance matrix, \code{nobs()} the number
#' of observations used, \code{formula()} the model formula, \code{fitted()}
#' the fitted values at the posterior mean coefficients, and
#' \code{residuals()} the corresponding residuals.
#'
#' @param object An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param x An object of class \code{"copregbayes"} (the argument name used
#'   by the \code{formula()} method).
#' @param ... Not used; present for S3 method consistency.
#'
#' @return For \code{coef.copregbayes}, a named numeric vector: the posterior
#'   mean of each regression coefficient.
#'
#'   For \code{vcov.copregbayes}, a numeric matrix: the posterior covariance
#'   of the regression coefficients.
#'
#'   For \code{nobs.copregbayes}, an integer: the number of observations.
#'
#'   For \code{formula.copregbayes}, the two-part model formula.
#'
#'   For \code{fitted.copregbayes}, a numeric vector of fitted values.
#'
#'   For \code{residuals.copregbayes}, a numeric vector of residuals
#'   (response minus fitted values, at the posterior mean coefficients).
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' coef(fit); vcov(fit); nobs(fit); formula(fit)
#' head(fitted(fit)); head(residuals(fit))
#'
#' @rdname copregbayes-extract
#' @export
coef.copregbayes    <- function(object, ...) object$coefficients
#' @rdname copregbayes-extract
#' @export
vcov.copregbayes    <- function(object, ...) object$vcov
#' @rdname copregbayes-extract
#' @export
nobs.copregbayes    <- function(object, ...) object$n
#' @rdname copregbayes-extract
#' @export
formula.copregbayes <- function(x, ...) x$formula
#' @rdname copregbayes-extract
#' @export
fitted.copregbayes  <- function(object, ...) object$fitted.values

#' @rdname copregbayes-extract
#' @export
residuals.copregbayes <- function(object, ...) object$residuals

#' @title Credible intervals for a copregbayes fit
#'
#' @description
#' Posterior credible intervals, taken directly as quantiles of the retained
#' draws in \code{object$draws} -- not a normal approximation to them, and no
#' asymptotic argument is involved anywhere.
#'
#' @param object An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param parm Optional character or integer vector selecting which
#'   parameters (columns of \code{object$draws}) to return intervals for. If
#'   missing, all parameters are returned.
#' @param level Credible level, e.g. \code{0.95} for a 95% interval.
#'   Defaults to \code{0.95}.
#' @param ... Not used; present for S3 method consistency.
#'
#' @return A numeric matrix with one row per parameter and two columns, the
#'   lower and upper quantiles of the posterior draws at the requested
#'   credible level.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' confint(fit, level = 0.9)
#'
#' @export
## Credible intervals, the posterior quantiles themselves rather than a normal
## approximation to them.
confint.copregbayes <- function(object, parm, level = 0.95, ...) {
  a <- (1 - level) / 2
  ci <- t(apply(object$draws, 2, stats::quantile, c(a, 1 - a)))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @title Predict from a copregbayes fit
#'
#' @description
#' Predictions from the structural regression only: the copula terms never
#' enter \code{predict()} or \code{fitted()} because they are endogeneity
#' controls, not part of the causal model being predicted.
#'
#' @param object An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param newdata An optional \code{data.frame} of new predictor values. If
#'   \code{NULL} (the default), the fitted values of the original data are
#'   returned.
#' @param ... Not used; present for S3 method consistency.
#'
#' @return A numeric vector of predictions, one per row of \code{newdata} (or
#'   of the original data if \code{newdata} is \code{NULL}), computed at the
#'   posterior mean coefficients.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' dat <- data.frame(y = y, z = z, x = x)
#' fit <- CopRegBAYES(y ~ z | x, data = dat,
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' predict(fit, newdata = dat[1:5, ])
#'
#' @export
predict.copregbayes <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  mt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(mt, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X  <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts)
  as.vector(X[, names(object$coefficients), drop = FALSE] %*% object$coefficients)
}


## -----------------------------------------------------------------------------
##  plot
## -----------------------------------------------------------------------------
##  'which' selects what to look at, by name or by position, because a model of
##  any size has far more parameters than fit on a screen; the default is the
##  endogenous regressors, which are the ones the correction is about.
##
##  type = "cdf" is different in kind from the other four.  It shows the
##  marginal distribution of a regressor as estimated by the model -- the
##  cumulative sums of the probability masses, with a pointwise credible band
##  across the draws.  That object simply does not exist in the frequentist
##  estimators, where the CDF is a fixed plug-in with no uncertainty attached.

#' @title Diagnostic plots for a copregbayes fit
#'
#' @description
#' Plots of the MCMC chain (or, for \code{type = "cdf"}, of an estimated
#' marginal distribution) for one or more parameters of a
#' \code{\link{CopRegBAYES}} fit.
#'
#' @param x An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param which Which parameters to plot, by name or by position. For
#'   \code{type} other than \code{"cdf"}, names are matched against
#'   \code{colnames(x$draws)} (regression coefficients, the residual
#'   variance \code{sigma2}, and the free entries of the copula correlation
#'   matrix); for \code{type = "cdf"}, names are matched against
#'   \code{names(x$lambda.draws)} (the regressors entering the copula).
#'   Defaults to \code{x$endo.names}, the endogenous regressors, since those
#'   are what the correction is about.
#' @param type One of \code{"trace"} (the sampled values against iteration
#'   number), \code{"acf"} or \code{"pacf"} (autocorrelation and partial
#'   autocorrelation of the draws), \code{"density"} (a kernel density
#'   estimate of the posterior with the mean marked), or \code{"cdf"}
#'   (the posterior of the regressor's marginal CDF, i.e. the cumulative
#'   sums of its probability masses, shown as a step function with a
#'   pointwise credible band across the draws -- an object that has no
#'   counterpart in the frequentist estimators, where the CDF is a fixed
#'   plug-in with no uncertainty attached). Defaults to \code{"trace"}.
#' @param level Credible level for the band shown when \code{type = "cdf"}.
#'   Defaults to \code{0.95}.
#' @param ask Whether to prompt between plots when more than one is drawn.
#'   Defaults to \code{TRUE} on an interactive, multi-plot device and
#'   \code{FALSE} otherwise.
#' @param ... Further arguments passed to the underlying \code{graphics::plot}
#'   (or \code{stats::acf} / \code{stats::pacf}) call.
#'
#' @return \code{x}, invisibly. Called for its side effect of drawing a plot.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' plot(fit, which = "z", type = "trace")
#' plot(fit, which = "z", type = "cdf")
#'
#' @export
plot.copregbayes <- function(x, which = NULL,
                             type = c("trace", "acf", "pacf", "density", "cdf"),
                             level = 0.95, ask = NULL, ...) {
  
  type <- match.arg(type)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  
  if (type == "cdf") {
    nm <- names(x$lambda.draws)
    if (is.null(which)) which <- x$endo.names
    if (is.numeric(which)) which <- nm[which]
    bad <- setdiff(which, nm)
    if (length(bad) > 0L)
      stop("No margin is estimated for: ", paste(bad, collapse = ", "),
           ".\n  Available: ", paste(nm, collapse = ", "), ".", call. = FALSE)
    a <- (1 - level) / 2
    if (is.null(ask)) ask <- length(which) > 1L && grDevices::dev.interactive()
    graphics::par(ask = ask)
    for (v in which) {
      L <- x$lambda.draws[[v]]
      F <- t(apply(L, 1, cumsum))
      q <- apply(F, 2, stats::quantile, c(a, 0.5, 1 - a))
      u <- x$margins[[v]]$uv
      graphics::plot(u, q[2, ], type = "s", ylim = c(0, 1), xlab = v,
                     ylab = "posterior CDF",
                     main = paste0("Estimated marginal CDF of ", v), ...)
      graphics::lines(u, q[1, ], type = "s", lty = 2, col = "grey40")
      graphics::lines(u, q[3, ], type = "s", lty = 2, col = "grey40")
      graphics::rug(x$margins[[v]]$uv, col = "grey60")
      graphics::legend("bottomright", bty = "n", lty = c(1, 2),
                       col = c("black", "grey40"),
                       legend = c("posterior median",
                                  paste0(format(100 * level), "% band")))
    }
    return(invisible(x))
  }
  
  nm <- colnames(x$draws)
  if (is.null(which)) which <- x$endo.names
  if (is.numeric(which)) which <- nm[which]
  bad <- setdiff(which, nm)
  if (length(bad) > 0L)
    stop("Not a parameter of this model: ", paste(bad, collapse = ", "),
         ".\n  Available: ", paste(nm, collapse = ", "), ".", call. = FALSE)
  
  if (is.null(ask)) ask <- length(which) > 1L && grDevices::dev.interactive()
  graphics::par(ask = ask)
  for (v in which) {
    d <- x$draws[, v]
    switch(type,
           trace   = graphics::plot(seq_along(d), d, type = "l",
                                    xlab = "draw", ylab = v,
                                    main = paste("Trace of", v), ...),
           acf     = stats::acf(d, main = paste("ACF of", v), ...),
           pacf    = stats::pacf(d, main = paste("PACF of", v), ...),
           density = { dd <- stats::density(d)
           graphics::plot(dd, xlab = v,
                          main = paste("Posterior of", v), ...)
           graphics::abline(v = mean(d), lty = 2) })
  }
  invisible(x)
}


## -----------------------------------------------------------------------------
##  validity
## -----------------------------------------------------------------------------
##  Two things at once, both of which the paper treats as prerequisites rather
##  than results.
##
##  Identification, Section 4.2.  The three known failure modes are normality
##  of the endogenous regressor, non-Gaussian regressor-error dependence, and a
##  misspecified structural error distribution.  Only the first is testable
##  from the data alone, and it is the one that matters most: if z is normal
##  then E(e|z) is linear and the copula cannot separate it from the regression
##  itself.  The paper's own reassurance is that the Bayesian version degrades
##  gracefully -- the posterior of rho(z,e) concentrates near zero rather than
##  producing confident nonsense -- so that posterior is reported alongside.
##
##  Convergence.  Geweke, effective sample size and the lag-one autocorrelation
##  come out of the chain that is already there.  Gelman-Rubin does not: it
##  needs several chains from dispersed starts, which means running the sampler
##  again.  Web Appendix E draws those starts from N(OLS, OLS se) for the
##  coefficients, LKJ for the correlation matrix and Dir(1) for the masses, and
##  that is what chains = TRUE does.  It costs as many runs as chains.

.bayes_geweke <- function(v, first = 0.1, last = 0.5) {
  n <- length(v)
  a <- v[seq_len(floor(first * n))]
  b <- v[(n - floor(last * n) + 1L):n]
  s <- function(u) {
    m <- max(1L, floor(10 * log10(length(u))))
    r <- stats::acf(u, lag.max = m, plot = FALSE)$acf[-1L]
    stats::var(u) * (1 + 2 * sum(r * (1 - seq_along(r) / (m + 1)))) / length(u)
  }
  se <- s(a) + s(b)
  if (!is.finite(se) || se <= 0) return(NA_real_)
  (mean(a) - mean(b)) / sqrt(se)
}

.bayes_ess <- function(v) {
  n <- length(v)
  m <- max(1L, floor(10 * log10(n)))
  r <- stats::acf(v, lag.max = m, plot = FALSE)$acf[-1L]
  ## Initial positive sequence: sum the autocorrelations up to the first
  ## negative one and stop there.  k == 1L, a chain whose lag-one
  ## autocorrelation is already negative, truncates to no lags at all and so
  ## to ESS = n; keeping the rest of the sequence instead would add up a tail
  ## of noise and could report several times as many effective draws as there
  ## are draws.
  k <- which(r < 0)[1L]
  if (!is.na(k)) r <- r[seq_len(k - 1L)]
  max(1, n / (1 + 2 * sum(r)))
}

## LKJ(1), the uniform distribution over correlation matrices, by the onion
## method; used only for dispersed starting values.
.bayes_rlkj <- function(d) {
  if (d < 2L) return(diag(d))
  b <- (d - 1) / 2
  r <- 2 * stats::rbeta(1, b, b) - 1
  R <- matrix(c(1, r, r, 1), 2L, 2L)
  if (d > 2L) for (m in 2:(d - 1L)) {
    b <- b - 0.5
    y <- stats::rbeta(1, m / 2, b)
    u <- stats::rnorm(m); u <- u / sqrt(sum(u * u))
    z <- t(chol(R)) %*% (sqrt(y) * u)
    R <- rbind(cbind(R, z), c(as.vector(z), 1))
  }
  R
}

#' @title Identification and convergence checks for a copregbayes fit
#'
#' @description
#' \code{validity()} for a \code{\link{CopRegBAYES}} fit reports two things
#' at once. First, the identification requirement of the model: the copula
#' correction is identified off the nonnormality of the endogenous
#' regressors, since under normality \eqn{E(e|z)} is linear and the copula
#' cannot separate regressor variation from error variation; the posterior of
#' the endogeneity correlations (the correlation between an endogenous
#' regressor's normal score and that of the structural error) is reported
#' alongside, since the Bayesian version is expected to degrade gracefully --
#' concentrating near zero rather than producing confident nonsense -- when
#' identification is weak. Second, convergence of the MCMC chain itself:
#' Geweke's statistic, the effective sample size and the lag-one
#' autocorrelation, all computed from the chain that is already there, plus
#' -- on request, since it needs several chains from dispersed starting
#' values -- the Gelman-Rubin statistic, which costs one further run of the
#' sampler per additional chain.
#'
#' \code{print()} on the result formats all of that.
#'
#' @param object An object of class \code{"copregbayes"}, as returned by
#'   \code{\link{CopRegBAYES}}.
#' @param chains Whether to also compute the Gelman-Rubin statistic from
#'   dispersed starting values. \code{FALSE} (the default) skips it.
#'   \code{TRUE} runs 4 additional chains; a number runs that many additional
#'   chains instead (at least 2). Each additional chain re-runs the sampler
#'   with the same call but a fresh, dispersed starting point drawn from
#'   \eqn{N}(coefficients, posterior sd) for the regression coefficients, an
#'   LKJ(1) draw for the copula correlation matrix, and a fresh Dir(1) draw
#'   per regressor for the probability masses.
#' @param power Statistical power used for the sample-size-dependent
#'   nonnormality thresholds (Becker, Proksch & Ringle 2022). Defaults to
#'   \code{0.8}.
#' @param verbose If \code{TRUE}, print a message before each additional
#'   chain required by \code{chains}. Defaults to \code{interactive()}.
#' @param x An object of class \code{"copregbayes.validity"}, as returned by
#'   \code{validity.copregbayes}.
#' @param digits Number of significant digits to print.
#' @param ... Not used; present for S3 method consistency.
#'
#' @return For \code{validity.copregbayes}, an object of class
#'   \code{"copregbayes.validity"}, a list with \code{nonnormality} (a table
#'   of skewness, excess kurtosis, Anderson-Darling and Cramer-von Mises
#'   statistics and a Kolmogorov-Smirnov p-value for each endogenous
#'   regressor, with pass/fail flags against two thresholds),
#'   \code{thresholds} (the sample-size-and-power-dependent thresholds
#'   used), \code{endogeneity} (posterior mean, 95% interval and
#'   \eqn{P(\rho > 0)} for each endogeneity correlation), \code{convergence}
#'   (a data frame of Geweke, effective sample size and lag-one
#'   autocorrelation for every retained parameter), \code{gelman.rubin} (a
#'   named vector of Gelman-Rubin statistics, or \code{NULL} if
#'   \code{chains = FALSE}), \code{acceptance}, \code{ndraws} and
#'   \code{endo.names}.
#'
#'   For \code{print.copregbayes.validity}, \code{x} is returned invisibly;
#'   called for its side effect of printing.
#'
#' @references
#' Haschka, R. E. (2025). Bayesian inference for joint estimation models
#' using copulas to handle endogenous regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}. \doi{10.1111/obes.70023}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n); z <- x + rnorm(n); y <- 1 + z + x + rnorm(n)
#' fit <- CopRegBAYES(y ~ z | x, data = data.frame(y, z, x),
#'                     iterations = 200, burnin = 50, thin = 5, verbose = FALSE)
#' validity(fit)
#'
#' @importFrom endogCopula validity
#' @rdname validity.copregbayes
#' @export
validity.copregbayes <- function(object, chains = FALSE, power = 0.8,
                                 verbose = interactive(), ...) {
  
  nz <- endogCopula::.validity_nonnormality(object$X[, object$endo.names, drop = FALSE],
                               object$endo.names, object$n, power)
  
  ## posterior of the endogeneity correlations, one per endogenous regressor
  ze <- object$rho.names[object$rho.names %in% colnames(object$draws)]
  rho <- t(apply(object$draws[, ze, drop = FALSE], 2, function(v)
    c(`P. Mean` = mean(v), `2.5%` = unname(stats::quantile(v, 0.025)),
      `97.5%` = unname(stats::quantile(v, 0.975)),
      `P(rho > 0)` = mean(v > 0))))
  
  conv <- data.frame(
    Geweke = apply(object$draws, 2, .bayes_geweke),
    ESS    = apply(object$draws, 2, .bayes_ess),
    `AC(1)` = apply(object$draws, 2, function(v)
      stats::acf(v, lag.max = 1, plot = FALSE)$acf[2L]),
    check.names = FALSE)
  
  gr <- NULL
  if (!isFALSE(chains)) {
    nc <- if (isTRUE(chains)) 4L else as.integer(chains)
    if (nc < 2L) stop("Gelman-Rubin needs at least two chains.", call. = FALSE)
    ## The extra chains re-evaluate the original call, so they have to be
    ## evaluated where that call makes sense: the frame validity() was called
    ## from.  Taking it here rather than inside the lapply() below matters,
    ## since parent.frame() there is lapply()'s own frame, from which the
    ## caller's data is reachable only if it happens to be global.
    where <- parent.frame()
    cl <- object$call
    cl$verbose <- verbose
    se <- object$std.error
    if (verbose)
      message("Gelman-Rubin needs one further run of the sampler per chain; ",
              nc, " to go.")
    reps <- lapply(seq_len(nc), function(i) {
      if (verbose) message("Chain ", i + 1L, " of ", nc + 1L, ":")
      cl$start <- list(
        coefficients = stats::rnorm(length(object$coefficients),
                                    object$coefficients, se),
        Sigma = .bayes_rlkj(length(object$copula.names)),
        lambda = lapply(object$margins, function(m)
        { g <- stats::rgamma(m$m, 1); g / sum(g) }))
      eval(cl, where)$draws
    })
    reps <- c(list(object$draws), reps)
    n <- nrow(reps[[1L]]); M <- length(reps)
    gr <- vapply(seq_len(ncol(reps[[1L]])), function(j) {
      x <- vapply(reps, function(r) r[, j], numeric(n))
      B <- n * stats::var(colMeans(x))
      W <- mean(apply(x, 2, stats::var))
      if (W <= 0) return(NA_real_)
      sqrt(((n - 1) / n * W + B / n) / W)
    }, numeric(1))
    names(gr) <- colnames(object$draws)
  }
  
  structure(list(nonnormality = nz$table, thresholds = nz$thresholds,
                 endogeneity = rho, convergence = conv, gelman.rubin = gr,
                 acceptance = object$acceptance, ndraws = object$ndraws,
                 endo.names = object$endo.names),
            class = "copregbayes.validity")
}


#' @rdname validity.copregbayes
#' @export
print.copregbayes.validity <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       ...) {
  i <- 0L; nxt <- function() paste0("[", i <<- i + 1L, "] ")
  cat("\nChecks for the Bayesian copula correction\n")
  
  cat("\n", nxt(), "Nonnormality of the endogenous regressors\n", sep = "")
  cat("    Identification comes from the nonlinearity of E(e|z), which is\n",
      "    present only when z is nonnormal; under normality the copula cannot\n",
      "    tell regressor variation from error variation.\n", sep = "")
  print(format(x$nonnormality, digits = digits))
  
  cat("\n", nxt(), "Posterior of the endogeneity correlations\n", sep = "")
  print(format(x$endogeneity, digits = digits), quote = FALSE)
  cat("    An interval covering zero says the data carry no evidence of\n",
      "    endogeneity, which is the honest reading; it is not a test decision.\n",
      sep = "")
  
  cat("\n", nxt(), "Convergence of the chain\n", sep = "")
  cat("    Geweke compares the first tenth with the last half of the draws and\n",
      "    is a z statistic, so |Geweke| > 2 is a warning sign. ESS is the\n",
      "    effective number of independent draws behind ", x$ndraws,
      " retained ones.\n", sep = "")
  print(format(x$convergence, digits = digits))
  cat("    Acceptance rates: ",
      paste0(names(x$acceptance), " ",
             formatC(x$acceptance, format = "f", digits = 3), collapse = ", "),
      "\n", sep = "")
  
  if (!is.null(x$gelman.rubin)) {
    cat("\n", nxt(), "Gelman-Rubin from dispersed starts\n", sep = "")
    print(format(x$gelman.rubin, digits = digits), quote = FALSE)
    cat("    Gelman et al. read convergence as GR <= 1.1.\n")
  } else {
    cat("\n", nxt(), "Gelman-Rubin was not computed\n", sep = "")
    cat("    It needs several chains from dispersed starts, so it costs one\n",
        "    further run of the sampler per chain. Pass chains = TRUE, or a\n",
        "    number, to compute it.\n", sep = "")
  }
  
  cat("\nSource: Haschka (2025), Oxford Bulletin of Economics and Statistics\n")
  invisible(x)
}


### REFERENCES ----------------------------------------------------------------
##
## Carvalho, C. M., N. G. Polson, and J. G. Scott (2010). The horseshoe
##   estimator for sparse signals. Biometrika 97(2), 465-480.
##
## Gamerman, D. (1997). Sampling from the posterior distribution in generalized
##   linear mixed models. Statistics and Computing 7(1), 57-68.
##
## Haschka, R. E. (2025). Bayesian inference for joint estimation models using
##   copulas to handle endogenous regressors. Oxford Bulletin of Economics and
##   Statistics.
##
## Ng, K. W., G. Tian, and M. Tang (2011). Dirichlet and Related Distributions.
##   Chichester: Wiley.
##
## Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##   estimation using copulas. Marketing Science 31(4), 567-586.
##
## Zhang, X., W. J. Boscardin, and T. R. Belin (2006). Sampling correlation
##   matrices in Bayesian models with correlated latent variables. Journal of
##   Computational and Graphical Statistics 15(4), 880-896.