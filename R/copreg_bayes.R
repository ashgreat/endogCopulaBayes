#' Bayesian Gaussian Copula Regression
#'
#' Internal helpers supporting the exported `CopRegBayes` sampler.
#' @noRd

# auxiliary functions
aux1 <- function(h, param, W1, data) {
  
  N <- nrow(data)
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  W01 <- W1
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  
  M01 <- matrix(ncol = 1, nrow = 3,
                data = c(qnorm(data$z1)[h], qnorm(data$x1)[h], 
                         qnorm(pnorm(data$y - a - b1*data$z - b2*data$x, sd = sqrt(param[4])))[h]))
  
  P01 <- t(M01)%*%(solve(W01) - diag(3))%*%M01
  return(P01)
  
}

# (log) posterior with normal priors on regression coefficients (and IG hyperpriors),
# IG for all variances (residual variance, variance of prior for regression
# coefficients), inverse Wishart for Copula correlation matrix, categorial 
# distribution for regressors, Dirichlet priors for probability masses 
# formalising PDF of regressors' distribution
post2 <- function(param, W, data) {
  
  N <- nrow(data)
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  s <- param[4]
  p1 <- param[5]
  p2 <- param[6]
  p3 <- param[7]
  
  sa <- param[(7 + 2*N) + 1]
  sb1 <- param[(7 + 2*N) + 2]
  sb2 <- param[(7 + 2*N) + 3]
  
  e <- data$y - a - b1*data$z - b2*data$x
  
  X_hat <- matrix(0, nrow = N, ncol = 3)
  X_hat[, 1] <- data$z1
  X_hat[, 2] <- data$x1
  X_hat[, 3] <- pnorm(q = e, mean = 0, sd = sqrt(s))
  
  Phi1 <- matrix(nrow = 3, ncol = 3, byrow = TRUE,
                 data = c(1, p1, p2,
                          p1, 1, p3,
                          p2, p3, 1))
  
  
  t1 <- sum(dCopula(copula = normalCopula(param = P2p(Phi1), dim = 3, 
                                           dispstr = "un"), u = X_hat, 
                     log = TRUE)) + 
           sum(dnorm(mean = 0, sd = sqrt(s), x = e, log = TRUE)) +
    sum(LaplacesDemon::dcat(x = c(table(data$z)), p = data$z1g, log = TRUE)) +
    sum(LaplacesDemon::dcat(x = c(table(data$x)), p = data$x1g, log = TRUE)) +
    dinvwishart(Sigma = W, nu = 3, S = diag(3), log = TRUE) + 
    invgamma::dinvgamma(x = s, shape = .001, rate = .001, log = TRUE) + 
    LaplacesDemon::ddirichlet(x = t(data$z1g), alpha = table(data$z) + rep(1, N), log = TRUE) +
    LaplacesDemon::ddirichlet(x = t(data$x1g), alpha = table(data$x) + rep(1, N), log = TRUE) +
    dnorm(x = a, mean = 0, sd = sqrt(sa), log = TRUE) + 
    dnorm(x = b1, mean = 0, sd = sqrt(sb1), log = TRUE) +
    dnorm(x = b2, mean = 0, sd = sqrt(sb2), log = TRUE) +
    invgamma::dinvgamma(x = sa, shape = .001, rate = .001, log = TRUE) +
    invgamma::dinvgamma(x = sb1, shape = .001, rate = .001, log = TRUE) +
    invgamma::dinvgamma(x = sb2, shape = .001, rate = .001, log = TRUE)
  return(t1)
  
}

# proposals
proposals1 <- function(param, W1, data) {
  
  N <- nrow(data)
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]

  P <- param[4]*solve(t(cbind(data$const, qnorm(data$z1), qnorm(data$x1))) %*% cbind(data$const, qnorm(data$z1), qnorm(data$x1)))
  
  # Proposals for beta
  p1 <- rmvnorm(n = 1, mean = c(a, b1, b2), sigma = P)
  # p1 <- rmvt(n = 1, delta = c(a, b1, b2), sigma = P, df = 1)

  return(p1)
  
}
proposals2 <- function(param, W1, data) {
  
  N <- nrow(data)
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  
  # Proposals for sigma2
  P01 <- sapply(X = 1:N, FUN = aux1, param = param, W1 = W1, data = data)
  
  p2 <- invgamma::rinvgamma(n = 1, shape = N/2, 
                            rate = abs(sum((data$y - a - b1*data$z - b2*data$x)^2)/2 - sum(P01)/2))
  
  return(p2)
  
}
d_proposals <- function(param, param1, W1, data) {
  
  N <- nrow(data)
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]

  a0 <- param1[1]
  b01 <- param1[2]
  b02 <- param1[3]

  P01 <- sapply(X = 1:N, FUN = aux1, param = param1, W1 = W1, data = data)
  
  p2 <- invgamma::dinvgamma(x = param[4], shape = N/2, log = TRUE,
                            rate = abs(sum((data$y - a - b1*data$z - b2*data$x)^2)/2 - sum(P01)/2))
  
  return(p2)
  
}

# MCMC algorithm
metropolis_Gibbs_MCMC1 <- function(startvalue, iterations, data) {
  
  dat1 <- data
  N <- nrow(dat1)

  chain <- matrix(data = NA, nrow = iterations + 1, ncol = (7 + 2*N) + 3)
  
  chain[1, ] <- startvalue
  
  W01 <- diag(3)
  W1 <- diag(3)
  
  
  for (i in 1:iterations) {
    
    
    ################################# Proposals ################################
    
    proposal <- chain[i, ]
    proposal[1:3] <- proposals1(param = chain[i, ], W1 = W01, data = dat1)
    
    
    # acceptance probability
    probab <- min(exp(post2(param = proposal, W = W1, data = dat1) - 
                        post2(param = chain[i, ], W = W1, data = dat1)), 1)
    if (is.nan(probab)) {probab <- 0}
    
    # Metropolis step betas
    if (runif(1) < probab) {
      
      chain[i + 1, ] <- proposal
      
      dat1$z1g <- proposal[8:(N+7)]
      dat1 <- dat1[order(dat1$z), ]
      dat1$z1 <- cumsum(dat1$z1g)*(length(dat1$z1g) / (length(dat1$z1g) + .01))
      
      dat1$x1g <- proposal[(N+8):(7 + 2*N)]
      dat1 <- dat1[order(dat1$x), ]
      dat1$x1 <- cumsum(dat1$x1g)*(length(dat1$x1g) / (length(dat1$x1g) + .01))
      
    } else {
      
      chain[i + 1, ] <- chain[i, ]
      
      dat1$z1g <- proposal[8:(N+7)]
      dat1 <- dat1[order(dat1$z), ]
      dat1$z1 <- cumsum(dat1$z1g)*(length(dat1$z1g) / (length(dat1$z1g) + .01))
      
      dat1$x1g <- proposal[(N+8):(7 + 2*N)]
      dat1 <- dat1[order(dat1$x), ]
      dat1$x1 <- cumsum(dat1$x1g)*(length(dat1$x1g) / (length(dat1$x1g) + .01))
      
    }
    
    
    # Metropolis step sigma2
    proposal <- chain[i + 1, ]
    proposal[4] <- proposals2(param = chain[i + 1, ], W1 = W01, data = dat1)
    
    probab <- min(exp(post2(param = proposal, W = W1, data = dat1) - 
                        post2(param = chain[i + 1, ], W = W1, data = dat1) + 
                        d_proposals(param = chain[i + 1, ], param1 = proposal, W1 = W01, data = dat1) - 
                        d_proposals(param = proposal, param1 = chain[i + 1, ], W1 = W01, data = dat1)), 1)
    if (is.nan(probab)) {probab <- 0}
    
    if (runif(1) < probab) {
      
      chain[i + 1, ] <- proposal
      
    } 
    
    
    ############################ Gibbs step Wishart ############################
    
    X01 <- matrix(ncol = 3, nrow = N, 
                  data = c(qnorm(dat1$z1),
                           qnorm(dat1$x1),
                           qnorm(pnorm(dat1$y - chain[i+1, 1] - 
                                         chain[i+1, 2]*dat1$z - chain[i+1, 3]*dat1$x, 
                                       sd = sqrt(chain[i+1, 4])))))
    
    W1 <- rinvwishart(nu = N + 3, S = diag(3) + t(X01)%*%X01)
    
    while (1 > 0) {
      
      W1 <- rinvwishart(nu = N + 3, S = diag(3) + t(X01)%*%X01) 
      W1[2, 3] <- 0
      W1[3, 2] <- 0
      
      if (min(eigen(W1)$values) > 0) { break }
      
    }
    
    W01 <- solve(sqrt(diag(W1))*diag(3))%*%W1%*%solve(sqrt(diag(W1))*diag(3))
    
    chain[i+1, 5:7] <- P2p(W01)
    # chain[i+1, 5] <- cor(qnorm(pobs(X[, 2])), qnorm(pobs(X[, 3])))
    # chain[i+1, 7] <- 0
    # W1[2, 3] <- 0
    # W1[3, 2] <- 0
    
    
    ########################### Gibbs step Dirichlet ###########################
    
    eps_n <- mvtnorm::rmvnorm(n = N, mean = rep(0, 3), sigma = W01, method = "eigen")
    
    pt_z <- qgamma(p = pnorm(eps_n[, 1]), shape = table(dat1$z) + rep(1, length(dat1$z)), rate = 1)
    chain[i+1, c(8:(N + 7))] <- pt_z/sum(pt_z)
    
    pt_x <- qgamma(p = pnorm(eps_n[, 2]), shape = table(dat1$x) + rep(1, length(dat1$x)), rate = 1)
    chain[i+1, c((N + 8):(2*N + 7))] <- pt_x/sum(pt_x)
    
    
    dat1$z1g <- NA
    dat1$z1 <- NA
    dat1$x1g <- NA
    dat1$x1 <- NA
    
    ### Hyperpriors
    
    chain[i + 1, (7 + 2*N) + 1] <- invgamma::rinvgamma(n = 1, shape = .001 + .5,
                                                       rate = .5*(chain[i + 1, 1]^2 + 2*.001))
    chain[i + 1, (7 + 2*N) + 2] <- invgamma::rinvgamma(n = 1, shape = .001 + .5,
                                                       rate = .5*(chain[i + 1, 2]^2 + 2*.001))
    chain[i + 1, (7 + 2*N) + 3] <- invgamma::rinvgamma(n = 1, shape = .001 + .5,
                                                       rate = .5*(chain[i + 1, 3]^2 + 2*.001))
    
  }
  
  return(chain)
  
}

#' Bayesian Gaussian Copula Sampler
#'
#' Metropolis-within-Gibbs sampler for the copula-based endogeneity
#' correction of Haschka. The chain columns are laid out as in the reference
#' implementation: columns 1--3 hold the regression coefficients (`beta_0`,
#' `beta_z`, `beta_x`), column 4 the residual variance (`sigma2`), columns
#' 5--7 the copula correlations (`rho_zx` between endogenous and exogenous
#' regressor, `rho_ze` between endogenous regressor and error, `rho_xe`
#' between exogenous regressor and error, the latter forced to zero), columns
#' `8:(N + 7)` the Dirichlet probability masses for the distribution of `z`,
#' columns `(N + 8):(2 * N + 7)` the Dirichlet probability masses for the
#' distribution of `x`, and the final three columns the hyperprior variances
#' of the normal priors on the regression coefficients (`hyper_a`,
#' `hyper_b1`, `hyper_b2`). The Dirichlet mass columns are left unnamed.
#'
#' @param data Data frame containing the columns `y`, `z`, and `x`.
#' @param iterations Total number of MCMC iterations.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after burn-in.
#' @param startvalue Optional numeric vector of starting values; if `NULL`, a
#'   default based on OLS estimates and Dirichlet draws is used.
#' @param seed Optional integer passed to [set.seed()] before any random
#'   number is drawn; if `NULL` (default) the current RNG state is used.
#' @return An object of class `endog_copula_bayes`: a list with components
#'   \describe{
#'     \item{`chain`}{Numeric matrix of dimension `(iterations + 1) x
#'       (2 * N + 10)` holding the full MCMC chain, one draw per row. The
#'       columns follow the layout described above: the seven named model
#'       parameters, then the `2 * N` unnamed Dirichlet probability masses
#'       (`N` for `z`, `N` for `x`), then the three hyperprior variances.}
#'     \item{`posterior`}{Numeric matrix with the same columns as `chain`,
#'       after discarding the first `burnin` rows and keeping every
#'       `thin`-th remaining row.}
#'     \item{`parameters`}{Character vector with the names of the ten tracked
#'       (named) parameters summarised by [summary.endog_copula_bayes()].}
#'     \item{`n`}{Number of complete observations used.}
#'     \item{`iterations`, `burnin`, `thin`, `seed`}{The sampler settings as
#'       supplied.}
#'   }
#' @seealso [summary.endog_copula_bayes()] and
#'   [print.endog_copula_bayes()] for posterior summaries.
#' @examples
#' set.seed(1)
#' n <- 60
#' z <- rlnorm(n)
#' x <- rnorm(n)
#' dat <- data.frame(y = 2 - 4 * z + 6 * x + rnorm(n, sd = sqrt(2)),
#'                   z = z, x = x)
#' fit <- CopRegBayes(dat, iterations = 100, burnin = 20, thin = 2, seed = 42)
#' fit
#' summary(fit)
#' @importFrom copula dCopula normalCopula P2p
#' @importFrom mvtnorm rmvnorm
#' @importFrom LaplacesDemon rinvwishart dinvwishart
#' @importFrom stats coef complete.cases dnorm lm median pnorm qgamma qnorm
#'   quantile residuals runif sd var
#' @export
CopRegBayes <- function(data, iterations = 10000, burnin = 2000, thin = 10,
                        startvalue = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dataset <- as.data.frame(data)
  required <- c("y", "z", "x")
  missing <- setdiff(required, names(dataset))
  if (length(missing) > 0) {
    stop(sprintf("Missing variables in data: %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  dataset <- dataset[, required]
  dataset <- dataset[stats::complete.cases(dataset), , drop = FALSE]
  N <- nrow(dataset)
  if (N == 0) {
    stop("No complete observations available.", call. = FALSE)
  }
  dataset$const <- 1
  if (is.null(startvalue)) {
    mod1 <- stats::lm(y ~ z + x, dataset)
    startvalue <- c(
      stats::coef(mod1),
      stats::var(stats::residuals(mod1)),
      rep(0, 3),
      c(LaplacesDemon::rdirichlet(1, rep(1, length(dataset$z)))),
      c(LaplacesDemon::rdirichlet(1, rep(1, length(dataset$x)))),
      rep(1000, 3)
    )
  }
  if (length(startvalue) != (7 + 2 * N) + 3) {
    stop("Length of 'startvalue' does not match expected dimension.", call. = FALSE)
  }
  chain <- metropolis_Gibbs_MCMC1(startvalue = startvalue, iterations = iterations,
                                  data = dataset)
  main_names <- c("beta_0", "beta_z", "beta_x", "sigma2",
                  "rho_zx", "rho_ze", "rho_xe")
  hyper_names <- c("hyper_a", "hyper_b1", "hyper_b2")
  colnames(chain) <- c(main_names, rep("", 2 * N), hyper_names)
  keep <- chain
  if (burnin >= nrow(chain)) {
    warning("Burn-in exceeds chain length; returning full chain.", call. = FALSE)
  } else {
    keep <- chain[(burnin + 1):nrow(chain), , drop = FALSE]
  }
  if (thin > 1 && nrow(keep) > 0) {
    keep <- keep[seq(1, nrow(keep), by = thin), , drop = FALSE]
  }
  structure(
    list(
      chain = chain,
      posterior = keep,
      parameters = c(main_names, hyper_names),
      n = N,
      iterations = iterations,
      burnin = burnin,
      thin = thin,
      seed = seed
    ),
    class = "endog_copula_bayes"
  )
}

#' Print a Bayesian Copula Regression Fit
#'
#' @param x An object of class `endog_copula_bayes`.
#' @param digits Number of significant digits to print.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @seealso [CopRegBayes()], [summary.endog_copula_bayes()]
#' @export
print.endog_copula_bayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     ...) {
  cat("Bayesian Gaussian copula regression (Haschka)\n")
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Iterations: %d (burn-in %d, thinning %d, %d posterior draws)\n",
              x$iterations, x$burnin, x$thin, nrow(x$posterior)))
  cat("\nPosterior means:\n")
  means <- colMeans(x$posterior[, x$parameters, drop = FALSE])
  print(round(means, digits))
  invisible(x)
}

#' Summarise a Bayesian Copula Regression Fit
#'
#' Posterior summaries (mean, standard deviation, median, and 95% credible
#' interval) of the named model parameters, computed from the thinned
#' posterior sample.
#'
#' @param object An object of class `endog_copula_bayes`.
#' @param ... Ignored.
#' @return An object of class `summary.endog_copula_bayes`: a list with the
#'   matrix of summary `statistics` (one row per named parameter, columns
#'   `mean`, `sd`, `median`, `2.5%`, and `97.5%`), the number of observations
#'   `n`, the number of retained posterior draws `n_draws`, and the
#'   `iterations`, `burnin`, and `thin` settings.
#' @seealso [CopRegBayes()]
#' @export
summary.endog_copula_bayes <- function(object, ...) {
  draws <- object$posterior[, object$parameters, drop = FALSE]
  stats_mat <- t(apply(draws, 2, function(v) {
    c(mean(v), stats::sd(v), stats::median(v),
      stats::quantile(v, probs = c(0.025, 0.975), names = FALSE))
  }))
  dimnames(stats_mat) <- list(object$parameters,
                              c("mean", "sd", "median", "2.5%", "97.5%"))
  structure(
    list(
      statistics = stats_mat,
      n = object$n,
      n_draws = nrow(draws),
      iterations = object$iterations,
      burnin = object$burnin,
      thin = object$thin
    ),
    class = "summary.endog_copula_bayes"
  )
}

#' Print a Summary of a Bayesian Copula Regression Fit
#'
#' @param x An object of class `summary.endog_copula_bayes`.
#' @param digits Number of significant digits to print.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.summary.endog_copula_bayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                             ...) {
  cat("Bayesian Gaussian copula regression (Haschka)\n")
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Iterations: %d (burn-in %d, thinning %d, %d posterior draws)\n",
              x$iterations, x$burnin, x$thin, x$n_draws))
  cat("\nPosterior summary:\n")
  print(round(x$statistics, digits))
  invisible(x)
}
