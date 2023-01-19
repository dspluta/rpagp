
#' Sample tau.
#'
#' @param y Matrix of observed trial data (n_time x n).
#' @param f Vector of f values (n_time).
#' @param current Named list of current parameter values.
#' @param hyperparam Named list of hyperparameter values.
#' @param K_f_inv Inverse covariance matrix of f (n_time x n_time).
sample_tau <- function(y, f, current, hyperparam, K_f_inv) {
  proposed <- current
  proposed$tau <- propose_tau(current$tau, hyperparam$tau_proposal_sd)

  lik_current <- likelihood(y, f, current, K_f_inv)
  prior_current <- prior$tau(current$tau, hyperparam$tau_prior_sd)

  lik_proposed <- likelihood(y, f, proposed, K_f_inv)
  prior_proposed <- prior$tau(proposed$tau, hyperparam$tau_prior_sd)

  prob <- exp(lik_proposed + prior_proposed - lik_current - prior_current)
  if (prob > runif(1)) {
    return(proposed$tau)
  } else {
    return(current$tau)
  }
}

#' Sample rho.
#'
#' @param y Matrix of observed trial data (n_time x n).
#' @param f Vector of f values (n_time).
#' @param current Named list of current parameter values.
#' @param hyperparam Named list of hyperparameter values.
#' @param K_f_inv Inverse covariance matrix of f (n_time x n_time).
sample_rho <- function(y, f, current, hyperparam, K_f_inv) {
  proposed <- current
  proposed$rho <- propose_rho(current$rho, hyperparam$rho_proposal_sd)

  lik_current <- likelihood(y, f, current, K_f_inv)
  prior_current <- prior$rho(current$rho, hyperparam$rho_prior_shape, hyperparam$rho_prior_scale)

  lik_proposed <- likelihood(y, f, proposed, K_f_inv)
  prior_proposed <- prior$rho(proposed$rho, hyperparam$rho_prior_shape, hyperparam$rho_prior_scale)

  prob <- exp(lik_proposed + prior_proposed - lik_current - prior_current)
  if (prob > runif(1)) {
    return(proposed$rho)
  } else {
    return(current$rho)
  }
}

#' Sample betas.
#'
#' @param y Matrix of observed trial data (n_time x n).
#' @param f Vector of f values (n_time).
#' @param current Named list of current parameter values.
#' @param K_f_inv Inverse covariance matrix of f (n_time x n_time).
#' @param y_hat
#' @param hyperparam Named list of hyperparameter values.
sample_betas <- function(y, f, current, K_f_inv, y_hat, hyperparam) {
  n <- ncol(y)
  betas <- c()
  for (i in 1:n) {
    betas[i] <- sample_blr(y[, i], y_hat[, i] / current$beta[i],
                           mu = hyperparam$beta_prior_mu, sigma = hyperparam$beta_prior_sd,
                           V = diag(1, 1), a = 1, b = 1)$beta
  }
  return(betas)
}


#' One draw from the posterior of a Bayesian Linear Regression.
#'
#' @param y Response vector.
#' @param X Design matrix.
#' @param mu Prior mean of coefficients.
#' @param sigma Prior standard deviation of coefficients.
#' @param V Prior covariance of coefficients.
#' @param a Hyperparameter for noise variance.
#' @param b Hyperparameter for noise variance.
sample_blr <- function(y, X, mu, sigma, V, a, b) {
  n <- length(y)
  V_post <- solve(solve(V) + t(X) %*% X)
  mu_post <- V_post %*% (solve(V) %*% mu + t(X) %*% y)
  a_post <- a + n / 2
  b_post <- b + 1 / 2 * (mu %*% solve(V) %*% mu + t(y) %*% y - t(mu_post) %*% solve(V_post) %*% mu_post)

  beta <- mvtnorm::rmvt(n = 1, sigma = (b_post[1] / a_post) * V_post,
                        df = 2 * a_post, delta = mu_post)
  return(list(beta = beta,
              sigma2 = invgamma::rinvgamma(1, shape = a_post, rate = b_post)))
}

#' Title
#'
#' @param z Matrix of ongoing activity residuals (n_time x n).
#' @param p Autoregressive order.
#' @param n_draws Number of draws.
get_ar_posterior <- function(z, p, n_draws) {
  V <- diag(1, p)
  mu <- rep(0, p)
  sigma <- 1
  a <- 1
  b <- 1

  dat <- format_AR_data(z, p, n_draws)

  chain <- list()
  for (iter in 1:n_draws) {
    chain[[iter]] <- sample_blr(dat$y, dat$X, mu, sigma, V, a, b)
  }

  return(chain)
}

#' Sample from posterior of f.
#'
#' @param y Matrix of observed trial data (n_time x n).
#' @param theta Named list of parameter values.
#' @param n_draws Numbver of draws.
#' @param nugget GP covariance nugget.
sample_f <- function(y, theta, n_draws, nugget = 1e-6) {
  n_time <- nrow(y)
  n <- ncol(y)
  chain_f <- vector(mode = "list", length = n_draws)
  x <- seq(0, 1, length.out = n_time)
  K_f <- sq_exp_kernel(x, theta$rho, nugget = nugget)
  K_f_inv <- solve(K_f)

  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    A <- K_f_inv
    b <- matrix(0, n_time)
    Sigma_i_inv <- diag(1 / theta$sigma^2, n_time)
    for (i in 1:n) {
      K_i <- get_K_i(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
      L <- K_i %*% K_f_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[, i] %*% G)
    }
    K_f_post <- solve(A)
    chain_f[[iter]] <- MASS::mvrnorm(n = 1, K_f_post %*% b, K_f_post)
  }
  f_draws <- matrix(unlist(lapply(chain_f, `[`)), nrow = n_time, ncol = n_draws)
  return(f_draws)
}
