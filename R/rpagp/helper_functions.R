#' Get Cov(y_i, f)
#'
#' @param x sequence of observation times
#' @param theta list with named entries rho and tau
get_K_i <- function(x, theta) {
  n_time <- length(x)
  K <- matrix(NA, n_time, n_time)
  for (i in 1:n_time) {
    for (j in 1:n_time) {
      K[i, j] <- exp(-theta$rho^2 / 2 * (x[i] - x[j] - theta$tau)^2)
    }
  }
  return(theta$beta * K)
}

#' Get y hat
#'
#' @param i Subject index (scalar, 1 < i < n).
#' @param f Vector of values for f.
#' @param theta Named list of parameter values.
#' @param K_f_inv Inverse covariance matrix of f.
get_y_hat <- function(i, f, theta, K_f_inv) {
  x <- seq(0, 1, length.out = length(f))
  K_i <- get_K_i(x, list(tau = theta$tau[i], beta = theta$beta[i], rho = theta$rho))
  mu <- K_i %*% K_f_inv %*% f
  return(mu)
}

#' Get y_hat for all trials, output in matrix form
#'
#' @param y Matrix of observed trial data.
#' @param f Vector of f values.
#' @param theta Parameter values.
#' @param K_f_inv Inverse covariance matrix of f.
get_y_hat_matrix <- function(y, f, theta, K_f_inv) {
  y_hat <- matrix(nrow = nrow(y), ncol = ncol(y))

  for (i in 1:ncol(y_hat)) {
    y_hat[, i] <- get_y_hat(i, f, theta, K_f_inv)
  }
  return(y_hat)
}


#' Format linear regression design matrix from autoregressive data.
#'
#' @param z Matrix of ongoing activity (n_time x n).
#' @param p Autoregressive order.
format_AR_data <- function(z, p) {
  n <- ncol(z)
  n_time <- nrow(z)

  X <- array(NA, dim = c(n_time - p, p, n))
  for (i in 1:n) {
    for (j in 1:p) {
      X[, j, i] <- z[j:(n_time - p + j - 1), i]
    }
  }

  y <- c(z[(p + 1):n_time, ])
  X <- apply(X, 2, c)
  return(list(y = y, X = X))
}

#' Generate covariance matrix for square exponential kernel.
#'
#' @param x Vector of time points.
#' @param rho Length scale.
#' @param alpha Amplitude.
#' @param nugget Covariance nugget.
sq_exp_kernel <- function(x, rho, alpha = 1, nugget = 0.0) {
  K <- toeplitz(alpha^2 * exp(-rho^2 / 2 * x^2))
  diag(K) <- diag(K) + nugget
  return(K)
}
