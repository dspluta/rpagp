#' Likelihood
#'
#' @param y Matrix of observed trial data (n_time x n).
#' @param f Vector of values for f (n_time).
#' @param theta Named list of parameter values.
likelihood <- function(y, f, theta) {
  n <- ncol(y)
  n_time <- nrow(y)
  p <- length(theta$phi)
  result <- 0
  z <- matrix(0, nrow = n_time, ncol = n)
  eps <- matrix(0, nrow = n_time, ncol = n)
  for (i in 1:n) {
    z[, i] <- y[, i] - get_y_hat(i, f, theta)
    for (j in (p + 1):n_time) {
      for (k in 1:p) {
        eps[j, k] <- eps[j, k] - theta$phi[k] * eps[j - k, i]
      }
    }
  }

  result <- sum(dnorm(eps, mean = 0, sd = theta$sigma, log = T))
  return(result)
}


