#' Generate data.
#'
#' @param n Number of trials.
#' @param n_time Number of time points.
#' @param theta Named list of parameter values.
generate_data <- function(n, n_time, theta) {
  x <- seq(0, 1, length.out = n_time)
  K_f <- sq_exp_kernel(x, theta$rho, nugget = 1e-6)
  K_f_inv <- solve(K_f)

  f <- MASS::mvrnorm(n = 1, rep(0, n_time), K_f)
  y <- matrix(NA, nrow = n_time, ncol = n)
  z <- matrix(NA, nrow = n_time, ncol = n)
  mu <- matrix(NA, nrow = n_time, ncol = n)
  for (i in 1:n) {
    K_i <- get_K_i(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
    mu[, i] <- K_i %*% K_f_inv %*% f
    z[, i] <- arima.sim(model = list(ar = theta$phi), sd = theta$sigma, n = n_time)
    y[, i] <- mu[, i] + z[, i]
  }
  return(list(y = y, f = f, z = z, mu = mu))
}

#' Title
#'
#' @param chain List of MCMC draws.
#' @param var_name Parameter name.
#' @param trial_categories Trial categories/conditions.
#' @param scaling Value to scale parameter samples by.
extract_parameter_chain <- function(chain, var_name, trial_categories, scaling = 1) {
  n_iter <- length(chain)
  n <- length(chain[[1]][var_name][[1]])
  chain_matrix <- matrix(nrow = n, ncol = n_iter)
  for (iter in 1:n_iter) {
    chain_matrix[, iter] <- scaling * chain[[iter]][var_name][[1]]
  }
  chain <- reshape2::melt(chain_matrix, varnames = c("trial", "iter")) %>%
    mutate(category = (trial_categories %$% category)[trial],
           trial = (trial_categories %$% trial)[trial])
  return(chain)
}
