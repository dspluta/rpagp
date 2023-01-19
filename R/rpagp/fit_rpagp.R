
#' Fit RPAGP Model
#'
#' @param y data
#' @param n_iter number of posterior samples
#' @param theta0 parameter initializations
#' @param hyperparam hyperparameters
#' @param pinned_point pinned point for identifiability of structural signal
#' @param pinned_value value of f at the pinned point
#'
fit_rpagp <- function(y, n_iter, theta0, hyperparam, pinned_point, pinned_value = 1) {
  chain <- vector(mode = "list", length = n_iter)
  chain_f <- vector(mode = "list", length = n_iter)
  chain_mu <- vector(mode = "list", length = n_iter)
  chain_z <- vector(mode = "list", length = n_iter)
  x <- seq(0, 1, length.out = nrow(y))
  z <- matrix(0, nrow = n_time, ncol = n)

  chain[[1]] <- theta0
  chain_f[[1]] <- sample_f(y, chain[[1]], n_draws = 1, nugget = 1e-6)
  start <- Sys.time()
  for (iter in 2:n_iter) {
    if (iter %% 100 == 0) cat(iter / 100)

    # Sample f and rescale
    f <- sample_f(y, chain[[iter - 1]], n_draws = 1)
    f <- pinned_value * f / f[pinned_point]
    K_f_inv <- solve(sq_exp_kernel(x, chain[[iter - 1]]$rho, nugget = 1e-9))

    # Update y hat
    y_hat <- get_y_hat_matrix(y, f, chain[[iter - 1]], K_f_inv)

    # Sample betas
    chain[[iter]]$beta <- sample_betas(y, f,
                                       chain[[iter - 1]], K_f_inv,
                                       y_hat, hyperparam)

    # Update y hat
    for (i in 1:n) y_hat[, i] <- chain[[iter]]$beta[i] * y_hat[, i] / chain[[iter - 1]]$beta[i]

    # Sample taus and rho
    chain[[iter]]$tau <- sample_tau(y, f, chain[[iter - 1]], hyperparam, K_f_inv)
    chain[[iter]]$rho <- sample_rho(y, f, chain[[iter - 1]], hyperparam, K_f_inv)

    # Update y hat
    y_hat <- get_y_hat_matrix(y, f, chain[[iter]], K_f_inv)

    # Compute residuals and sample residual parameters
    z <- y - y_hat
    ar_post <- get_ar_posterior(z, p = 2, n_draws = 1)

    # Record draws from current iteration
    chain_f[[iter]] <- f
    chain[[iter]]$phi <- ar_post[[1]][[1]]
    chain[[iter]]$sigma <- sqrt(ar_post[[1]][[2]])
    chain_mu[[iter]] <- y_hat
    chain_z[[iter]] <- z
  }
  cat("\n")
  end <- Sys.time()
  print(end - start)
  return(list(chain = chain,
              chain_f = chain_f,
              chain_mu = chain_mu,
              chain_z = chain_z))
}
