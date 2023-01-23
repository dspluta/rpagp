library(dplyr)
library(ggplot2)
files.sources = list.files("R/rpagp/", full.names = T)
sapply(files.sources, source)

################
# Set data generation parameters and generate data
n <- 20
n_time <- 50
n_sims <- 50
n_iter <- 5000
seed_start <- 123

theta <- list(rho = 10,
              beta = rep(c(1, 2), n / 2),
              tau = rnorm(n, 0, 0.1),
              phi = c(0.5, 0.1),
              sigma = 0.25)

seed <- 123 + i
set.seed(seed)
dat <- generate_data(n, n_time, theta)
dat_trials <- reshape2::melt(dat$y, varnames = c("time", "trial"))

###########
# Run RPAGP
hyperparam <- list(tau_prior_sd = 0.2, tau_proposal_sd = 0.0001,
                   rho_prior_shape = 10, rho_prior_scale = 1,
                   rho_proposal_sd = 1, beta_prior_mu = 1, beta_prior_sd = 0.2)

theta0 <- list(rho = 10,
               beta = c(rep(1, n)),
               tau = rep(0, n),
               phi = c(0.01, 0),
               sigma = 0.5)

results <- fit_rpagp(y = dat$y, n_iter = n_iter,
                     theta0 = theta0, hyperparam = hyperparam,
                     pinned_point = 25, pinned_value = 2)

#########
# Prepare f posterior data frame, scaled by mean of betas
dat_f_post <- reshape2::melt(matrix(unlist(lapply(results$chain_f, `[`)), nrow = n_time, ncol = n_iter), varnames = c("s", "iter"))
dat_f_post$time <- dat_f_post$s

beta_means <- extract_parameter_chain(results$chain, "beta",
                                      data.frame(trial = 1:n,
                                                 category = rep(c("A", "B"), n / 2))) |>
  group_by(iter) |>
  summarize(value = mean(value))

dat_f_post$value <- dat_f_post$value * beta_means$value[dat_f_post$iter]

#########
# Figure A1
ggplot(dat_trials) +
  geom_line(aes(x = time, y = value, group = factor(trial)), alpha = 0.25) +
  geom_line(data = dat_trials |> group_by(time) |>
              summarize(value = mean(value)),
            aes(x = time, y = value), linetype = 3, size = 1.2) +
  geom_line(data = data.frame(time = 1:50, value = 1.5 * sim_results[[8]]$f),
            aes(x = time, y = value), size = 1.2) +
  geom_line(data = dat_f_post %>% filter(iter > burn_in) %>% group_by(time) %>%
              summarize(value = median(value)),
            aes(x = time, y = value), size = 1.2, linetype = 2) +
  labs(y = "", x = "Time") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(glue::glue("figures/supplement/figureA1.png"), height = 5, width = 8, dpi = 640)
