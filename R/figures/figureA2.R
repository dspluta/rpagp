library(dplyr)
library(ggplot2)
files.sources = list.files("R/rpagp/", full.names = T)
sapply(files.sources, source)

# Calculate SNR
calc_SNR <- function(tau_sigma, sigma, beta_sigma, delta_beta = 0.5) {
  return(delta_beta / (beta_sigma + 2 * tau_sigma))
}

# Load results for single subject analysis
load("data/dat_tests.RData")

# Calculate SNR for each simulation
sim_SNR <- list()
sim_SNR[[1]] <- calc_SNR(tau_sigma = 0,
                         sigma = 0.1,
                         beta_sigma = 0.01)
sim_SNR[[3]] <- calc_SNR(tau_sigma = 0,
                         sigma = 0.1,
                         beta_sigma = 0.1)
sim_SNR[[5]] <- calc_SNR(tau_sigma = 0,
                         sigma = 0.1,
                         beta_sigma = 0.2)
sim_SNR[[6]] <- calc_SNR(tau_sigma = 0.1,
                         sigma = 0.1,
                         beta_sigma = 0.01)
sim_SNR[[8]] <- calc_SNR(tau_sigma = 0.1,
                         sigma = 0.1,
                         beta_sigma = 0.1)
sim_SNR[[10]] <- calc_SNR(tau_sigma = 0.1,
                          sigma = 0.1,
                          beta_sigma = 0.2)
sim_SNR[[11]] <- calc_SNR(tau_sigma = 0.2,
                          sigma = 0.1,
                          beta_sigma = 0.01)
sim_SNR[[13]] <- calc_SNR(tau_sigma = 0.2,
                          sigma = 0.1,
                          beta_sigma = 0.1)
sim_SNR[[15]] <- calc_SNR(tau_sigma = 0.2,
                          sigma = 0.1,
                          beta_sigma = 0.2)

# Set up data for simulation power plot
dat_tests_plt <- dat_tests

dat_tests_plt$SNR <- NA
for (i in 1:nrow(dat_tests_plt)) {
  dat_tests_plt$SNR[i] <- sim_SNR[[dat_tests_plt$sim_id[i]]]
}

dat_tests_plt <- dat_tests_plt %>% pivot_longer(cols = c("RPAGP", "Empirical", "anova"))

# Figure A2: Simulation Power Plot
ggplot(dat_tests_plt %>% group_by(n, sim_id, SNR, name) %>%
         filter(n != 4, sim_id %in% c(1, 3, 5, 6, 8, 10, 11, 13, 15)) %>%
         summarize(power = mean(value))) +
  geom_point(aes(x = 1/SNR, y = power, shape = name), size = 3.1, stroke = 1.1) +
  scale_shape_manual(limits = c("RPAGP", "Empirical"), values = c(1, 4)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "Power", x = "1/SNR", shape = "") +
  facet_grid(~factor(n)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14))
ggsave("figures/supplement/figureA2_sim_power_plot_SNR.png",
       dpi = 300, width = 8, height = 6)
