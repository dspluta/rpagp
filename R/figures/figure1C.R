library(reshape2)
library(ggplot2)
library(dplyr)
library(readr)
library(magrittr)
library(png)
source("R/plot_functions.R")
source("R/LPP_analysis_functions.R")
options(dplyr.summarise.inform = FALSE)

###################
# Set data variables
n_iter <- 5000
n <- 166
n_time <- 99
subj_id <- 30763

#############
dat_raw <- read_delim("data/Only_WEIGHT_PicByPic_ERPsnocigPresentatioOrder.txt",
                      delim = "\t",
                      col_names = T)

selected_raw_categories <- c("PH", "UH")
dat <- dat_raw %>% filter(SubjectID == subj_id,
                          CATOKnumFD %in% selected_raw_categories) %>%
  tidyr::pivot_longer(cols = 8:232) %>%
  mutate(time = as.numeric(name),
         category = recode_category(CATOKnumFD)) %>%
  dplyr::select(-name, -EXPNUM, -NewNamePicOK, -PLENUEUNP, -CATOKnumFD, -SubjectID)

dat$value[dat$time < 50] <- dat$value[dat$time < 50] * seq(0, 1, length.out = sum(unique(dat$time) < 50))

y <- stats::filter(dat |> group_by(time) |> summarize(value = mean(value)) %$% value, filter = rep(1/15, 15))
dat_plot <- data.frame(y = y, time = unique(dat$time))

ggplot(dat_plot) +
  geom_line(aes(x = time, y = y), size = 2) +
  ylim(-2.5, 8.5) +
  xlim(-100, 700) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave(ggsave("figures/figure1/figure1C_f.png", height = 5, width = 7, dpi = 640))

ggplot(dat_plot) +
  geom_line(aes(x = time, y = 1.3 * y), size = 2) +
  geom_line(aes(x = time, y = y), size = 1, alpha = 0.5) +
  ylim(-2.5, 9.5) +
  xlim(-100, 700) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave(ggsave("figures/figure1/figure1C_beta_f.png", height = 5, width = 7, dpi = 640))

ggplot(dat_plot) +
  geom_line(aes(x = time - 32, y = 1.3 * y), size = 2) +
  geom_line(aes(x = time, y = 1.3 * y), size = 1, alpha = 0.5) +
  geom_line(aes(x = time, y = y), size = 1, alpha = 0.5) +
  ylim(-2.5, 9.5) +
  xlim(-100, 700) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave(ggsave("figures/figure1/figure1C_beta_f_tau.png", height = 5, width = 7, dpi = 640))

ggplot(dat_plot) +
  geom_line(aes(x = time - 32, y = 1.3 * y + 0.2 * arima.sim(n = nrow(dat_plot), model = list(ar = c(0.5, 0.1)))), size = 2) +
  geom_line(aes(x = time - 32, y = 1.3 * y), size = 1, alpha = 0.5) +
  geom_line(aes(x = time, y = 1.3 * y), size = 1, alpha = 0.5) +
  geom_line(aes(x = time, y = y), size = 1, alpha = 0.5) +
  ylim(-2.5, 9.5) +
  xlim(-100, 700) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave(ggsave("figures/figure1/figure1C_y.png", height = 5, width = 7, dpi = 640))



