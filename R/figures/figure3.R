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

# Filter categories of interest, recode to 3 levels
# Reformat to include column for time
# Select needed columns
dat <- tidyr::pivot_longer(dat_raw, cols = 8:232) %>%
  filter(CATOKnumFD %in% c("PH", "pl", "ul", "UH", "np")) |>
  mutate(category = recode_category(CATOKnumFD),
         time = as.numeric(name)) |>
  dplyr::select(SubjectID, prognum, category, time, value, PresentationOrder) |>
  arrange(SubjectID, prognum, time, PresentationOrder)

trial_info <- dat |> filter(SubjectID == subj_id, time == -100) |>
  dplyr::select(-time, value) |>
  arrange(prognum)
trial_info$trial <- 1:nrow(trial_info)

# EEG scalp map for plots
img <- readPNG(source = "figures/EEG_scalp_map.PNG")
g <- rasterGrob(img, interpolate=TRUE)

# Load results for single subject analysis
load("data/dat_CI.RData")
load(glue::glue("data/subj{subj_id}/RPAGP_results_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_bs_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_f_post_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_model_predictions_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_LPP_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_post_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_beta_{subj_id}.RData"))

# Set up CI interval plot data
dat_CI$snr_EMP <- abs(dat_CI$est_EMP) / dat_CI$width_EMP
dat_CI$snr_RPAGP <- abs(dat_CI$est_RPAGP) / dat_CI$width_RPAGP
dat_CI$lwr_EMP <- dat_CI$est_EMP - 0.5 * dat_CI$width_EMP
dat_CI$upr_EMP <- dat_CI$est_EMP + 0.5 * dat_CI$width_EMP
dat_CI$lwr_RPAGP <- dat_CI$est_RPAGP - 0.5 * dat_CI$width_RPAGP
dat_CI$upr_RPAGP <- dat_CI$est_RPAGP + 0.5 * dat_CI$width_RPAGP
dat_CI$Subject_Order <- rank(-dat_CI$est_RPAGP)

#### Fig 3 ####
# Figure 3A (right): Density of Group Differences (RPAGP)
ggplot(dat_LPP) +
  geom_density(aes(x = value, fill = name), alpha = 0.5) +
  scale_fill_manual(values = lpp_palette[c("H_N", "H_L", "L_N")],
                    breaks = c("H_N", "H_L", "L_N"),
                    labels = c("High - Neutral", "Low - Neutral", "High - Low")) +
  labs(fill = "", y = "", x = "Difference of LPP Means") +
  xlim(0, 7) +
  my_theme +
  theme(legend.position = c(0.2, 0.85),
        axis.text = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 22))
ggsave(ggsave("figures/figure3/figure3A_RPAGP_LPP_mean_differences_density_30763.png",
              height = 5, width = 7, dpi = 640))

# Figure 3A (left): Density of LPP Means (RPAGP)
ggplot(dat_post) +
  geom_density(aes(x = value, fill = category), alpha = 0.5) +
  labs(x = "LPP Means of Trial Predictions", fill = '', y = "Density") +
  scale_fill_manual(values = lpp_palette[c(1, 2, 3)], breaks = c("High", "Low", "Neutral")) +
  my_theme +
  theme(legend.position = c(0.75, 0.85),
        axis.text = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 24))
ggsave(ggsave("figures/figure3/figure3A_RPAGP_LPP_mean_densities_30763.png",
              height = 5, width = 7, dpi = 640))

# Figure 3B (right): High vs Neutral LPP Means
ggplot(dat_beta |> mutate(order = rank(PresentationOrder)) |> filter(category != "Low")) +
  geom_smooth(aes(x = order, y = value, color = category), alpha = 0.4, fill = "grey") +
  geom_errorbar(aes(x = order, ymin = lwr, ymax = upr, group = category), size = 1) +
  geom_point(aes(x = order, y = value, fill = category, pch = category), size = 3) +
  scale_shape_discrete(palette = function(x) {return(c(23, 21))},
                       limits = c("High", "Neutral")) +
  scale_x_discrete(limits = c(1, 20, 40, 60, 77)) +
  scale_color_manual(values = lpp_palette[c(1, 3)], breaks = c("High",  "Neutral")) +
  scale_fill_manual(values = lpp_palette[c(1, 3)], breaks = c("High", "Neutral")) +
  labs(x = "Trial", y = "RPAGP LPP Mean", fill = "") +
  ylim(c(-1.1, 2.5)) +
  my_theme +
  theme(legend.position = 'none',
        axis.text = element_text(size = 26, color = "black", hjust = 0.8),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 24))
ggsave(ggsave("figures/figure3/figure3B_RPAGP_amplitude_CIs_longitudinal_H_N_30763.png",
              height = 5, width = 7, dpi = 640))

# Figure 3B (left): Low vs Neutral LPP Means Longitudinal
ggplot(dat_beta |> mutate(order = rank(PresentationOrder)) |> filter(category != "High")) +
  geom_smooth(aes(x = order, y = value, color = category), alpha = 0.4, fill = "grey") +
  geom_errorbar(aes(x = order, ymin = lwr, ymax = upr, group = category), size = 1) +
  geom_point(aes(x = order, y = value, fill = category, pch = category), size = 3) +
  scale_x_discrete(limits = c(2, 20, 40, 60, 76)) +
  scale_shape_discrete(palette = function(x) {return(c(24, 21))}, limits = c("Low", "Neutral")) +
  scale_color_manual(values = lpp_palette[c(1, 2)], breaks = c("Low",  "Neutral")) +
  scale_fill_manual(values = lpp_palette[c(1 ,2)], breaks = c("Low", "Neutral")) +
  labs(x = "Trial", y = "RPAGP LPP Mean", fill = "") +
  ylim(c(-1.1, 2.5)) +
  my_theme +
  theme(axis.text = element_text(size = 26, color = "black", hjust = 0.8),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        legend.position = 'none')
ggsave(ggsave("figures/figure3/figure3B_RPAGP_amplitude_CIs_longitudinal_L_N_30763.png",
              height = 5, width = 7, dpi = 640))
