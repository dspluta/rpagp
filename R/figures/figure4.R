library(reshape2)
library(ggplot2)
library(dplyr)
library(readr)
library(magrittr)
library(png)
source("R/plot_functions.R")
source("R/LPP_analysis_functions.R")

###################
# Set data variables
n_iter <- 5000
n <- 166
n_time <- 99
subj_id <- 30763

# Read raw data
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

# Set up data frame of trial info for selected subject
trial_info <- dat |> filter(SubjectID == subj_id, time == -100) |>
  dplyr::select(-time, value) |>
  arrange(prognum)
trial_info$trial <- 1:nrow(trial_info)

# Load results for single subject analysis
load("data/dat_CI.RData")
load(glue::glue("data/subj{subj_id}/RPAGP_results_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_bs_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_f_post_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_model_predictions_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_LPP_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_post_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_beta_{subj_id}.RData"))

dat_CI$snr_EMP <- abs(dat_CI$est_EMP) / dat_CI$width_EMP
dat_CI$snr_RPAGP <- abs(dat_CI$est_RPAGP) / dat_CI$width_RPAGP
dat_CI$lwr_EMP <- dat_CI$est_EMP - 0.5 * dat_CI$width_EMP
dat_CI$upr_EMP <- dat_CI$est_EMP + 0.5 * dat_CI$width_EMP
dat_CI$lwr_RPAGP <- dat_CI$est_RPAGP - 0.5 * dat_CI$width_RPAGP
dat_CI$upr_RPAGP <- dat_CI$est_RPAGP + 0.5 * dat_CI$width_RPAGP
dat_CI$Subject_Order <- rank(-dat_CI$est_RPAGP)

#### Fig 4 ####
# Figure 4A (right): Point estimates scatterplot
ggplot(dat_CI) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(x = est_EMP, y = est_RPAGP, fill = factor(color)),
             size = 4.1, pch = 21, alpha = 0.5) +
  labs(x = "EMP Difference of LPP Means", y = "RPAGP Difference of LPP Means",
       fill = "Significant Difference") +
  scale_fill_manual(values = lpp_palette, limits = c("RPAGP & EMP", "RPAGP Only", "Neither"),
                     breaks = c("RPAGP & EMP", "RPAGP Only", "Neither")) +
  ylim(-5.5, 8.5) +
  xlim(-5.5, 8.5) +
  my_theme +
  theme(legend.position = c(0.27, 0.89),
        axis.text = element_text(size = 22, color = "black"),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        axis.ticks = element_blank())
ggsave("figures/figure4/figure4A_point_estimates_scatterplot.png", height = 7, width = 7, dpi = 640)

# Figure 4A (right): SNR scatterplot
ggplot(dat_CI) +
  geom_point(aes(x = snr_EMP, y = snr_RPAGP, fill = color), pch = 21, size = 4.1,
             alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "EMP SNR", y = "RPAGP SNR",
       fill = "Significant Difference") +
  scale_fill_manual(values = lpp_palette, limits = c("RPAGP & EMP", "RPAGP Only", "Neither"),
                    breaks = c("RPAGP & EMP", "RPAGP Only", "Neither")) +
  my_theme +
  theme(legend.position = c(0.27, 0.89),
        axis.text = element_text(size = 22, color = "black"),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        axis.ticks = element_blank())
ggsave("figures/figure4/figure4A_snr_scatterplot.png", height = 7, width = 7, dpi = 640)

# Figure 4B (right): RPAGP CI Widths of Group Differences
ggplot(dat_CI) +
  geom_errorbar(aes(x = -Subject_Order, ymin = lwr_RPAGP, ymax = upr_RPAGP, color = color), size = 2) +
  geom_point(aes(x = -Subject_Order, y = est_RPAGP)) +
  scale_color_manual(values = lpp_palette, limits = c("RPAGP & EMP", "RPAGP Only", "Neither"),
                     breaks = c("RPAGP & EMP", "RPAGP Only", "Neither")) +
  coord_flip() +
  labs(color = "", x = "", y = "RPAGP Group Differences") +
  my_theme +
  theme(legend.position = c(0.2, 0.85),
        axis.text = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.2, 'cm'))
ggsave("figures/figure4/figure4B_RPAGP_estimated_group_differences.png", height = 7, width = 7, dpi = 640)

# Figure 4B (left): EMP CI Widths of Group Differences
ggplot(dat_CI) +
  geom_errorbar(aes(x = -Subject_Order, ymin = lwr_EMP, ymax = upr_EMP, color = color), size = 2,
                alpha = 0.5) +
  geom_point(aes(x = -Subject_Order, y = est_EMP)) +
  scale_color_manual(values = lpp_palette, limits = c("RPAGP Only", "RPAGP & EMP", "Neither"),
                     breaks = c("RPAGP Only", "RPAGP & EMP", "Neither")) +
  coord_flip() +
  labs(color = "", x = "Rank (RPAGP LPP Category Difference)", y = "EMP Category Differences") +
  my_theme +
  theme(legend.position = 'none',
        axis.text = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
ggsave("figures/figure4/figure4B_EMP_estimated_group_differences.png", height = 7, width = 7, dpi = 640)
