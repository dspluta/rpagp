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
load(glue::glue("data/subj{subj_id}/RPAGP_results_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_bs_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_f_post_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_model_predictions_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_LPP_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_post_{subj_id}.RData"))
load(glue::glue("data/subj{subj_id}/dat_beta_{subj_id}.RData"))

############
# Fig 2

# Figure 2A (left): Bootstrap confidence intervals for subject means
ggplot() +
  geom_ribbon(data = dat_bs %>% group_by(category, time) %>% summarize(lwr = quantile(value, 0.05),
                                                                       upr = quantile(value, 0.95)),
              aes(x = time, ymin = lwr, ymax = upr, group = category), alpha = 0.3) +
  geom_line(data = dat_bs %>% group_by(category, time) %>% summarize(value = mean(value)),
            aes(x = time, y = value, color = category), size = 2) +

  labs(x = "Time (ms)", y = "Average ERP (\U003BCV)", linetype = "", color = "") +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  ylim(c(-5, 9.5)) +
  xlim(c(300, 700)) +
  my_theme + theme(panel.grid = element_blank(),
                   legend.position = 'none',
                   axis.title = element_text(size = 22),
                   legend.text = element_text(size = 18),
                   legend.title = element_text(size = 20),
                   axis.text = element_text(size = 20, color = "black"),
                   panel.border = element_blank(),
                   axis.line.x = element_line(),
                   axis.line.y = element_line())
ggsave(glue::glue("figures/figure2/figure2A_bs_category_means_{subj_id}.png"), height = 5, width = 5, dpi = 640)

# Figure 2A (right): Raw Data Trials Plot
ggplot(dat |> filter(SubjectID == subj_id) |> group_by(category, time, prognum) |> summarize(value = mean(value))) +
  geom_line(aes(x = time, y = value, color = category, group = factor(prognum)), size = 1, alpha = 0.8) +
  xlim(300, 700) +
  ylim(-17, 28) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(color = "", x = "Time (ms)", y = "Trial-Level ERP (\u03BCV)") +
  my_theme + theme(panel.grid = element_blank(),
                     legend.position = 'none',
                     axis.title = element_text(size = 22),
                     legend.text = element_text(size = 18),
                     legend.title = element_text(size = 20),
                     axis.text = element_text(size = 20, color = "black"),
                     panel.border = element_blank(),
                     axis.line.x = element_line(),
                     axis.line.y = element_line())
ggsave(glue::glue("figures/figure2/figure2A_raw_data_30763.png"), height = 5, width = 5, dpi = 640)


# Figure 2B (right): RPAGP Posterior Median ERPs
ggplot(dat_f_post %>% filter(iter > n_iter / 2) %>% group_by(time, category) %>% summarize(med = median(value),
                                                                                           lwr = quantile(value, 0.025),
                                                                                           upr = quantile(value, 0.975))) +
  geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, group = category), alpha = 0.35) +
  geom_line(aes(x = time, y = med, color = category), size = 2) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(x = "Time (ms)", y = "Posterior Median ERP (\u03BCV)", linetype = "", color = "") +
  ylim(c(-5, 9.5)) +
  xlim(c(300, 700)) +
  my_theme + theme(panel.grid = element_blank(),
                   legend.position = 'none',
                   axis.title = element_text(size = 22),
                   legend.text = element_text(size = 18),
                   legend.title = element_text(size = 20),
                   axis.text = element_text(size = 20, color = "black"),
                   panel.border = element_blank(),
                   axis.line.x = element_line(),
                   axis.line.y = element_line())
ggsave(glue::glue("figures/figure2/figure2B_RPAGP_category_posterior_30763.png"), height = 5, width = 5, dpi = 640)

# Figure 2B (left): RPAGP Trial Predictions
ggplot(dat_subj) +
  geom_line(aes(x = time, y = pred, group = prognum, color = category), alpha = 0.8, size = 1.1) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  scale_fill_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(x = "Time (ms)", y = "RPAGP Trial Predictions (\U003BCV)", color = "") +
  ylim(-17, 28) +
  xlim(c(300, 700)) +
  my_theme + theme(panel.grid = element_blank(),
                   legend.position = 'none',
                   axis.title = element_text(size = 22),
                   legend.text = element_text(size = 18),
                   legend.title = element_text(size = 20),
                   axis.text = element_text(size = 20, color = "black"),
                   panel.border = element_blank(),
                   axis.line.x = element_line(),
                   axis.line.y = element_line())
ggsave(glue::glue("figures/figure2/figure2B_RPAGP_trial_predictions_30763.png"), height = 5, width = 5, dpi = 640)
