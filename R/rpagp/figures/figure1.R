library(reshape2)
library(ggplot2)
library(dplyr)
library(readr)
library(magrittr)
library(png)
library(grid)
source("functions_AR.R")
source("LPP_analysis_functions.R")
source("plot_functions.R")

###################
n_iter <- 5000
n <- 166
n_time <- 99
subj_id <- 30763

# Read raw data
dat_raw <- read_delim("data/Only_WEIGHT_PicByPic_ERPsnocigPresentatioOrder.txt")
colnames(dat_raw)

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
img <- readPNG(source = "fig/EEG_scalp_map.PNG")
g <- rasterGrob(img, interpolate=TRUE)

# Load results for single subject analysis
load(file = glue::glue("data/RPAGP_results_{subj_id}.RData"))
load(glue::glue("data/dat_bs_subj_{subj_id}.RData"))
load(glue::glue("data/subj_{subj_id}_dat_f_post.RData"))
load("data/dat_CI.RData")
load("data/dat_subj_model_predictions_30763.RData")
load("data/dat_LPP_30763.RData")
load("data/dat_post_30763.RData")
load("data/dat_f_post_30763.RData")
load(glue::glue("data/dat_beta_{subj_id}.RData"))

dat_CI$snr_EMP <- abs(dat_CI$est_EMP) / dat_CI$width_EMP
dat_CI$snr_RPAGP <- abs(dat_CI$est_RPAGP) / dat_CI$width_RPAGP
dat_CI$lwr_EMP <- dat_CI$est_EMP - 0.5 * dat_CI$width_EMP
dat_CI$upr_EMP <- dat_CI$est_EMP + 0.5 * dat_CI$width_EMP
dat_CI$lwr_RPAGP <- dat_CI$est_RPAGP - 0.5 * dat_CI$width_RPAGP
dat_CI$upr_RPAGP <- dat_CI$est_RPAGP + 0.5 * dat_CI$width_RPAGP
dat_CI$Subject_Order <- rank(-dat_CI$est_RPAGP)

############################
#### Fig 1 ####

# Subj 30763 High Category Plot
ggplot(dat %>% filter(SubjectID == subj_id, category == "High") |> group_by(category, time) %>% summarize(value = mean(value))) +
  geom_line(data = dat |> filter(SubjectID == subj_id, category == "High"),
            aes(x = time, y = value, group = factor(prognum)), alpha = 0.5) +
  geom_line(aes(x = time, y = value, color = category), linewidth = 2.9) +
  geom_vline(xintercept = c(300, 700), linetype = 2) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(x = "", y = "", color = "Stimulus Category") +
  ylim(-25, 30) +
  my_theme +
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 26, color = "black"),
        axis.title.y = element_text(size = 30))
ggsave(glue::glue("fig/subject_{subj_id}_means_high.png"), height = 5, width = 8, dpi = 640)

# Subj 30763 Low Category Plot
ggplot(dat %>% filter(SubjectID == subj_id, category == "Low") |> group_by(category, time) %>% summarize(value = mean(value))) +
  geom_line(data = dat |> filter(SubjectID == subj_id, category == "Low"),
            aes(x = time, y = value, group = factor(prognum)), alpha = 0.5) +
  geom_line(aes(x = time, y = value, color = category), size = 2.9) +
  geom_vline(xintercept = c(300, 700), linetype = 2) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(x = "", y = "Trial-level ERP (\U003BCV)", color = "Stimulus Category") +
  ylim(-25, 30) +
  my_theme +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 18),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 26, color = "black"),
        axis.title.y = element_text(size = 30))
ggsave(glue::glue("fig/subject_{subj_id}_means_low.png"), height = 5, width = 8, dpi = 640)

# Subj 30763 Neutral Category Plot
ggplot(dat %>% filter(SubjectID == subj_id, category == "Neutral") |> group_by(category, time) %>% summarize(value = mean(value))) +
  geom_line(data = dat |> filter(SubjectID == subj_id, category == "Neutral"),
            aes(x = time, y = value, group = factor(prognum)), alpha = 0.5) +
  geom_line(aes(x = time, y = value, color = category), size = 2.9) +
  geom_vline(xintercept = c(300, 700), linetype = 2) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(x = "Time (ms)", y = "", color = "Stimulus Category") +
  ylim(-25, 30) +
  my_theme +
  theme(legend.position = 'none',
        axis.text = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30))
ggsave(glue::glue("fig/subject_{subj_id}_means_neutral.png"), height = 5, width = 8, dpi = 640)

# Subj 30763 Category Means
ggplot(dat |> filter(SubjectID == subj_id) |> group_by(category, time) |> summarize(value = mean(value))) +
  geom_line(aes(x = time, y = value, color = category), size = 2.2) +
  geom_vline(xintercept = c(300, 700), linetype = 2, size = 1.1) +
  annotate(geom = "text", label = "LPP ROI", x = 550, y = -2.9, size = 7) +
  scale_color_manual(values = lpp_palette, breaks = c("High", "Low", "Neutral")) +
  labs(color = "", x = "Time (ms)", y = "Average ERP (\u03BCV)") +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = 'none',
                     legend.text = element_text(size = 16),
                     legend.title = element_text(size = 18),
                     panel.border = element_blank(),
                     axis.line.x = element_line(),
                     axis.line.y = element_line(),
                     axis.text = element_text(size = 26, color = "black"),
                     axis.title = element_text(size = 30))
ggsave(glue::glue("fig/category_means_30763.png"), height = 5, width = 7, dpi = 640)
