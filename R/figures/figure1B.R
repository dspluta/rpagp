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

# EEG scalp map for plots
img <- readPNG(source = "figures/EEG_scalp_map.PNG")
g <- rasterGrob(img, interpolate=TRUE)

############################
#### Fig 1 ####

# Subj 30763 High Category Plot
selected_category <- "High"
ggplot(dat %>% filter(SubjectID == subj_id, category == selected_category) |> group_by(category, time) %>% summarize(value = mean(value))) +
  geom_line(data = dat |> filter(SubjectID == subj_id, category == "High"),
            aes(x = time, y = value, group = factor(prognum)), alpha = 0.5) +
  geom_line(aes(x = time, y = value, color = category), size = 2.9) +
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
ggsave(glue::glue("figures/figure1/figure1B_subject_{subj_id}_means_{selected_category}.png"), height = 5, width = 8, dpi = 640)

# Subj 30763 Low Category Plot
selected_category <- "Low"
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
ggsave(glue::glue("figures/figure1/figure1B_subject_{subj_id}_means_{selected_category}.png"), height = 5, width = 8, dpi = 640)

# Subj 30763 Neutral Category Plot
selected_category <- "Neutral"
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
ggsave(glue::glue("figures/figure1/figure1B_subject_{subj_id}_means_{selected_category}.png"), height = 5, width = 8, dpi = 640)

# Subj 30763 Category Means
selected_category <- "Low"
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
ggsave(glue::glue("figures/figure1/figure1B_subject_{subj_id}_category_means.png"), height = 5, width = 8, dpi = 640)
