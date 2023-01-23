library(tidyr)
library(dplyr)
options(dplyr.summarise.inform = FALSE)

###################
# Set data variables
subj_id <- 30763

# Load results for single subject analysis
load(glue::glue("data/subj{subj_id}/dat_LPP_{subj_id}.RData"))

# Compute ests and CIs for RPAGP
dat_LPP |> group_by(name) |>
  summarize(lwr = quantile(value, 0.05),
            upr = quantile(value, 0.95),
            value = mean(value))

# Load results for bootstrap
load(glue::glue("data/subj{subj_id}/dat_bs_{subj_id}.RData"))

# Compute estimates and CIs for EMP
dat_bs |> filter(time > 300, time < 700) |>
  pivot_wider(names_from = category) |>
  group_by(b) |>
  summarize(H_L = mean(High - Low),
            H_N = mean(High - Neutral),
            L_N = mean(Low - Neutral)) |>
  ungroup() |>
  pivot_longer(cols = 2:4) |>
  group_by(name) |>
  summarize(lwr = quantile(value, 0.05),
            upr = quantile(value, 0.95),
            mean_value = mean(value))

