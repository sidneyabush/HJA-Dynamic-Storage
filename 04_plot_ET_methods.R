# 04_plot_ET_methods.R: Plot ET methods and export figures

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

input_dir <- "../outputs"
plot_dir <- "../outputs/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

data <- read_csv(file.path(input_dir, "daily_water_balance_all_ET_methods_1997_present.csv"))

# Example: Long format for comparison plot
et_long <- data %>%
  filter(DATE >= as.Date("1997-01-01")) %>%
  select(DATE, SITECODE, ET_PT_zhang, ET_Hamon_uncalibrated, ET_Hamon_pt_zhang_cal) %>%
  pivot_longer(
    cols = c(ET_PT_zhang, ET_Hamon_uncalibrated, ET_Hamon_pt_zhang_cal),
    names_to = "Method",
    values_to = "ET_mm_day"
  ) %>%
  mutate(Method = factor(Method, levels = c("ET_PT_zhang", "ET_Hamon_pt_zhang_cal", "ET_Hamon_uncalibrated"),
                         labels = c("Zhang et al. (2024)", "Hamon-Zhang Calibrated", "Hamon Uncalibrated")))

# Plot
p <- ggplot(et_long, aes(x = as.Date(DATE), y = ET_mm_day, color = Method)) +
  geom_line(size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(title = "ET Methods Comparison", x = "Date", y = "ET (mm/day)", color = "Method") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")

ggsave(file.path(plot_dir, "ET_methods_comparison.png"), p, width = 14, height = 10)
