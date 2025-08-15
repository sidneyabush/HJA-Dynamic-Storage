# 01_prep_PT_methods.R: Calculate Priestley-Taylor ET methods and export

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

rm(list = ls())

# --- SET DIRECTORIES ---
input_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/MET/data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"
plot_dir <- file.path(output_dir, "plots")
alpha_plot_dir <- file.path(plot_dir, "PT_alpha")
et_plot_dir <- file.path(plot_dir, "ET_methods_comparison")

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(alpha_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(et_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Import meteorological data
all_watersheds_data <- read_csv(file.path(input_dir, "watersheds_met_data_q.csv"))

# PT alpha calculation functions
calculate_alpha_zhang <- function(temp_celsius, rh_percent, pressure_kpa = 101.325) {
  q_specific <- 0.622 * 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3)) * (rh_percent/100) / (pressure_kpa - 0.378 * 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3)) * (rh_percent/100))
  delta <- 4098 * (0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))) / ((temp_celsius + 237.3)^2)
  gamma <- 0.067
  Bo <- 1 # For full formula, use your function
  alpha <- (1 + Bo) * (delta / (delta + gamma))
  return(alpha)
}
calculate_alpha_szilagyi <- function(temp_celsius) {
  alpha <- -3.89e-6 * temp_celsius^3 + 4.78e-4 * temp_celsius^2 - 2.54e-2 * temp_celsius + 1.64
  pmax(0.8, pmin(1.4, alpha))
}
calculate_et_pt <- function(alpha, net_radiation_wm2, temp_celsius, rh_percent) {
  net_radiation_mjm2d <- net_radiation_wm2 * 0.0864
  delta <- 4098 * 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3)) / ((temp_celsius + 237.3)^2)
  gamma <- 0.067
  lambda_v <- 2.501 - (0.002361 * temp_celsius)
  et_pt <- alpha * (delta / (delta + gamma)) * (net_radiation_mjm2d / lambda_v)
  pmax(0, et_pt)
}

# Calculate alphas
all_watersheds_data$alpha_zhang <- mapply(
  calculate_alpha_zhang,
  temp_celsius = all_watersheds_data$T_C,
  rh_percent = all_watersheds_data$RH_d_pct
)
all_watersheds_data$alpha_szilagyi <- sapply(all_watersheds_data$T_C, calculate_alpha_szilagyi)

# Calculate ET PT
all_watersheds_data$ET_PT_zhang <- calculate_et_pt(
  alpha = all_watersheds_data$alpha_zhang,
  net_radiation_wm2 = all_watersheds_data$NR_Wm2_d,
  temp_celsius = all_watersheds_data$T_C,
  rh_percent = all_watersheds_data$RH_d_pct
)
all_watersheds_data$ET_PT_szilagyi <- calculate_et_pt(
  alpha = all_watersheds_data$alpha_szilagyi,
  net_radiation_wm2 = all_watersheds_data$NR_Wm2_d,
  temp_celsius = all_watersheds_data$T_C,
  rh_percent = all_watersheds_data$RH_d_pct
)

# --- Export for next step ---
write_csv(all_watersheds_data, file.path(output_dir, "PT_ET_methods_timeseries.csv"))

# --- Alpha plots: long format ---
alpha_long <- all_watersheds_data %>%
  filter(DATE >= as.Date("1997-01-01")) %>%
  select(DATE, SITECODE, alpha_zhang, alpha_szilagyi) %>%
  pivot_longer(
    cols = starts_with("alpha_"),
    names_to = "Method",
    values_to = "Alpha"
  ) %>%
  mutate(Method = factor(Method,
                         levels = c("alpha_zhang", "alpha_szilagyi"),
                         labels = c("Zhang et al. (2024)", "Szilagyi et al. (2014)"))
  )

# Colors for alpha methods
alpha_colors <- c("Zhang et al. (2024)" = "#0173B2", "Szilagyi et al. (2014)" = "#DE8F05")

# Plot alpha through time for each site
for (site in unique(alpha_long$SITECODE)) {
  p <- ggplot(filter(alpha_long, SITECODE == site), aes(x = DATE, y = Alpha, color = Method)) +
    geom_line(size = 0.8, alpha = 0.98) +
    scale_color_manual(values = alpha_colors) +
    labs(title = paste0("Alpha (α) Through Time: ", site), x = "Date", y = "Alpha (α)", color = "Method") +
    theme_bw(base_size = 15) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  
  ggsave(file.path(alpha_plot_dir, paste0("Alpha_time_series_", site, ".png")), p, width = 10, height = 5)
}

# Grid plot for each method (all sites)
for (method in unique(alpha_long$Method)) {
  df <- alpha_long %>% filter(Method == method)
  method_name <- as.character(method)
  p_grid <- ggplot(df, aes(x = DATE, y = Alpha)) +
    geom_line(color = alpha_colors[method_name], size = 0.5, alpha = 0.8) +
    facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
    labs(title = paste(method_name, "Alpha Values - All Sites"),
         x = "Date", y = "Alpha (α)") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold"))
  
  ggsave(file.path(alpha_plot_dir, paste0("Alpha_", gsub(" ", "_", method_name), "_all_sites_grid.png")), 
         p_grid, width = 14, height = 10)
}
