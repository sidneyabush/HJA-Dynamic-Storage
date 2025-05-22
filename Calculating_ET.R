library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

rm(list = ls())

# Set data directories
input_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/MET/data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Import all_watersheds_data
all_watersheds_data <- read_csv(file.path(input_dir, "watersheds_met_data_q.csv"))

# =============== THEORETICAL ALPHA FUNCTIONS (Zhang et al. 2016) ===============

calculate_specific_humidity <- function(temp_celsius, rh_percent, pressure_kpa = 101.325) {
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  e <- es * (rh_percent / 100)
  q <- 0.622 * e / (pressure_kpa - 0.378 * e)
  return(q)
}

calculate_delta <- function(temp_celsius) {
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  delta <- 4098 * es / ((temp_celsius + 237.3)^2)
  return(delta)
}

calculate_bowen_ratio_simplified <- function(temp_celsius, q_specific) {
  cp <- 1005  # J/kg/K
  lambda_v <- (2.501 - 0.002361 * temp_celsius) * 1e6  # J/kg
  epsilon_a <- lambda_v / cp
  k1 <- q_specific
  denominator <- 2 + (1 - k1 * epsilon_a)
  Bo <- (k1 * epsilon_a) / denominator
  Bo <- pmax(0.1, pmin(2.0, Bo))
  return(Bo)
}

calculate_alpha_theoretical <- function(temp_celsius, rh_percent, pressure_kpa = 101.325) {
  q_specific <- calculate_specific_humidity(temp_celsius, rh_percent, pressure_kpa)
  delta <- calculate_delta(temp_celsius)
  gamma <- 0.067  # kPa/°C at sea level
  Bo <- calculate_bowen_ratio_simplified(temp_celsius, q_specific)
  alpha <- (1 + Bo) * (delta / (delta + gamma))
  alpha <- pmax(0.9, pmin(1.4, alpha))
  return(alpha)
}

# =========== SZILAGYI et al. (2014) α(T) FUNCTION ===========
calculate_alpha_szilagyi <- function(temp_celsius) {
  alpha <- -3.89e-6 * temp_celsius^3 +
    4.78e-4 * temp_celsius^2 -
    2.54e-2 * temp_celsius + 1.64
  alpha <- pmax(0.8, pmin(1.4, alpha))  # enforce reasonable bounds
  return(alpha)
}

# ==================== PRIESTLEY-TAYLOR ET ====================
calculate_et_pt <- function(alpha, net_radiation_wm2, temp_celsius, rh_percent) {
  G <- 0
  net_radiation_mjm2d <- net_radiation_wm2 * 0.0864
  G_mjm2d <- G * 0.0864
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  delta <- 4098 * es / ((temp_celsius + 237.3)^2)
  gamma <- 0.067
  lambda_v <- 2.501 - (0.002361 * temp_celsius)
  available_energy <- net_radiation_mjm2d - G_mjm2d
  et_pt <- alpha * (delta / (delta + gamma)) * (available_energy / lambda_v)
  return(pmax(0, et_pt))
}

# ========== WEEKLY AVERAGING FOR THEORETICAL ALPHA ==========
calculate_weekly_theoretical_alpha <- function(data) {
  data$DATE <- as.Date(data$DATE)
  data$week_start <- floor_date(data$DATE, "week")
  weekly_stats <- data %>%
    group_by(SITECODE, week_start) %>%
    summarise(
      T_C_weekly = mean(T_C, na.rm = TRUE),
      RH_weekly = mean(RH_d_pct, na.rm = TRUE),
      n_days = sum(!is.na(T_C) & !is.na(RH_d_pct)),
      .groups = 'drop'
    ) %>%
    filter(n_days >= 4)
  weekly_stats$alpha_theoretical_weekly <- mapply(
    calculate_alpha_theoretical,
    temp_celsius = weekly_stats$T_C_weekly,
    rh_percent = weekly_stats$RH_weekly
  )
  data <- data %>%
    left_join(
      weekly_stats %>% 
        select(SITECODE, week_start, alpha_theoretical_weekly),
      by = c("SITECODE", "week_start")
    )
  return(data)
}

# ==================== MAIN WORKFLOW ====================
results <- all_watersheds_data

# Calculate weekly theoretical alpha
results <- calculate_weekly_theoretical_alpha(results)

# Calculate alpha_szilagyi
results$alpha_szilagyi <- sapply(results$T_C, calculate_alpha_szilagyi)

# Calculate ET for all methods
results$ET_PT_alpha_1.26 <- calculate_et_pt(
  alpha = 1.26,
  net_radiation_wm2 = results$NR_Wm2_d,
  temp_celsius = results$T_C,
  rh_percent = results$RH_d_pct
)

results$ET_PT_alpha_0.9 <- calculate_et_pt(
  alpha = 0.9,
  net_radiation_wm2 = results$NR_Wm2_d,
  temp_celsius = results$T_C,
  rh_percent = results$RH_d_pct
)

results$ET_PT_theoretical <- calculate_et_pt(
  alpha = results$alpha_theoretical_weekly,
  net_radiation_wm2 = results$NR_Wm2_d,
  temp_celsius = results$T_C,
  rh_percent = results$RH_d_pct
)

results$ET_PT_szilagyi <- calculate_et_pt(
  alpha = results$alpha_szilagyi,
  net_radiation_wm2 = results$NR_Wm2_d,
  temp_celsius = results$T_C,
  rh_percent = results$RH_d_pct
)

# Filter for complete cases across all four ET estimates
results_complete <- results %>%
  filter(complete.cases(select(., ET_PT_alpha_1.26, ET_PT_alpha_0.9, ET_PT_theoretical, ET_PT_szilagyi)))

# Pivot to long format for plotting
et_long <- results_complete %>%
  select(DATE, SITECODE, ET_PT_alpha_1.26, ET_PT_alpha_0.9, ET_PT_theoretical, ET_PT_szilagyi) %>%
  pivot_longer(
    cols = starts_with("ET_PT_"),
    names_to = "Method",
    values_to = "ET_mm_day"
  ) %>%
  mutate(Method = case_when(
    Method == "ET_PT_alpha_1.26" ~ "Fixed α = 1.26",
    Method == "ET_PT_alpha_0.9" ~ "Fixed α = 0.9",
    Method == "ET_PT_theoretical" ~ "Theoretical α (weekly)",
    Method == "ET_PT_szilagyi" ~ "Szilagyi et al. (2014) α(T)"
  ))

# ==================== PLOTTING ====================
p_et_facet <- ggplot(et_long, aes(x = DATE, y = ET_mm_day, color = Method)) +
  geom_line(alpha = 0.9, size = 0.6) +
  labs(
    title = "Daily ET Estimates by Method for Each Site",
    x = "Date",
    y = "ET (mm/day)",
    color = "Method"
  ) +
  facet_wrap(~ SITECODE, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "ET_facetgrid_all_methods_szilagyi.png"), p_et_facet, width = 14, height = 9)

cat("\n========== DAILY ET FACET GRID PLOT SAVED ==========\n")
cat("\nFigure saved to:", file.path(output_dir, "ET_facetgrid_all_methods_szilagyi.png"), "\n")
cat("\n====================================================\n")
