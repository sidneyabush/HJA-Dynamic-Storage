library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

rm(list = ls())

# Set data directories
input_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/MET/data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"
plot_dir <- file.path(output_dir, "plots")
alpha_plot_dir <- file.path(plot_dir, "PT_alpha")
et_plot_dir <- file.path(plot_dir, "ET_methods_comparison")

# Create all directories
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(alpha_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(et_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Import all_watersheds_data
all_watersheds_data <- read_csv(file.path(input_dir, "watersheds_met_data_q.csv"))

# HAMON COEFFICIENT CALIBRATION FUNCTIONS
calibrate_hamon_coefficient <- function(data, pt_method = "ET_PT_zhang", 
                                        calibration_period = c("2013-01-01", "2019-12-31"),
                                        by_day = FALSE, by_week = TRUE) {
  
  cal_data <- data %>%
    filter(DATE >= as.Date(calibration_period[1]) & 
             DATE <= as.Date(calibration_period[2])) %>%
    filter(!is.na(.data[[pt_method]]) & !is.na(ET_Hamon_uncalibrated))
  
  if (by_day) {
    # Daily calibration coefficients
    daily_coefs <- cal_data %>%
      group_by(SITECODE, DATE) %>%
      summarise(
        coefficient = .data[[pt_method]] / ET_Hamon_uncalibrated,
        .groups = "drop"
      ) %>%
      filter(!is.na(coefficient) & is.finite(coefficient))
    
    return(daily_coefs)
    
  } else if (by_week) {
    # Weekly calibration coefficients
    weekly_coefs <- cal_data %>%
      mutate(week_start = floor_date(DATE, "week")) %>%
      group_by(SITECODE, week_start) %>%
      summarise(
        coefficient = mean(.data[[pt_method]] / ET_Hamon_uncalibrated, na.rm = TRUE),
        n_obs = n(),
        rmse_uncalibrated = sqrt(mean((.data[[pt_method]] - ET_Hamon_uncalibrated)^2, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      filter(n_obs >= 4)
    
    return(weekly_coefs)
    
  } else {
    # Annual calibration coefficient
    annual_coefs <- cal_data %>%
      group_by(SITECODE) %>%
      summarise(
        coefficient = mean(.data[[pt_method]] / ET_Hamon_uncalibrated, na.rm = TRUE),
        n_obs = n(),
        rmse_uncalibrated = sqrt(mean((.data[[pt_method]] - ET_Hamon_uncalibrated)^2, na.rm = TRUE)),
        .groups = "drop"
      )
    
    return(annual_coefs)
  }
}

apply_hamon_calibration <- function(data, calibration_coefs, by_day = FALSE, by_week = TRUE) {
  
  if (by_day) {
    # Apply daily coefficients
    data <- data %>%
      left_join(calibration_coefs %>% select(SITECODE, DATE, coefficient), 
                by = c("SITECODE", "DATE")) %>%
      mutate(
        ET_Hamon_calibrated = ET_Hamon_uncalibrated * coefficient,
        ET_Hamon_calibrated = ifelse(is.na(coefficient), ET_Hamon_uncalibrated * 0.1651, ET_Hamon_calibrated)
      ) %>%
      select(-coefficient)
    
  } else if (by_week) {
    # Apply weekly coefficients
    data <- data %>%
      mutate(week_start = floor_date(DATE, "week")) %>%
      left_join(calibration_coefs %>% select(SITECODE, week_start, coefficient), 
                by = c("SITECODE", "week_start")) %>%
      mutate(
        ET_Hamon_calibrated = ET_Hamon_uncalibrated * coefficient,
        ET_Hamon_calibrated = ifelse(is.na(coefficient), ET_Hamon_uncalibrated * 0.1651, ET_Hamon_calibrated)
      ) %>%
      select(-week_start, -coefficient)
    
  } else {
    # Apply annual coefficient
    data <- data %>%
      left_join(calibration_coefs %>% select(SITECODE, coefficient), 
                by = "SITECODE") %>%
      mutate(
        ET_Hamon_calibrated = ET_Hamon_uncalibrated * coefficient,
        ET_Hamon_calibrated = ifelse(is.na(coefficient), ET_Hamon_uncalibrated * 0.1651, ET_Hamon_calibrated)
      ) %>%
      select(-coefficient)
  }
  
  return(data)
}

# Function to create linear interpolation relationship between Zhang and Hamon
create_linear_interpolation <- function(data, calibration_period = c("2013-01-01", "2019-12-31")) {
  
  # Filter to calibration period where both methods exist
  cal_data <- data %>%
    filter(DATE >= as.Date(calibration_period[1]) & 
             DATE <= as.Date(calibration_period[2])) %>%
    filter(!is.na(ET_PT_zhang) & !is.na(ET_Hamon_uncalibrated))
  
  # Create linear model for each site: Hamon_uncalibrated = slope * Zhang + intercept
  linear_models <- cal_data %>%
    group_by(SITECODE) %>%
    summarise(
      model_data = list(lm(ET_Hamon_uncalibrated ~ ET_PT_zhang, data = cur_data())),
      slope = sapply(model_data, function(x) coef(x)[2]),
      intercept = sapply(model_data, function(x) coef(x)[1]),
      r_squared = sapply(model_data, function(x) summary(x)$r.squared),
      .groups = "drop"
    ) %>%
    select(-model_data)
  
  return(linear_models)
}

# Function to apply linear interpolation to estimate Hamon from Zhang
apply_linear_interpolation <- function(data, linear_models) {
  
  data <- data %>%
    left_join(linear_models, by = "SITECODE") %>%
    mutate(
      # Calculate Hamon from Zhang using linear relationship
      ET_Hamon_linear = slope * ET_PT_zhang + intercept,
      # Use original uncalibrated Hamon where Zhang doesn't exist, linear where it does
      ET_Hamon_linear = ifelse(!is.na(ET_PT_zhang), ET_Hamon_linear, ET_Hamon_uncalibrated),
      # Clean up temporary variables
      ET_Hamon_linear = pmax(0, ET_Hamon_linear, na.rm = TRUE)
    ) %>%
    select(-slope, -intercept, -r_squared)
  
  return(data)
}

# HAMON (1963) ET CALCULATION
calculate_et_hamon <- function(temp_celsius, date_vector, latitude = 44.24, coefficient = 0.1651) {
  day_of_year <- yday(as.Date(date_vector))
  declination <- 0.4093 * sin((2 * pi * day_of_year / 365) - 1.405)
  lat_rad <- latitude * pi / 180
  cos_sunset <- -tan(lat_rad) * tan(declination)
  cos_sunset <- pmax(-1, pmin(1, cos_sunset))
  sunset_hour_angle <- acos(cos_sunset)
  daylight_hours <- (24 / pi) * sunset_hour_angle
  es_kpa <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  temp_kelvin <- temp_celsius + 273.16
  vapor_density <- 216.7 * es_kpa / temp_kelvin
  et_hamon <- coefficient * (daylight_hours / 12) * vapor_density
  return(pmax(0, et_hamon))
}

# ALPHA FUNCTIONS
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
  cp <- 1005
  lambda_v <- (2.501 - 0.002361 * temp_celsius) * 1e6
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
  gamma <- 0.067
  Bo <- calculate_bowen_ratio_simplified(temp_celsius, q_specific)
  alpha <- (1 + Bo) * (delta / (delta + gamma))
  return(alpha)
}

calculate_alpha_szilagyi <- function(temp_celsius) {
  alpha <- -3.89e-6 * temp_celsius^3 +
    4.78e-4 * temp_celsius^2 -
    2.54e-2 * temp_celsius + 1.64
  alpha <- pmax(0.8, pmin(1.4, alpha))
  return(alpha)
}

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

# MAIN WORKFLOW
results <- all_watersheds_data
results <- calculate_weekly_theoretical_alpha(results)
results$alpha_szilagyi <- sapply(results$T_C, calculate_alpha_szilagyi)

# Calculate PT methods
results$ET_PT_zhang <- calculate_et_pt(
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

# Calculate Hamon
study_latitude <- (44.20127400 + 44.28226000) / 2

results$ET_Hamon_uncalibrated <- calculate_et_hamon(
  temp_celsius = results$T_C,
  date_vector = results$DATE,
  latitude = study_latitude,
  coefficient = 0.1651
)

# Calibrate Hamon
daily_calibration <- calibrate_hamon_coefficient(
  data = results,
  pt_method = "ET_PT_zhang",
  calibration_period = c("2013-01-01", "2019-12-31"),
  by_day = TRUE
)

sites_with_daily_cal <- daily_calibration %>%
  group_by(SITECODE) %>%
  summarise(days_available = n()) %>%
  filter(days_available >= 100)

if (nrow(sites_with_daily_cal) > 0) {
  results <- apply_hamon_calibration(
    data = results,
    calibration_coefs = daily_calibration,
    by_day = TRUE
  )
} else {
  weekly_calibration <- calibrate_hamon_coefficient(
    data = results,
    pt_method = "ET_PT_zhang",
    calibration_period = c("2013-01-01", "2019-12-31"),
    by_week = TRUE
  )
  
  sites_with_weekly_cal <- weekly_calibration %>%
    group_by(SITECODE) %>%
    summarise(weeks_available = n()) %>%
    filter(weeks_available >= 26)
  
  if (nrow(sites_with_weekly_cal) > 0) {
    results <- apply_hamon_calibration(
      data = results,
      calibration_coefs = weekly_calibration,
      by_day = FALSE,
      by_week = TRUE
    )
  } else {
    annual_calibration <- calibrate_hamon_coefficient(
      data = results,
      pt_method = "ET_PT_zhang",
      calibration_period = c("2013-01-01", "2019-12-31"),
      by_day = FALSE,
      by_week = FALSE
    )
    results <- apply_hamon_calibration(
      data = results,
      calibration_coefs = annual_calibration,
      by_day = FALSE,
      by_week = FALSE
    )
  }
}

results$ET_Hamon <- results$ET_Hamon_calibrated

# Create linear interpolation relationship between Zhang and Hamon uncalibrated
linear_models <- create_linear_interpolation(results, calibration_period = c("2013-01-01", "2019-12-31"))

# Apply linear interpolation to create Hamon linear estimates
results <- apply_linear_interpolation(results, linear_models)

# Filter for complete cases
results_complete <- results %>%
  filter(complete.cases(select(., ET_PT_zhang, ET_PT_szilagyi, ET_Hamon)))

# Prepare data for plotting - Include Hamon calibrated, uncalibrated, and linear
et_long <- results_complete %>%
  select(DATE, SITECODE, ET_PT_zhang, ET_PT_szilagyi, ET_Hamon, ET_Hamon_uncalibrated, ET_Hamon_linear) %>%
  pivot_longer(
    cols = c(ET_PT_zhang, ET_PT_szilagyi, ET_Hamon, ET_Hamon_uncalibrated, ET_Hamon_linear),
    names_to = "Method",
    values_to = "ET_mm_day"
  ) %>%
  mutate(Method = factor(Method, 
                         levels = c("ET_PT_zhang", "ET_PT_szilagyi", "ET_Hamon", "ET_Hamon_uncalibrated", "ET_Hamon_linear"),
                         labels = c("Zhang et al. (2024)", "Szilagyi et al. (2014)", "Hamon Calibrated", "Hamon Uncalibrated", "Hamon Linear")
  ))

alpha_long <- results_complete %>%
  mutate(
    alpha_zhang = alpha_theoretical_weekly,
    alpha_szilagyi = alpha_szilagyi
  ) %>%
  select(DATE, SITECODE, alpha_zhang, alpha_szilagyi) %>%
  pivot_longer(
    cols = starts_with("alpha_"),
    names_to = "Method",
    values_to = "Alpha"
  ) %>%
  mutate(Method = factor(Method,
                         levels = c("alpha_zhang", "alpha_szilagyi"),
                         labels = c("Zhang et al. (2024)", "Szilagyi et al. (2014)")
  ))

# PLOTTING

# Colors - Colorblind friendly palette with 5 methods
method_colors <- c("Zhang et al. (2024)" = "#0173B2",        # Blue
                   "Szilagyi et al. (2014)" = "#DE8F05",     # Orange  
                   "Hamon Calibrated" = "#029E73",           # Green
                   "Hamon Uncalibrated" = "#CC78BC",         # Pink
                   "Hamon Linear" = "#56B4E9")               # Light Blue
alpha_colors <- c("Zhang et al. (2024)" = "#0173B2", "Szilagyi et al. (2014)" = "#DE8F05")
comparison_colors <- c("Zhang et al. (2024)" = "#0173B2", "Hamon Calibrated" = "#029E73", "Hamon Uncalibrated" = "#CC78BC")

# Individual site plots - Now with 5 methods
site_list <- unique(et_long$SITECODE)

for (site in site_list) {
  p <- ggplot(filter(et_long, SITECODE == site), aes(x = DATE, y = ET_mm_day, color = Method)) +
    geom_line(size = 0.7, alpha = 0.93) +  
    scale_color_manual(values = method_colors) +
    labs(title = paste0("Daily ET Comparison at ", site), x = "Date", y = "ET (mm/day)", color = "Method") +
    theme_bw(base_size = 15) +  
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank()) +
    guides(color = guide_legend(ncol = 2))
  
  ggsave(file.path(et_plot_dir, paste0("ET_methods_comparison_", site, ".png")), p, width = 12, height = 8)
}

# Grid plots for each method
zhang_data <- et_long %>% filter(Method == "Zhang et al. (2024)")
p_zhang_grid <- ggplot(zhang_data, aes(x = DATE, y = ET_mm_day)) +
  geom_line(color = method_colors["Zhang et al. (2024)"], size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(title = "Zhang et al. (2024) Method - All Sites", subtitle = "Priestley-Taylor (2013-2019)",
       x = "Date", y = "ET (mm/day)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

ggsave(file.path(et_plot_dir, "ET_Zhang_all_sites_grid.png"), p_zhang_grid, width = 14, height = 10)

szilagyi_data <- et_long %>% filter(Method == "Szilagyi et al. (2014)")
p_szilagyi_grid <- ggplot(szilagyi_data, aes(x = DATE, y = ET_mm_day)) +
  geom_line(color = method_colors["Szilagyi et al. (2014)"], size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(title = "Szilagyi et al. (2014) Method - All Sites", subtitle = "Priestley-Taylor (2013-2019)", 
       x = "Date", y = "ET (mm/day)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

ggsave(file.path(et_plot_dir, "ET_Szilagyi_all_sites_grid.png"), p_szilagyi_grid, width = 14, height = 10)

hamon_data <- et_long %>% filter(Method == "Hamon Calibrated")

if (nrow(hamon_data) > 0) {
  p_hamon_grid <- ggplot(hamon_data, aes(x = DATE, y = ET_mm_day)) +
    geom_line(color = method_colors["Hamon Calibrated"], size = 0.5, alpha = 0.8) +
    facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
    labs(title = "Hamon Calibrated Method - All Sites (2013-2019)", 
         subtitle = "Temperature-based PET (calibrated to PT-Zhang)", x = "Date", y = "ET (mm/day)") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11))
  
  ggsave(file.path(et_plot_dir, "ET_Hamon_all_sites_grid_2013_2019.png"), p_hamon_grid, width = 14, height = 10)
}

# Full time series Hamon plot
hamon_full_data <- results %>%
  filter(DATE >= as.Date("1997-01-01") & !is.na(ET_Hamon)) %>%
  select(DATE, SITECODE, ET_Hamon)

p_hamon_full_grid <- ggplot(hamon_full_data, aes(x = DATE, y = ET_Hamon)) +
  geom_line(color = method_colors["Hamon Calibrated"], size = 0.4, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(title = "Hamon Calibrated Method - All Sites (Full Time Series)",
       subtitle = "Temperature-based PET (1997-present, calibrated to PT-Zhang 2013-2019)",
       x = "Date", y = "ET (mm/day)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

ggsave(file.path(et_plot_dir, "ET_Hamon_all_sites_grid_full_timeseries.png"), p_hamon_full_grid, width = 16, height = 12)

# COMPARISON PLOTS: Zhang vs Hamon Calibrated vs Hamon Uncalibrated

# Overlap period comparison (2013-2019)
comparison_overlap_data <- results_complete %>%
  select(DATE, SITECODE, ET_PT_zhang, ET_Hamon, ET_Hamon_uncalibrated) %>%
  pivot_longer(cols = c(ET_PT_zhang, ET_Hamon, ET_Hamon_uncalibrated),
               names_to = "Method", values_to = "ET_mm_day") %>%
  filter(!is.na(ET_mm_day)) %>%
  mutate(Method = factor(Method,
                         levels = c("ET_PT_zhang", "ET_Hamon", "ET_Hamon_uncalibrated"),
                         labels = c("Zhang et al. (2024)", "Hamon Calibrated", "Hamon Uncalibrated")))

p_comparison_overlap <- ggplot(comparison_overlap_data, aes(x = DATE, y = ET_mm_day, color = Method)) +
  geom_line(size = 0.4, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  scale_color_manual(values = comparison_colors) +
  labs(title = "ET Method Comparison - All Sites (2013-2019)",
       subtitle = "Zhang vs Hamon Calibrated vs Hamon Uncalibrated (overlap period)",
       x = "Date", y = "ET (mm/day)", color = "Method") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom")

ggsave(file.path(et_plot_dir, "ET_comparison_Zhang_Hamon_all_sites_2013_2019.png"), p_comparison_overlap, width = 16, height = 14)

# Full time series comparison (1997-present)
comparison_full_data <- results %>%
  filter(DATE >= as.Date("1997-01-01")) %>%
  select(DATE, SITECODE, ET_PT_zhang, ET_Hamon, ET_Hamon_uncalibrated) %>%
  pivot_longer(cols = c(ET_PT_zhang, ET_Hamon, ET_Hamon_uncalibrated),
               names_to = "Method", values_to = "ET_mm_day") %>%
  filter(!is.na(ET_mm_day)) %>%
  mutate(Method = factor(Method,
                         levels = c("ET_PT_zhang", "ET_Hamon", "ET_Hamon_uncalibrated"),
                         labels = c("Zhang et al. (2024)", "Hamon Calibrated", "Hamon Uncalibrated")))

p_comparison_full <- ggplot(comparison_full_data, aes(x = DATE, y = ET_mm_day, color = Method)) +
  geom_line(size = 0.3, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  scale_color_manual(values = comparison_colors) +
  labs(title = "ET Method Comparison - All Sites (Full Time Series)",
       subtitle = "Zhang vs Hamon Calibrated vs Hamon Uncalibrated (1997-present, where data available)",
       x = "Date", y = "ET (mm/day)", color = "Method") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom")

ggsave(file.path(et_plot_dir, "ET_comparison_Zhang_Hamon_all_sites_full_timeseries.png"), p_comparison_full, width = 18, height = 14)

# Alpha plots
for (site in unique(alpha_long$SITECODE)) {
  p <- ggplot(filter(alpha_long, SITECODE == site), aes(x = DATE, y = Alpha, color = Method)) +
    geom_line(size = 0.8, alpha = 0.98) +
    labs(title = paste0("Alpha (α) Through Time: ", site), x = "Date", y = "Alpha (α)", color = "Method") +
    theme_bw(base_size = 15) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  
  ggsave(file.path(alpha_plot_dir, paste0("Alpha_time_series_", site, ".png")), p, width = 10, height = 5)
}

# Alpha grid plots
zhang_alpha_data <- alpha_long %>% filter(Method == "Zhang et al. (2024)")
p_zhang_alpha_grid <- ggplot(zhang_alpha_data, aes(x = DATE, y = Alpha)) +
  geom_line(color = alpha_colors["Zhang et al. (2024)"], size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(title = "Zhang et al. (2024) Alpha Values - All Sites", 
       subtitle = "Temperature and humidity dependent alpha coefficient", x = "Date", y = "Alpha (α)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

ggsave(file.path(alpha_plot_dir, "Alpha_Zhang_all_sites_grid.png"), p_zhang_alpha_grid, width = 14, height = 10)

szilagyi_alpha_data <- alpha_long %>% filter(Method == "Szilagyi et al. (2014)")
p_szilagyi_alpha_grid <- ggplot(szilagyi_alpha_data, aes(x = DATE, y = Alpha)) +
  geom_line(color = alpha_colors["Szilagyi et al. (2014)"], size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(title = "Szilagyi et al. (2014) Alpha Values - All Sites",
       subtitle = "Temperature dependent alpha coefficient", x = "Date", y = "Alpha (α)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

ggsave(file.path(alpha_plot_dir, "Alpha_Szilagyi_all_sites_grid.png"), p_szilagyi_alpha_grid, width = 14, height = 10)

# EXPORT SINGLE WATER BALANCE CSV
water_balance_all_methods <- results %>%
  filter(DATE >= as.Date("1997-01-01")) %>%
  transmute(
    DATE,
    SITECODE,
    P_mm_day = P_mm_d,
    Q_mm_day = Q_mm_d,
    ET_hamon_calibrated = ET_Hamon,
    ET_hamon_uncalibrated = ET_Hamon_uncalibrated,
    ET_hamon_linear = ET_Hamon_linear,
    ET_pt_zhang = ET_PT_zhang,
    ET_pt_szilagyi = ET_PT_szilagyi
  )

write_csv(water_balance_all_methods, file.path(output_dir, "daily_water_balance_all_ET_methods_1997_present.csv"))