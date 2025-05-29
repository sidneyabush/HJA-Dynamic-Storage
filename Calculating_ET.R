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

# THEORETICAL ALPHA FUNCTIONS (Zhang et al. 2024): 
# https://hess.copernicus.org/articles/28/4349/2024/hess-28-4349-2024.html 
calculate_specific_humidity <- function(temp_celsius, rh_percent, pressure_kpa = 101.325) {
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  e <- es * (rh_percent / 100)
  q <- 0.622 * e / (pressure_kpa - 0.378 * e)
  return(q)
}

#  Equation 17 
calculate_delta <- function(temp_celsius) {
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  delta <- 4098 * es / ((temp_celsius + 237.3)^2)
  return(delta)
}

# Equation 9
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

# Equation 12
calculate_alpha_theoretical <- function(temp_celsius, rh_percent, pressure_kpa = 101.325) {
  q_specific <- calculate_specific_humidity(temp_celsius, rh_percent, pressure_kpa)
  delta <- calculate_delta(temp_celsius)
  gamma <- 0.067  # kPa/°C at sea level
  Bo <- calculate_bowen_ratio_simplified(temp_celsius, q_specific)
  alpha <- (1 + Bo) * (delta / (delta + gamma))
  # alpha <- pmax(0.9, pmin(1.4, alpha)) 
  return(alpha)
}

# This is a formula I tried, but am no longer using. Using Zhang method.
# SZILAGYI et al. (2014) alpha(T) FUNCTION
calculate_alpha_szilagyi <- function(temp_celsius) {
  alpha <- -3.89e-6 * temp_celsius^3 +
    4.78e-4 * temp_celsius^2 -
    2.54e-2 * temp_celsius + 1.64
  alpha <- pmax(0.8, pmin(1.4, alpha))
  return(alpha)
}

# HAMON COEFFICIENT CALIBRATION FUNCTIONS
# Calibrate Hamon coefficient using Priestley-Taylor as reference

# Function to calibrate Hamon coefficient against PT reference
calibrate_hamon_coefficient <- function(data, pt_method = "ET_PT_zhang", 
                                        calibration_period = c("2013-01-01", "2019-12-31"),
                                        by_month = TRUE) {
  
  # Filter to calibration period
  cal_data <- data %>%
    filter(DATE >= as.Date(calibration_period[1]) & 
             DATE <= as.Date(calibration_period[2])) %>%
    filter(!is.na(.data[[pt_method]]) & !is.na(ET_Hamon_uncalibrated))
  
  if (by_month) {
    # Monthly calibration coefficients
    monthly_coefs <- cal_data %>%
      mutate(month = month(DATE)) %>%
      group_by(SITECODE, month) %>%
      summarise(
        coefficient = mean(.data[[pt_method]] / ET_Hamon_uncalibrated, na.rm = TRUE),
        n_obs = n(),
        rmse_uncalibrated = sqrt(mean((.data[[pt_method]] - ET_Hamon_uncalibrated)^2, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      filter(n_obs >= 10)  # Require at least 10 observations per month
    
    return(monthly_coefs)
    
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

# Function to apply calibrated coefficients
apply_hamon_calibration <- function(data, calibration_coefs, by_month = TRUE) {
  
  if (by_month) {
    # Apply monthly coefficients
    data <- data %>%
      mutate(month = month(DATE)) %>%
      left_join(calibration_coefs %>% select(SITECODE, month, coefficient), 
                by = c("SITECODE", "month")) %>%
      mutate(
        ET_Hamon_calibrated = ET_Hamon_uncalibrated * coefficient,
        ET_Hamon_calibrated = ifelse(is.na(coefficient), ET_Hamon_uncalibrated * 0.1651, ET_Hamon_calibrated)
      ) %>%
      select(-month, -coefficient)
    
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
# Based on: https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/evaporation-and-transpiration/hamon-method
# ET_o = c * (N/12) * P_t
# where:
# - c = coefficient (0.1651 mm/g/m³ default)
# - N = daylight hours 
# - P_t = saturated vapor density (g/m³)
calculate_et_hamon <- function(temp_celsius, date_vector, latitude = 44.24, coefficient = 0.1651) {
  # Calculate day of year
  day_of_year <- yday(as.Date(date_vector))
  
  # Calculate solar declination (radians) - Allen et al. 1998
  declination <- 0.4093 * sin((2 * pi * day_of_year / 365) - 1.405)
  
  # Convert latitude to radians
  lat_rad <- latitude * pi / 180
  
  # Calculate sunset hour angle (radians) - Allen et al. 1998
  cos_sunset <- -tan(lat_rad) * tan(declination)
  cos_sunset <- pmax(-1, pmin(1, cos_sunset))  # Constrain to [-1, 1]
  sunset_hour_angle <- acos(cos_sunset)
  
  # Calculate daylight hours N - Allen et al. 1998
  # N = (24/π) * ω_s
  daylight_hours <- (24 / pi) * sunset_hour_angle  # This is N
  
  # Calculate saturated vapor pressure es(T) - Allen et al. 1998
  # e_s(T) = 0.6108 * exp[17.27*T / (T + 237.3)] (kPa)
  es_kpa <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  
  # Calculate saturated vapor density P_t - Wiederhold 1997
  # P_t = 216.7 * e_s(T) / (T + 273.16) (g/m³)
  temp_kelvin <- temp_celsius + 273.16
  vapor_density <- 216.7 * es_kpa / temp_kelvin
  
  # Hamon equation - HEC-HMS implementation
  # ET_o = c * (N/12) * P_t (mm/day)
  et_hamon <- coefficient * (daylight_hours / 12) * vapor_density
  
  return(pmax(0, et_hamon))
}

# PRIESTLEY-TAYLOR ET
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

# WEEKLY AVERAGING FOR THEORETICAL ALPHA (ZHANG)
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

# Calculate all ET methods
results$ET_PT_fixed_1.26 <- calculate_et_pt(
  alpha = 1.26,
  net_radiation_wm2 = results$NR_Wm2_d,
  temp_celsius = results$T_C,
  rh_percent = results$RH_d_pct
)

results$ET_PT_fixed_0.9 <- calculate_et_pt(
  alpha = 0.9,
  net_radiation_wm2 = results$NR_Wm2_d,
  temp_celsius = results$T_C,
  rh_percent = results$RH_d_pct
)

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

# ADD HAMON METHOD - Using study area coordinates
# Latitude: Center of your bounding box (44.20127400 to 44.28226000)
study_latitude <- (44.20127400 + 44.28226000) / 2  # 44.24176700

# Calculate uncalibrated Hamon first
results$ET_Hamon_uncalibrated <- calculate_et_hamon(
  temp_celsius = results$T_C,
  date_vector = results$DATE,
  latitude = study_latitude,  # ~44.24° N
  coefficient = 0.1651  # HEC-HMS default (mm/g/m³)
)

# CALIBRATE HAMON COEFFICIENT using Priestley-Taylor Zhang method as reference
print("Calibrating Hamon coefficient using PT-Zhang (2013-2019)...")

# Try monthly calibration first
monthly_calibration <- calibrate_hamon_coefficient(
  data = results,
  pt_method = "ET_PT_zhang",
  calibration_period = c("2013-01-01", "2019-12-31"),
  by_month = TRUE
)

# Check if we have sufficient data for monthly calibration
sites_with_monthly_cal <- monthly_calibration %>%
  group_by(SITECODE) %>%
  summarise(months_available = n()) %>%
  filter(months_available >= 6)  # At least 6 months of data

if (nrow(sites_with_monthly_cal) > 0) {
  print("Using monthly calibration coefficients")
  
  # Apply monthly calibration
  results <- apply_hamon_calibration(
    data = results,
    calibration_coefs = monthly_calibration,
    by_month = TRUE
  )
  
  # Print calibration summary
  print("Monthly calibration coefficients summary:")
  cal_summary <- monthly_calibration %>%
    group_by(month) %>%
    summarise(
      mean_coef = mean(coefficient, na.rm = TRUE),
      min_coef = min(coefficient, na.rm = TRUE),
      max_coef = max(coefficient, na.rm = TRUE),
      .groups = "drop"
    )
  print(cal_summary)
  
} else {
  print("Insufficient data for monthly calibration, using annual calibration")
  
  # Fall back to annual calibration
  annual_calibration <- calibrate_hamon_coefficient(
    data = results,
    pt_method = "ET_PT_zhang",
    calibration_period = c("2013-01-01", "2019-12-31"),
    by_month = FALSE
  )
  
  # Apply annual calibration
  results <- apply_hamon_calibration(
    data = results,
    calibration_coefs = annual_calibration,
    by_month = FALSE
  )
  
  print("Annual calibration coefficients:")
  print(annual_calibration)
}

# Set final calibrated Hamon ET
results$ET_Hamon <- results$ET_Hamon_calibrated

# Filter for complete cases (now including Hamon)
results_complete <- results %>%
  filter(complete.cases(select(., ET_PT_fixed_1.26, ET_PT_fixed_0.9, ET_PT_zhang, ET_PT_szilagyi, ET_Hamon)))

# Pivot to long format for ET (excluding fixed alpha methods)
et_long <- results_complete %>%
  select(DATE, SITECODE, ET_PT_zhang, ET_PT_szilagyi, ET_Hamon) %>%
  pivot_longer(
    cols = c(ET_PT_zhang, ET_PT_szilagyi, ET_Hamon),
    names_to = "Method",
    values_to = "ET_mm_day"
  ) %>%
  mutate(Method = factor(Method, 
                         levels = c("ET_PT_zhang", "ET_PT_szilagyi", "ET_Hamon"),
                         labels = c("Zhang et al. (2024)", "Szilagyi et al. (2014)", "Hamon (1963)")
  ))

# Add the alpha values (for plotting) - excluding fixed values
results_complete <- results_complete %>%
  mutate(
    alpha_zhang = alpha_theoretical_weekly,
    alpha_szilagyi = alpha_szilagyi
  )

# Pivot to long format for alpha (excluding fixed methods and Hamon)
alpha_long <- results_complete %>%
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

# Export one plot per SITECODE (now including Hamon)
site_list <- unique(et_long$SITECODE)

for (site in site_list) {
  p <- ggplot(filter(et_long, SITECODE == site),
              aes(x = DATE, y = ET_mm_day, color = Method)) +
    geom_line(size = 0.7, alpha = 0.93) +  
    labs(
      title = paste0("Daily ET Comparison at ", site),
      x = "Date",
      y = "ET (mm/day)",
      color = "Method"
    ) +
    theme_bw(base_size = 15) +  
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(ncol = 2))  # Organize legend in 2 columns due to more methods
  
  ggsave(file.path(plot_dir, paste0("ET_methods_comparison_", site, ".png")), p, width = 12, height = 7)
}

# Export alpha plots per SITECODE (now with 2 methods only) - to PT_alpha folder
for (site in unique(alpha_long$SITECODE)) {
  p <- ggplot(filter(alpha_long, SITECODE == site),
              aes(x = DATE, y = Alpha, color = Method)) +
    geom_line(size = 0.8, alpha = 0.98) +
    labs(
      title = paste0("Alpha (α) Through Time: ", site),
      x = "Date",
      y = "Alpha (α)",
      color = "Method"
    ) +
    theme_bw(base_size = 15) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(alpha_plot_dir, paste0("Alpha_time_series_", site, ".png")), p, width = 10, height = 5)
}

# CREATE ALPHA GRID PLOTS - to PT_alpha folder
alpha_colors <- c("Zhang et al. (2024)" = "#2E8B57",      # Sea Green
                  "Szilagyi et al. (2014)" = "#4169E1")   # Royal Blue

# Zhang alpha - All sites
zhang_alpha_data <- alpha_long %>% filter(Method == "Zhang et al. (2024)")

p_zhang_alpha_grid <- ggplot(zhang_alpha_data, aes(x = DATE, y = Alpha)) +
  geom_line(color = alpha_colors["Zhang et al. (2024)"], size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(
    title = "Zhang et al. (2024) Alpha (α) Values - All Sites",
    subtitle = "Temperature and humidity dependent alpha coefficient",
    x = "Date",
    y = "Alpha (α)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

ggsave(file.path(alpha_plot_dir, "Alpha_Zhang_all_sites_grid.png"), p_zhang_alpha_grid, width = 14, height = 10)

# Szilagyi alpha - All sites
szilagyi_alpha_data <- alpha_long %>% filter(Method == "Szilagyi et al. (2014)")

p_szilagyi_alpha_grid <- ggplot(szilagyi_alpha_data, aes(x = DATE, y = Alpha)) +
  geom_line(color = alpha_colors["Szilagyi et al. (2014)"], size = 0.5, alpha = 0.8) +
  facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
  labs(
    title = "Szilagyi et al. (2014) Alpha (α) Values - All Sites",
    subtitle = "Temperature dependent alpha coefficient", 
    x = "Date",
    y = "Alpha (α)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

ggsave(file.path(alpha_plot_dir, "Alpha_Szilagyi_all_sites_grid.png"), p_szilagyi_alpha_grid, width = 14, height = 10)

# Add year and month columns for grouping
et_long <- et_long %>%
  mutate(
    year = year(DATE),
    month = month(DATE, label = TRUE, abbr = TRUE)
  )

# Summarize mean daily ET by month, method, and site (now including Hamon)
et_monthly_summary <- et_long %>%
  group_by(SITECODE, Method, month) %>%
  summarise(
    mean_ET_mm_day = mean(ET_mm_day, na.rm = TRUE),
    sd_ET_mm_day = sd(ET_mm_day, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  )

# Preview result
print(head(et_monthly_summary, 15))

# Final dataset: Zhang method + Hamon, all variables retained
export_data <- results_complete %>%
  mutate(
    zhang_alpha = alpha_theoretical_weekly,      
    ET_PT_zhang = ET_PT_zhang,
    ET_Hamon = ET_Hamon
  ) %>%
  select(DATE, SITECODE, everything())

# Export to CSV
write_csv(export_data, file.path(output_dir, "daily_ET_watersheds_all_methods.csv"))

# ---- Export daily water balance file with both Zhang and Hamon ----
water_balance_export <- results_complete %>%
  transmute(
    DATE,
    SITECODE,
    P_mm_day = P_mm_d,
    Q_mm_day = Q_mm_d,
    ET_Zhang_mm_day = ET_PT_zhang,
    ET_Hamon_mm_day = ET_Hamon
  )

write_csv(water_balance_export, file.path(output_dir, "daily_water_balance_zhang_hamon.csv"))

# ---- Export FULL HISTORICAL SERIES (1997-present) with Hamon ----
# This exports all available data back to 1997, not just the complete cases period

# Full Hamon dataset: All variables retained, all years available
hamon_full_export <- results %>%
  filter(DATE >= as.Date("1997-01-01")) %>%  # Start from 1997
  filter(!is.na(ET_Hamon)) %>%  # Only need Hamon ET to be available
  mutate(
    hamon_coefficient_used = ifelse(!is.na(ET_Hamon_calibrated), 
                                    ET_Hamon / ET_Hamon_uncalibrated, 
                                    0.1651)  # Show which coefficient was used
  ) %>%
  select(DATE, SITECODE, everything())

write_csv(hamon_full_export, file.path(output_dir, "daily_ET_watersheds_hamon_full_series.csv"))

# ---- Export COMPREHENSIVE dataset with ALL METHODS for ALL YEARS (1997-present) ----
# This includes ALL ET methods with their original column names, plus supporting variables

comprehensive_export <- results %>%
  filter(DATE >= as.Date("1997-01-01")) %>%
  transmute(
    # Date and location
    DATE,
    SITECODE,
    
    # Meteorological variables
    T_C,                    # Temperature (°C)
    RH_d_pct,              # Relative humidity (%)
    NR_Wm2_d,              # Net radiation (W/m²)
    P_mm_d,                # Precipitation (mm/day)
    Q_mm_d,                # Discharge (mm/day)
    
    # Alpha coefficients
    alpha_zhang = alpha_theoretical_weekly,
    alpha_szilagyi = alpha_szilagyi,
    
    # All ET methods (with original naming scheme)
    et_pt_fixed_1.26 = ET_PT_fixed_1.26,
    et_pt_fixed_0.9 = ET_PT_fixed_0.9, 
    et_pt_zhang = ET_PT_zhang,
    et_pt_szilagyi = ET_PT_szilagyi,
    et_hamon = ET_Hamon
  ) %>%
  arrange(SITECODE, DATE)

write_csv(comprehensive_export, file.path(output_dir, "daily_ET_all_methods_1997_present.csv"))
water_balance_full_export <- results %>%
  filter(DATE >= as.Date("1997-01-01")) %>%
  filter(!is.na(ET_Hamon)) %>%  # Require Hamon ET
  transmute(
    DATE,
    SITECODE,
    P_mm_day = P_mm_d,
    Q_mm_day = Q_mm_d,
    ET_mm_day = ET_Hamon
  )

write_csv(water_balance_full_export, file.path(output_dir, "daily_water_balance_hamon_1997_present.csv"))

# Create summary statistics comparing all methods
method_comparison <- et_long %>%
  group_by(SITECODE, Method) %>%
  summarise(
    mean_ET = mean(ET_mm_day, na.rm = TRUE),
    median_ET = median(ET_mm_day, na.rm = TRUE),
    sd_ET = sd(ET_mm_day, na.rm = TRUE),
    min_ET = min(ET_mm_day, na.rm = TRUE),
    max_ET = max(ET_mm_day, na.rm = TRUE),
    .groups = "drop"
  )

# Export calibration results for future use
if (exists("monthly_calibration") && nrow(monthly_calibration) > 0) {
  write_csv(monthly_calibration, file.path(output_dir, "hamon_monthly_calibration_coefficients.csv"))
}
if (exists("annual_calibration") && nrow(annual_calibration) > 0) {
  write_csv(annual_calibration, file.path(output_dir, "hamon_annual_calibration_coefficients.csv"))
}

# Create calibration validation plot - TO ET_METHODS_COMPARISON FOLDER
if (exists("monthly_calibration") || exists("annual_calibration")) {
  
  # Validation data (2013-2019)
  validation_data <- results %>%
    filter(DATE >= as.Date("2013-01-01") & DATE <= as.Date("2019-12-31")) %>%
    filter(!is.na(ET_PT_zhang) & !is.na(ET_Hamon)) %>%
    select(DATE, SITECODE, ET_PT_zhang, ET_Hamon_uncalibrated, ET_Hamon) %>%
    pivot_longer(
      cols = c(ET_Hamon_uncalibrated, ET_Hamon),
      names_to = "Hamon_Type",
      values_to = "ET_Hamon_value"
    ) %>%
    mutate(Hamon_Type = factor(Hamon_Type,
                               levels = c("ET_Hamon_uncalibrated", "ET_Hamon"),
                               labels = c("Hamon (Uncalibrated)", "Hamon (Calibrated)")))
  
  # Create validation plots for each site
  for (site in unique(validation_data$SITECODE)) {
    
    site_data <- filter(validation_data, SITECODE == site)
    
    p <- ggplot(site_data, aes(x = ET_PT_zhang, y = ET_Hamon_value, color = Hamon_Type)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
      labs(
        title = paste0("Hamon Calibration Validation: ", site),
        subtitle = "2013-2019 Calibration Period",
        x = "PT-Zhang ET (mm/day)",
        y = "Hamon ET (mm/day)",
        color = "Method"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      ) +
      coord_equal(xlim = c(0, max(site_data$ET_PT_zhang, site_data$ET_Hamon_value, na.rm = TRUE)),
                  ylim = c(0, max(site_data$ET_PT_zhang, site_data$ET_Hamon_value, na.rm = TRUE)))
    
    ggsave(file.path(et_plot_dir, paste0("Hamon_calibration_validation_", site, ".png")), 
           p, width = 8, height = 6)
  }
}