# 02_fill_Hamon_methods.R

library(readr)
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)

# --- Directories ---
input_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"
output_dir <- input_dir

# --- 1. Load PT‐ET timeseries & compute uncalibrated Hamon ---
pt <- read_csv(file.path(input_dir, "PT_ET_methods_timeseries.csv"),
               show_col_types = FALSE)

calculate_et_hamon <- function(temp_c, date_vec, lat, coeff) {
  doy      <- yday(date_vec)
  δ        <- 0.4093 * sin((2*pi*doy/365) - 1.405)
  φ_rad    <- lat * pi/180
  cos_sun  <- pmax(-1, pmin(1, -tan(φ_rad)*tan(δ)))
  ω        <- acos(cos_sun)
  dayl     <- (24/pi)*ω
  es       <- 0.6108 * exp(17.27*temp_c/(temp_c+237.3))
  ρ_v      <- 216.7 * es/(temp_c+273.16)
  et_raw   <- coeff * (dayl/12) * ρ_v
  ifelse(is.na(temp_c), NA_real_, pmax(0, et_raw))
}

study_lat <- (44.20127400 + 44.28226000)/2

data <- pt %>%
  mutate(
    ET_Hamon_uncalibrated = calculate_et_hamon(T_C, DATE, study_lat, 0.1651),
    MONTH                 = month(DATE)
  )

# --- 2. Derive monthly median coefficients (2013–2019) ---
cal_window_raw <- data %>%
  filter(DATE >= as.Date("2013-01-01"),
         DATE <= as.Date("2019-12-31"))

monthly_coef_zhang <- cal_window_raw %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_uncalibrated)) %>%
  group_by(SITECODE, MONTH) %>%
  summarise(coef_zhang_month = median(ET_PT_zhang / ET_Hamon_uncalibrated),
            .groups="drop")

monthly_coef_szilagyi <- cal_window_raw %>%
  filter(!is.na(ET_PT_szilagyi), !is.na(ET_Hamon_uncalibrated)) %>%
  group_by(SITECODE, MONTH) %>%
  summarise(coef_szilagyi_month = median(ET_PT_szilagyi / ET_Hamon_uncalibrated),
            .groups="drop")

# --- 3. Fallback: site‐mean & global median coefficients ---
site_mean_zhang    <- monthly_coef_zhang    %>% group_by(SITECODE) %>% summarise(coef_zhang_site    = mean(coef_zhang_month),    .groups="drop")
site_mean_szilagyi <- monthly_coef_szilagyi %>% group_by(SITECODE) %>% summarise(coef_szilagyi_site = mean(coef_szilagyi_month), .groups="drop")

global_coef_zhang    <- mean(monthly_coef_zhang$coef_zhang_month,    na.rm=TRUE)
global_coef_szilagyi <- mean(monthly_coef_szilagyi$coef_szilagyi_month, na.rm=TRUE)

# --- 4. Merge in coefficients & compute monthly‐calibrated Hamon ---
data <- data %>%
  left_join(monthly_coef_zhang,    by = c("SITECODE","MONTH")) %>%
  left_join(site_mean_zhang,       by = "SITECODE") %>%
  left_join(monthly_coef_szilagyi, by = c("SITECODE","MONTH")) %>%
  left_join(site_mean_szilagyi,    by = "SITECODE") %>%
  mutate(
    coef_zhang_final    = coalesce(coef_zhang_month,    coef_zhang_site,    global_coef_zhang),
    coef_szilagyi_final = coalesce(coef_szilagyi_month, coef_szilagyi_site, global_coef_szilagyi),
    ET_Hamon_pt_zhang_monthly   = ifelse(is.na(ET_Hamon_uncalibrated),
                                         NA_real_,
                                         pmax(0, ET_Hamon_uncalibrated * coef_zhang_final)),
    ET_Hamon_pt_szilagyi_monthly = ifelse(is.na(ET_Hamon_uncalibrated),
                                          NA_real_,
                                          pmax(0, ET_Hamon_uncalibrated * coef_szilagyi_final))
  )

# --- 5. Re‐define calibration window on the _new_ columns ---
cal_window <- data %>%
  filter(DATE >= as.Date("2013-01-01"),
         DATE <= as.Date("2019-12-31"))

# --- 6. Fit per‐site regressions (require ≥10 pts) ---
reg_sites_zhang    <- cal_window %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_pt_zhang_monthly)) %>%
  pull(SITECODE) %>% unique()

reg_sites_szilagyi <- cal_window %>%
  filter(!is.na(ET_PT_szilagyi), !is.na(ET_Hamon_pt_szilagyi_monthly)) %>%
  pull(SITECODE) %>% unique()

regressions_zhang <- map(reg_sites_zhang, function(site) {
  df <- cal_window %>% filter(SITECODE == site)
  if(nrow(df) >= 10) lm(ET_Hamon_pt_zhang_monthly ~ ET_PT_zhang, data = df) else NULL
}) %>% set_names(reg_sites_zhang)

regressions_szilagyi <- map(reg_sites_szilagyi, function(site) {
  df <- cal_window %>% filter(SITECODE == site)
  if(nrow(df) >= 10) lm(ET_Hamon_pt_szilagyi_monthly ~ ET_PT_szilagyi, data = df) else NULL
}) %>% set_names(reg_sites_szilagyi)

# --- 7. Predict via regression where monthly‐calibrated is NA but PT exists ---
pred_zhang <- map_dfr(reg_sites_zhang, function(site) {
  m <- regressions_zhang[[site]]
  df <- data %>%
    filter(SITECODE == site,
           is.na(ET_Hamon_pt_zhang_monthly),
           !is.na(ET_PT_zhang))
  if(!is.null(m) && nrow(df)>0) {
    df %>% transmute(SITECODE, DATE,
                     ET_Hamon_pt_zhang_pred = pmax(0, predict(m, newdata = df)))
  } else {
    tibble(SITECODE=character(), DATE=as.Date(character()), ET_Hamon_pt_zhang_pred=double())
  }
})

pred_szilagyi <- map_dfr(reg_sites_szilagyi, function(site) {
  m <- regressions_szilagyi[[site]]
  df <- data %>%
    filter(SITECODE == site,
           is.na(ET_Hamon_pt_szilagyi_monthly),
           !is.na(ET_PT_szilagyi))
  if(!is.null(m) && nrow(df)>0) {
    df %>% transmute(SITECODE, DATE,
                     ET_Hamon_pt_szilagyi_pred = pmax(0, predict(m, newdata = df)))
  } else {
    tibble(SITECODE=character(), DATE=as.Date(character()), ET_Hamon_pt_szilagyi_pred=double())
  }
})

data <- data %>%
  left_join(pred_zhang,    by = c("SITECODE","DATE")) %>%
  left_join(pred_szilagyi, by = c("SITECODE","DATE"))

# --- 8. Create initial interpolated‐by‐coefficient‐and‐regression series ---
data <- data %>%
  mutate(
    ET_Hamon_pt_zhang_interpolated   = coalesce(ET_Hamon_pt_zhang_monthly,
                                                ET_Hamon_pt_zhang_pred),
    ET_Hamon_pt_szilagyi_interpolated = coalesce(ET_Hamon_pt_szilagyi_monthly,
                                                 ET_Hamon_pt_szilagyi_pred)
  )

# load zoo for the time interpolation helpers
library(zoo)

# --- 9. Fill any remaining NA by time‐series interpolation within each site ---
data <- data %>%
  arrange(SITECODE, DATE) %>%
  group_by(SITECODE) %>%
  mutate(
    # 9a) linear interpolate interior gaps
    ET_Hamon_pt_zhang_interp_full   = na.approx(ET_Hamon_pt_zhang_interpolated,
                                                x = DATE, na.rm = FALSE),
    ET_Hamon_pt_szilagyi_interp_full = na.approx(ET_Hamon_pt_szilagyi_interpolated,
                                                 x = DATE, na.rm = FALSE),
    
    # 9b) carry‐forward to fill any remaining leading NAs
    ET_Hamon_pt_zhang_interp_full   = na.locf(ET_Hamon_pt_zhang_interp_full,
                                              na.rm = FALSE),
    ET_Hamon_pt_szilagyi_interp_full = na.locf(ET_Hamon_pt_szilagyi_interp_full,
                                               na.rm = FALSE),
    
    # 9c) carry‐backward to fill any remaining trailing NAs
    ET_Hamon_pt_zhang_interp_full   = na.locf(ET_Hamon_pt_zhang_interp_full,
                                              fromLast = TRUE),
    ET_Hamon_pt_szilagyi_interp_full = na.locf(ET_Hamon_pt_szilagyi_interp_full,
                                               fromLast = TRUE)
  ) %>%
  ungroup()

# --- 10. Drop all the temporary coefficient & prediction columns ---
data <- data %>%
  select(-contains("coef"),
         -contains("alpha"),
         -contains("interpolated"),
         -contains("pred"),
         -MONTH
         )

# --- 11. Export the fully‐filled daily water‐balance record (no Hamon gaps) ---
write_csv(
  data,
  file.path(output_dir, "daily_MET_water_balance_all_ET_methods_1997_present.csv")
)

# --- 12. Drop all the temporary coefficient & prediction columns ---
WB_data <- data %>%
  select(-RH_d_pct, -NR_Wm2_d,-VPD_kPa
  )

# --- 13. Export the fully‐filled daily water‐balance record (no Hamon gaps) ---
write_csv(
  WB_data,
  file.path(output_dir, "daily_water_balance_all_ET_methods_1997_present.csv")
)

# --- 14. Drop all the temporary coefficient & prediction columns ---
WB_data_final_ET <- data %>%
  dplyr::select(-RH_d_pct, -NR_Wm2_d,-VPD_kPa, 
         -contains("szilagyi"),
         -contains("monthly"),
         -contains("uncalibrated"), 
         -ET_PT_zhang) %>%
  dplyr::rename(ET_mm_d = ET_Hamon_pt_zhang_interp_full)

# --- 15. Export the fully‐filled daily water‐balance record (no Hamon gaps) ---
write_csv(
  WB_data_final_ET,
  file.path(output_dir, "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv")
)
