# ──────────────────────────────────────────────────────────────────────────────
# Dynamic Storage: Overall + Annual per Site
# (excluding COLD, LONGER, GSWSMC, GSMACK; outputs in mm and m³)
# ──────────────────────────────────────────────────────────────────────────────

# 0) Load libraries & clear environment
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
rm(list = ls())

# 1) Directories & input
data_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data/DynamicStorage"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"
input_file <- file.path(data_dir, "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv")
stopifnot(file.exists(input_file))

# 2) Read data and drop unwanted sites
#    Rename columns and add calendar year
df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(DATE = as.Date(DATE)) %>%
  rename(
    date = DATE,
    site = SITECODE,
    P    = P_mm_d,
    Q    = Q_mm_d,
    ET   = ET_mm_d
  ) %>%
  filter(!site %in% c("COLD", "LONGER", "GSWSMC", "GSMACK")) %>%
  arrange(site, date) %>%
  mutate(year = year(date))

# 3) Helper functions
compute_fdc <- function(Q) {
  Qpos <- Q[Q > 0]
  Qs   <- sort(Qpos, decreasing = TRUE)
  n    <- length(Qs)
  r    <- seq_len(n)
  exc  <- (r - 0.44) / (n + 0.12) * 100
  tibble(exceedance = exc, Q = Qs)
}

get_Q <- function(fdc, P_exc) {
  if (nrow(fdc) < 2) return(NA_real_)
  approx(fdc$exceedance, fdc$Q, xout = P_exc)$y
}

deltaS <- function(Qu, Ql, k, p) {
  (Qu^(2 - p) - Ql^(2 - p)) / (k * (2 - p))
}

# 4) Site-year analysis function with annual-mean extremes
analyze <- function(df_sub) {
  # 4.1 Flow-duration thresholds
  fdc  <- compute_fdc(df_sub$Q)
  Q99  <- get_Q(fdc, 99)
  Q50  <- get_Q(fdc, 50)
  Q01  <- get_Q(fdc, 1)
  if (any(is.na(c(Q99, Q50, Q01)))) {
    return(tibble(
      k = NA_real_, p = NA_real_,
      Q_max = NA_real_, Q_min = NA_real_,
      Q99 = Q99, Q50 = Q50, Q01 = Q01,
      S_annual_mm = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm = NA_real_
    ))
  }
  
  # 4.2 Recession data (rain-free falling limbs)
  tmp <- df_sub %>%
    arrange(date) %>%
    mutate(
      dt = as.numeric(difftime(date, lag(date), "days")),
      dQ = (lag(Q) - Q) / dt,
      is_rain = P > 0
    )
  rec <- tmp %>% filter(!is_rain, !is.na(dQ), dQ > 0, Q > 0)
  if (nrow(rec) < 10) {
    return(tibble(
      k = NA_real_, p = NA_real_,
      Q_max = NA_real_, Q_min = NA_real_,
      Q99 = Q99, Q50 = Q50, Q01 = Q01,
      S_annual_mm = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm = NA_real_
    ))
  }
  
  # 4.3 Fit recession law
  fit   <- lm(log(dQ) ~ log(Q), data = rec)
  p_est <- coef(fit)["log(Q)"]
  k_est <- exp(coef(fit)["(Intercept)"])
  
  # 4.4 Compute annual-mean extremes
  yearly <- df_sub %>%
    group_by(year) %>%
    summarize(
      yr_max = if (all(is.na(Q))) NA_real_ else max(Q, na.rm = TRUE),
      yr_min = if (all(Q <= 0 | is.na(Q))) NA_real_ else min(Q[Q > 0], na.rm = TRUE)
    )
  Qmax_avg <- mean(yearly$yr_max, na.rm = TRUE)
  Qmin_avg <- mean(yearly$yr_min, na.rm = TRUE)
  
  # 4.5 Compute dynamic storage depths in mm
  S_ann_mm   <- deltaS(Qmax_avg, Qmin_avg, k_est, p_est)
  S_hm_mm    <- deltaS(Q01, Q50, k_est, p_est)
  S_ml_mm    <- deltaS(Q50, Q99, k_est, p_est)
  
  tibble(
    k             = k_est,
    p             = p_est,
    Q_max         = Qmax_avg,
    Q_min         = Qmin_avg,
    Q99           = Q99,
    Q50           = Q50,
    Q01           = Q01,
    S_annual_mm   = S_ann_mm,
    S_high_med_mm = S_hm_mm,
    S_med_low_mm  = S_ml_mm
  )
}

# 5) Compute overall results
overall <- bind_rows(
  df %>%
    group_by(site) %>%
    group_map(~ analyze(.x) %>% mutate(site = .y$site), .keep = TRUE)
)

# 6) Compute annual results
annual <- bind_rows(
  df %>%
    group_by(site, year) %>%
    group_map(~ analyze(.x) %>% mutate(site = .y$site, year = .y$year), .keep = TRUE)
)

# 7) Add catchment areas and convert depths to volumes
areas_ha <- tibble(
  SITECODE = c("GSWS01","GSWS02","GSWS03","GSWS06","GSWS07",
               "GSWS08","GSWS09","GSWS10","GSLOOK_FULL"),
  area_ha = c(96, 60, 101, 13.0, 15.4, 21.4, 8.5, 10.2, 6242)
)
add_vol <- function(df_in) {
  df_in %>%
    left_join(areas_ha, by = c("site" = "SITECODE")) %>%
    mutate(
      area_m2 = area_ha * 10000,
      S_annual_m3   = (S_annual_mm   / 1000) * area_m2,
      S_high_med_m3 = (S_high_med_mm / 1000) * area_m2,
      S_med_low_m3  = (S_med_low_mm  / 1000) * area_m2
    )
}

overall_vol <- add_vol(overall)
annual_vol  <- add_vol(annual)

# 8) Save outputs
write.csv(overall_vol, file.path(output_dir, "storage_overall_per_site.csv"), row.names = FALSE)
write.csv(annual_vol,  file.path(output_dir, "storage_annual_per_site_year.csv"), row.names = FALSE)


# 9) Quick plot of annual dynamic storage depths
#    Plots S_annual_mm over years for each site
library(ggplot2)
ggplot(annual_vol, aes(x = year, y = S_annual_mm, color = site, group = site)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Annual Dynamic Storage (Depth) per Site",
    x = "Year",
    y = "Dynamic Storage (mm)",
    color = "Site"
  ) +
  theme_minimal()
