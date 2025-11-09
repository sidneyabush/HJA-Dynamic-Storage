# =============================================================================
# Dynamic Storage: Overall + Annual per Site 
# =============================================================================

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

# 2) Read data, drop unwanted sites, parse date & compute water-year
library(lubridate)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    # robustly parse DATE into a Date object
    date = parse_date_time(
      DATE,
      orders = c("Ymd", "Y-m-d", "mdy", "dmy")
    ) %>% as.Date(),
    
    # water-year: Oct–Dec → next calendar year, Jan–Sep → same year
    wateryear = if_else(
      month(date) >= 10,
      year(date) + 1,
      year(date)
    )
  ) %>%
  rename(
    site = SITECODE,
    P    = P_mm_d,
    Q    = Q_mm_d,
    ET   = ET_mm_d
  ) %>%
  filter(!site %in% c("COLD", "LONGER")) %>%
  arrange(site, date)


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

# 4) Site-water-year analysis function
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
      S_annual_mm   = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm  = NA_real_
    ))
  }
  
  # 4.2 Recession data (rain-free falling limbs)
  tmp <- df_sub %>%
    arrange(date) %>%
    mutate(
      dt      = as.numeric(difftime(date, lag(date), "days")),
      dQ      = (lag(Q) - Q) / dt,
      is_rain = P > 0
    )
  rec <- tmp %>% filter(!is_rain, !is.na(dQ), dQ > 0, Q > 0)
  if (nrow(rec) < 10) {
    return(tibble(
      k = NA_real_, p = NA_real_,
      Q_max = NA_real_, Q_min = NA_real_,
      Q99 = Q99, Q50 = Q50, Q01 = Q01,
      S_annual_mm   = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm  = NA_real_
    ))
  }
  
  # 4.3 Fit recession law
  fit   <- lm(log(dQ) ~ log(Q), data = rec)
  p_est <- coef(fit)["log(Q)"]
  k_est <- exp(coef(fit)["(Intercept)"])
  
  # 4.4 Compute annual-mean extremes by water-year
  yearly <- df_sub %>%
    group_by(wateryear) %>%
    summarize(
      yr_max = if (all(is.na(Q))) NA_real_ else max(Q, na.rm = TRUE),
      yr_min = if (all(Q <= 0 | is.na(Q))) NA_real_ else min(Q[Q > 0], na.rm = TRUE),
      .groups = "drop"
    )
  Qmax_avg <- mean(yearly$yr_max, na.rm = TRUE)
  Qmin_avg <- mean(yearly$yr_min, na.rm = TRUE)
  
  # 4.5 Compute dynamic storage depths
  tibble(
    k             = k_est,
    p             = p_est,
    Q_max         = Qmax_avg,
    Q_min         = Qmin_avg,
    Q99           = Q99,
    Q50           = Q50,
    Q01           = Q01,
    S_annual_mm   = deltaS(Qmax_avg, Qmin_avg, k_est, p_est),
    S_high_med_mm = deltaS(Q01, Q50, k_est, p_est),
    S_med_low_mm  = deltaS(Q50, Q99, k_est, p_est)
  )
}

# 5) Compute overall results (across all water-years)
overall <- df %>%
  group_by(site) %>%
  group_map(~ analyze(.x) %>% mutate(site = .y$site), .keep = TRUE) %>%
  bind_rows()

# 6) Compute per-water-year results
annual <- df %>%
  group_by(site, wateryear) %>%
  group_map(~ analyze(.x) %>% 
              mutate(site = .y$site, wateryear = .y$wateryear),
            .keep = TRUE) %>%
  bind_rows()

# 7) Add catchment areas & convert mm → m³
areas_ha <- tibble(
  SITECODE = c("GSWS01","GSWS02","GSWS03","GSWS06","GSWS07",
               "GSWS08","GSWS09","GSWS10", "GSMACK", "GSLOOK_FULL"),
  area_ha  = c(96, 60, 101, 13.0, 15.4, 21.4, 8.5, 10.2, 580, 6242)
)
add_vol <- function(df_in) {
  df_in %>%
    left_join(areas_ha, by = c("site" = "SITECODE")) %>%
    mutate(
      area_m2        = area_ha * 10000,
      S_annual_m3    = (S_annual_mm   / 1000) * area_m2,
      S_high_med_m3  = (S_high_med_mm / 1000) * area_m2,
      S_med_low_m3   = (S_med_low_mm  / 1000) * area_m2
    )
}

overall_vol <- add_vol(overall)
annual_vol  <- add_vol(annual)

# 8) Save outputs (with wateryear in filenames for clarity)
write.csv(overall_vol,
          file = file.path(output_dir, "storage_overall_per_site.csv"),
          row.names = FALSE)
write.csv(annual,
          file = file.path(output_dir,
                           "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"),
          row.names = FALSE)
write.csv(annual_vol,
          file = file.path(output_dir,
                           "storage_discharge_method_annual_vol_per_site_wateryear.csv"),
          row.names = FALSE)

# 9) Quick plot of annual dynamic storage depths vs. water-year
library(ggplot2)
ggplot(annual_vol, aes(x = wateryear, y = S_annual_mm, color = site, group = site)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Annual Dynamic Storage Depth by Water-Year",
    x     = "Water Year",
    y     = "Storage Depth (mm)",
    color = "Site"
  ) +
  theme_minimal()

# reshape to long format
annual_long <- annual_vol %>%
  select(site, wateryear, S_annual_mm, S_high_med_mm, S_med_low_mm) %>%
  pivot_longer(
    cols = starts_with("S_"),
    names_to  = "metric",
    values_to = "storage_mm"
  )

# facet by metric
ggplot(annual_long, aes(x = wateryear, y = storage_mm, color = site)) +
  geom_line() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Water Year", y = "Storage (mm)", color = "Site") +
  theme_minimal()

ggplot(annual_long, aes(x = wateryear, y = storage_mm, 
                        linetype = metric, color = metric)) +
  geom_line() +
  facet_wrap(~ site) +
  labs(x = "Water Year", y = "Storage (mm)") +
  theme_classic()

