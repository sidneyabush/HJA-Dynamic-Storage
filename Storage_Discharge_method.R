# ──────────────────────────────────────────────────────────────────────────────
# Dynamic Storage: Overall + Annual per Site
# (excluding COLD, LONGER, MACK; outputs in mm and m³)
# ──────────────────────────────────────────────────────────────────────────────

# 0) Load libraries & clear environment
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
rm(list = ls())

# 1) Directories & input
input_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data/DynamicStorage"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"
input_file <- file.path(input_dir, "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv")
stopifnot(file.exists(input_file))

# 2) Read data and drop unwanted sites
#    Rename columns and add calendar year
df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(DATE = as.Date(DATE)) %>%
  rename(
    date = DATE,
    site = SITECODE,
    P    = P_mm_d,    # precipitation (mm/day)
    Q    = Q_mm_d,    # discharge    (mm/day)
    ET   = ET_mm_d    # evapotransp. (mm/day)
  ) %>%
  filter(!site %in% c("COLD", "LONGER", "GSMACK")) %>%
  arrange(site, date) %>%
  mutate(year = year(date))

# 3) Helper functions
#  3.1 Flow-Duration Curve (Gringorten)
compute_fdc <- function(Q) {
  Qpos  <- Q[Q > 0]
  Qs    <- sort(Qpos, decreasing = TRUE)
  n     <- length(Qs)
  r     <- seq_len(n)
  exc   <- (r - 0.44) / (n + 0.12) * 100
  tibble(exceedance = exc, Q = Qs)
}

#  3.2 Safe extraction of Q thresholds
get_Q <- function(fdc, P_exc) {
  if (nrow(fdc) < 2) return(NA_real_)
  approx(fdc$exceedance, fdc$Q, xout = P_exc)$y
}

#  3.3 Analytical dynamic storage (Kirchner 2009)
deltaS <- function(Qu, Ql, k, p) {
  (Qu^(2 - p) - Ql^(2 - p)) / (k * (2 - p))
}

#  3.4 Site-year analysis function
analyze <- function(df_sub) {
  fdc  <- compute_fdc(df_sub$Q)
  Q99  <- get_Q(fdc, 99); Q50 <- get_Q(fdc, 50); Q01 <- get_Q(fdc, 1)
  # if thresholds invalid, return NAs
  if (any(is.na(c(Q99, Q50, Q01)))) {
    return(tibble(
      k             = NA_real_,
      p             = NA_real_,
      Q_max         = NA_real_,
      Q_min         = NA_real_,
      Q99           = Q99,
      Q50           = Q50,
      Q01           = Q01,
      S_annual_mm   = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm  = NA_real_
    ))
  }
  tmp <- df_sub %>%
    arrange(date) %>%
    mutate(
      dt      = as.numeric(difftime(date, lag(date), "days")),
      dQ      = (lag(Q) - Q) / dt,
      is_rain = P > 0
    )
  rec <- tmp %>% filter(!is_rain, !is.na(dQ), dQ > 0, Q > 0)
  # if too few recession points
  if (nrow(rec) < 10) {
    return(tibble(
      k             = NA_real_,
      p             = NA_real_,
      Q_max         = max(tmp$Q, na.rm = TRUE),
      Q_min         = min(tmp$Q[tmp$Q > 0], na.rm = TRUE),
      Q99           = Q99,
      Q50           = Q50,
      Q01           = Q01,
      S_annual_mm   = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm  = NA_real_
    ))
  }
  fit   <- lm(log(dQ) ~ log(Q), data = rec)
  p_est <- coef(fit)["log(Q)"]
  k_est <- exp(coef(fit)["(Intercept)"])
  Qmax  <- max(df_sub$Q, na.rm = TRUE)
  Qmin  <- min(df_sub$Q[df_sub$Q > 0], na.rm = TRUE)
  # compute depth in mm
  S_ann_mm   <- deltaS(Qmax, Qmin, k_est, p_est)
  S_hm_mm    <- deltaS(Q01, Q50, k_est, p_est)
  S_ml_mm    <- deltaS(Q50, Q99, k_est, p_est)
  tibble(
    k             = k_est,
    p             = p_est,
    Q_max         = Qmax,
    Q_min         = Qmin,
    Q99           = Q99,
    Q50           = Q50,
    Q01           = Q01,
    S_annual_mm   = S_ann_mm,
    S_high_med_mm = S_hm_mm,
    S_med_low_mm  = S_ml_mm
  )
}

# 4) Compute overall and annual results
overall <- df %>%
  group_by(site) %>%
  nest() %>%
  mutate(res = map(data, analyze)) %>%
  select(site, res) %>%
  unnest(res)

annual <- df %>%
  group_by(site, year) %>%
  nest() %>%
  mutate(res = map(data, analyze)) %>%
  select(site, year, res) %>%
  unnest(res)

# 5) Add areas and convert depths (mm) to volumes (m³)
areas_ha <- tibble(
  SITECODE = c("GSWS01","GSWS02","GSWS03","GSWS06","GSWS07",
               "GSWS08","GSWS09","GSWS10","GSLOOK_FULL"),
  area_ha   = c(96, 60, 101, 13.0, 15.4, 21.4,  8.5, 10.2, 6242)
)
add_vol <- function(df_in) {
  df_in %>%
    left_join(areas_ha, by = c("site" = "SITECODE")) %>%
    mutate(
      area_m2           = area_ha * 10000,
      S_annual_m3       = (S_annual_mm      / 1000) * area_m2,
      S_high_med_m3     = (S_high_med_mm    / 1000) * area_m2,
      S_med_low_m3      = (S_med_low_mm     / 1000) * area_m2
    )
}

overall_vol <- add_vol(overall)
annual_vol  <- add_vol(annual)

# 6) Save outputs
write.csv(
  overall_vol,
  file.path(output_dir, "storage_discharge_overall_per_site.csv"),
  row.names = FALSE
)
write.csv(
  annual_vol,
  file.path(output_dir, "storage_discharge_annual_per_site_year.csv"),
  row.names = FALSE
)
