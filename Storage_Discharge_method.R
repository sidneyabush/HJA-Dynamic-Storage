# ──────────────────────────────────────────────────────────────────────────────
# Full R script: Dynamic Storage via SD‐method per SITECODE
# (using daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv)
# ──────────────────────────────────────────────────────────────────────────────

# 0) Load libraries & clear environment
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
rm(list = ls())

# 1) Set directories
input_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data/DynamicStorage"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

# 2) Read in the CSV with base R (no stringsAsFactors here)
input_file <- file.path(
  input_dir,
  "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv"
)
if (!file.exists(input_file)) stop("File not found:\n", input_file)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  # Convert DATE (or whatever the header is) into Date
  mutate(DATE = as.Date(DATE)) %>%
  # Rename exactly the columns you have in the file!
  rename(
    date        = DATE,          # your date column
    site        = SITECODE,      # your site code column
    P           = P_mm_d,      # your precipitation column
    Q           = Q_mm_d,      # your discharge column
    ET          = ET_mm_d      # your ET column
  ) %>%
  arrange(site, date)

# Quick sanity check: you should see date, site, P, Q, ET
str(df)
head(df)


# 3) Helper functions

# 3.1 Flow‐Duration Curve (Gringorten plotting positions)
compute_fdc <- function(Q){
  Qpos   <- Q[Q > 0]
  Qsort  <- sort(Qpos, decreasing = TRUE)
  n      <- length(Qsort)
  r      <- seq_len(n)
  exceed <- (r - 0.44) / (n + 0.12) * 100
  tibble(exceedance = exceed, Q = Qsort)
}

# 3.2 Safe interpolator for FDC thresholds
get_Q <- function(fdc, P_exc){
  # if fdc has fewer than 2 points, return NA
  if (nrow(fdc) < 2) return(NA_real_)
  approx(fdc$exceedance, fdc$Q, xout = P_exc)$y
}

# 3.3 Analytical ΔS from Kirchner (2009) eqn (3)
deltaS <- function(Q_up, Q_lo, k, p){
  (Q_up^(2 - p) - Q_lo^(2 - p)) / (k * (2 - p))
}


# 4) Site‐level analysis
analyze_site <- function(df_site){
  # A) Build and check FDC
  fdc <- compute_fdc(df_site$Q)
  Q99 <- get_Q(fdc, 99)
  Q50 <- get_Q(fdc, 50)
  Q01 <- get_Q(fdc, 1)
  
  # B) If FDC is invalid, return NA row
  if (is.na(Q99) || is.na(Q50) || is.na(Q01)) {
    return(tibble(
      site           = unique(df_site$site),
      k              = NA_real_,
      p              = NA_real_,
      Q_max          = NA_real_,
      Q_min          = NA_real_,
      Q99            = Q99,
      Q50            = Q50,
      Q01            = Q01,
      S_annual       = NA_real_,
      S_low_to_med   = NA_real_,
      S_med_to_high  = NA_real_
    ))
  }
  
  # C) Prepare recession data (rain‐free falling limbs)
  df_site <- df_site %>%
    mutate(
      dt      = as.numeric(difftime(date, lag(date), units = "days")),
      dQ      = (lag(Q) - Q) / dt,
      is_rain = P > 0
    )
  
  rec <- df_site %>%
    filter(
      !is_rain,
      !is.na(dQ),
      dQ > 0,
      Q  > 0
    )
  
  # D) If too few recession points, return partial results
  if (nrow(rec) < 10) {
    return(tibble(
      site           = unique(df_site$site),
      k              = NA_real_,
      p              = NA_real_,
      Q_max          = max(df_site$Q, na.rm = TRUE),
      Q_min          = min(df_site$Q[df_site$Q > 0], na.rm = TRUE),
      Q99            = Q99,
      Q50            = Q50,
      Q01            = Q01,
      S_annual       = NA_real_,
      S_low_to_med   = NA_real_,
      S_med_to_high  = NA_real_
    ))
  }
  
  # E) Fit the power‐law recession: -dQ/dt = k Q^p
  fit   <- lm(log(dQ) ~ log(Q), data = rec)
  p_est <- coef(fit)["log(Q)"]
  k_est <- exp(coef(fit)["(Intercept)"])
  
  # F) Compute dynamic storage
  Q_max     <- max(df_site$Q, na.rm = TRUE)
  Q_min_pos <- min(df_site$Q[df_site$Q > 0], na.rm = TRUE)
  
  S_annual     <- deltaS(Q_max,      Q_min_pos, k_est, p_est)
  S_low_to_med <- deltaS(Q50,         Q99,       k_est, p_est)
  S_med_to_high<- deltaS(Q01,         Q50,       k_est, p_est)
  
  tibble(
    site           = unique(df_site$site),
    k              = k_est,
    p              = p_est,
    Q_max          = Q_max,
    Q_min          = Q_min_pos,
    Q99            = Q99,
    Q50            = Q50,
    Q01            = Q01,
    S_annual       = S_annual,
    S_low_to_med   = S_low_to_med,
    S_med_to_high  = S_med_to_high
  )
}


# 5) Run the analysis for every site
results <- df %>%
  group_by(site) %>%
  nest() %>%
  mutate(analysis = map(data, analyze_site)) %>%
  select(-data) %>%
  unnest(analysis)

# 6) Write out
write.csv(
  results,
  file.path(output_dir, "dynamic_storage_results_per_site.csv"),
  row.names = FALSE
)
# ──────────────────────────────────────────────────────────────────────────────
