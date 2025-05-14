library(readr)
library(dplyr)
library(lubridate)

# ── Clear environment ─────────────────────────────────────────────────────────
rm(list = ls())

# ── Directories & site list ─────────────────────────────────────────────────
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
sites_keep <- c(
  "GSLOOK","GSWS01","GSWS02","GSWS03","GSWS09",
  "GSWS10","GSWSMC","GSWS06","GSWS07","GSWS08"
)

# ── 1) Read & parse 5-min data ────────────────────────────────────────────────
Q_5min <- read_csv(
  file.path(base_dir, "Q", "HF00401_v17.txt"),
  show_col_types = FALSE
) %>%
  filter(WATERYEAR > 1997, SITECODE %in% sites_keep) %>%
  mutate(
    DATE_TIME  = as.POSIXct(DATE_TIME,
                            format = "%Y-%m-%d %H:%M:%S",
                            tz     = "UTC"),
    DATE       = as.Date(DATE_TIME),
    HOUR       = floor_date(DATE_TIME, "hour"),
    Q          = INST_Q * 0.02831683199881,    # convert to m³/s
    WATERYEAR  = if_else(month(DATE) >= 10,
                         year(DATE) + 1,
                         year(DATE))
  )

EC_5min <- read_csv(
  file.path(base_dir, "EC", "CF01201_v3.txt"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE_TIME  = as.POSIXct(DATE_TIME,
                            format = "%Y-%m-%d %H:%M:%S",
                            tz     = "UTC"),
    DATE       = as.Date(DATE_TIME),
    HOUR       = floor_date(DATE_TIME, "hour"),
    WATERYEAR  = if_else(month(DATE) >= 10,
                         year(DATE) + 1,
                         year(DATE))) %>%
  filter(WATERYEAR > 1997, SITECODE %in% sites_keep)
  

# Hourly summaries ──────────────────────────────────────────────────────
hourly_Q <- Q_5min %>%
  group_by(SITECODE, HOUR, WATERYEAR) %>%
  summarize(
    MEAN_Q        = mean(Q,   na.rm = TRUE),
    SD_Q          = sd(Q,     na.rm = TRUE),
    N_READINGS_Q  = n(),
    .groups       = "drop"
  )

hourly_EC_T <- EC_5min %>%
  group_by(SITECODE, HOUR, WATERYEAR) %>%
  summarize(
    MEAN_EC       = mean(EC_INST,        na.rm = TRUE),
    SD_EC         = sd(EC_INST,          na.rm = TRUE),
    MEAN_TEMP     = mean(WATERTEMP_INST, na.rm = TRUE),
    SD_TEMP       = sd(WATERTEMP_INST,   na.rm = TRUE),
    N_READINGS_T  = n(),
    .groups       = "drop"
  )

hourly_data <- hourly_Q %>%
  left_join(hourly_EC_T, by = c("SITECODE", "HOUR", "WATERYEAR"))

# Daily summaries ───────────────────────────────────────────────────────
daily_Q <- Q_5min %>%
  group_by(SITECODE, DATE, WATERYEAR) %>%
  summarize(
    MEAN_Q        = mean(Q,   na.rm = TRUE),
    SD_Q          = sd(Q,     na.rm = TRUE),
    N_READINGS_Q  = n(),
    .groups       = "drop"
  )

daily_EC_T <- EC_5min %>%
  group_by(SITECODE, DATE, WATERYEAR) %>%
  summarize(
    MEAN_EC       = mean(EC_INST,        na.rm = TRUE),
    SD_EC         = sd(EC_INST,          na.rm = TRUE),
    MEAN_TEMP     = mean(WATERTEMP_INST, na.rm = TRUE),
    SD_TEMP       = sd(WATERTEMP_INST,   na.rm = TRUE),
    N_READINGS_T  = n(),
    .groups       = "drop"
  )

daily_data <- daily_Q %>%
  left_join(daily_EC_T, by = c("SITECODE", "DATE", "WATERYEAR"))


