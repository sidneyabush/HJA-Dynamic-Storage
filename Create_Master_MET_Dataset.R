library(readr)
library(dplyr)
library(lubridate)
library(tidyr)

rm(list = ls())

met_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data/all_met"

parse_my_date <- function(d) {
  parse_date_time(d,
                  orders = c("Ymd", "Y-m-d", "mdy", "m/d/Y", "dmy", "d/m/Y")
  ) %>% as_date()
}

make_inter_long <- function(fname, var) {
  read_csv(file.path(met_dir, fname), show_col_types = FALSE) %>%
    mutate(DATE = parse_my_date(DATE)) %>%
    pivot_longer(
      cols         = ends_with("_inter"),
      names_to     = "SITECODE",
      names_pattern= "(.*)_inter$",     # ← capture everything before “_inter”
      values_to    = var                # ← e.g. "Temp" or "Precip"
    )
}

Temp   <- make_inter_long("Temperature_original_&_filled_1979_2023_v2.csv", "Temp")
Temp <- Temp %>%
  dplyr::select(DATE, SITECODE, Temp) %>%
  dplyr::rename(T_C = Temp)

Precip <- make_inter_long("Precipitation_original_&_filled_1979_2023.csv",  "Precip")
Precip <- Precip %>%
  dplyr::select(DATE, SITECODE, Precip) %>%
  dplyr::rename(P_mm_d = Precip)

# For the other two, same trick:
RH <- read_csv(file.path(met_dir, "MS00102_v9.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  dplyr::select(SITECODE, DATE, RELHUM_MEAN_DAY) %>%
  dplyr::rename(RH_d_pct = RELHUM_MEAN_DAY)

NetRad <- read_csv(file.path(met_dir, "MS05025_v3.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  select(SITECODE, DATE, NR_TOT_MEAN_DAY) %>%
  dplyr::rename(NR_Wm2_d = NR_TOT_MEAN_DAY)

combined_met <- Temp %>%
  full_join(Precip, by = c("DATE", "SITECODE")) %>%
  full_join(RH,     by = c("DATE", "SITECODE")) %>%
  full_join(NetRad, by = c("DATE", "SITECODE")) %>%
  filter(DATE >= ymd("2013-04-17")) %>%
  arrange(SITECODE, DATE)

# 1) pivot your master into long form
long_met <- combined_met %>%
  pivot_longer(
    cols      = c(P_mm_d, T_C, RH_d_pct, NR_Wm2_d),
    names_to  = "var",
    values_to = "value"
  )

combined_unique <- combined_met %>%
  group_by(DATE, SITECODE) %>%
  summarise(
    P_mm_d    = mean(P_mm_d,    na.rm = TRUE),
    T_C       = mean(T_C,       na.rm = TRUE),
    RH_d_pct  = mean(RH_d_pct,  na.rm = TRUE),
    NR_Wm2_d  = mean(NR_Wm2_d,  na.rm = TRUE),
    .groups = "drop"
  )

# 2) pivot wide → long
long_met <- combined_unique %>%
  pivot_longer(
    cols      = c(P_mm_d, T_C, RH_d_pct, NR_Wm2_d),
    names_to  = "var",
    values_to = "value"
  )
#  long_met: DATE | SITECODE | var     | value

# 3) define your station-selection mapping:
#    for each watershed×variable, which site to pull from
mapping <- tribble(
  ~watershed, ~var,       ~SITECODE,
  "WS01",     "P_mm_d",   "PRIMET",
  "WS01",     "T_C",      "PRIMET",
  "WS01",     "NR_Wm2_d", "PRIMET",
  "WS01",     "RH_d_pct", "CENMET",
  # add rows here for WS02, WS03, etc...
)

# 4) join mapping to the long data
ws_long <- mapping %>%
  left_join(long_met, by = c("SITECODE","var"))
#  ws_long: watershed | var | SITECODE | DATE | value

# 5) pivot long → wide per watershed
combined_by_watershed <- ws_long %>%
  select(watershed, DATE, var, value) %>%
  pivot_wider(
    names_from  = var,
    values_from = value
  ) %>%
  arrange(watershed, DATE)