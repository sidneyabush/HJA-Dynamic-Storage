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

map_df <- tribble(
  ~SITECODE, ~watershed,
  "PRIMET",  "GSWS01",
  "PRIMET",  "GSWS02",
  "PRIMET",  "GSWS03"
  # add more mappings for other groups here
)

# 2) Join & expand
combined_expanded <- combined_met %>%
  # join on SITECODE; each PRIMET row will now match three rows in map_df
  left_join(map_df, by = "SITECODE") %>%
  # if you only want the ones you mapped, drop the rest:
  filter(!is.na(watershed)) %>%
  arrange(DATE, SITECODE, watershed)