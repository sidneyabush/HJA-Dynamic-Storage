library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)

rm(list = ls())

base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

da_df <- read_csv(file.path(base_dir, "Q", "drainage_area.csv"))
discharge <- read_csv(file.path(base_dir, "Q", "HF00402_v14.csv")) %>%
  filter(WATERYEAR > 1997) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    Q    = MEAN_Q * 0.02831683199881      # convert cfs to m³/s
  ) %>%
  arrange(SITECODE, Date)

calc_recession <- function(df) {
  tmp <- df %>%
    arrange(Date) %>%
    mutate(
      dQ           = Q - lag(Q),
      dQ_dt        = dQ / as.numeric(Date - lag(Date)),
      change_ratio = Q / lag(Q)
    ) %>%
    filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
    mutate(recession_slope = -dQ_dt)
  fit <- lm(recession_slope ~ Q, data = tmp)
  sm  <- summary(fit)
  tibble(
    slope     = coef(fit)[2],
    p_value   = sm$coefficients[2, 4],
    r_squared = sm$r.squared
  )
}

calc_rbfi <- function(df) {
  tmp <- df %>%
    arrange(Date) %>%
    mutate(
      dQ     = Q - lag(Q),
      abs_dQ = abs(dQ)
    ) %>%
    filter(!is.na(abs_dQ))
  tibble(RBFI = sum(tmp$abs_dQ, na.rm = TRUE) / sum(tmp$Q, na.rm = TRUE))
}

annual_metrics <- discharge %>%
  group_by(SITECODE, WATERYEAR) %>%
  group_map(~ bind_cols(
    tibble(site = .y$SITECODE, year = .y$WATERYEAR),
    calc_recession(.x),
    calc_rbfi(.x)
  )) %>%
  bind_rows()

p_rbfi <- ggplot(annual_metrics, aes(x = year, y = RBFI)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ site, ncol = 3, scales = "fixed") +
  labs(
    x = "Water Year",
    y = "Richards–Baker Flashiness Index"
  ) +
  theme_classic()

ggsave(
  filename = "rbfi_by_site_wy.png",
  plot     = p_rbfi,
  path     = output_dir,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)

recession_data_all <- discharge %>%
  arrange(SITECODE, Date) %>%
  mutate(
    dQ           = Q - lag(Q),
    dQ_dt        = dQ / as.numeric(Date - lag(Date)),
    change_ratio = Q / lag(Q)
  ) %>%
  filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
  mutate(
    recession_slope = -dQ_dt,
    site            = SITECODE
  )

p_recession <- ggplot(recession_data_all, aes(x = Q, y = recession_slope)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ site, ncol = 3, scales = "free") +
  labs(
    x = expression(Discharge~(m^3/s)),
    y = expression(Recession~Slope~(-dQ/dt))
  ) +
  theme_classic()

ggsave(
  filename = "recession_curve_by_site.png",
  plot     = p_recession,
  path     = output_dir,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)
