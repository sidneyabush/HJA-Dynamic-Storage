library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(scales)

rm(list = ls())

base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# define which sites to include
sites_keep <- c("GSLOOK","GSWS01","GSWS02","GSWS03","GSWSMC","GSWS06","GSWS07","GSWS08")

# Read in drainage area and discharge, but drop GSWS09/10 immediately
da_df <- read_csv(file.path(base_dir, "Q", "drainage_area.csv"))
discharge <- read_csv(file.path(base_dir, "Q", "HF00402_v14.csv")) %>%
  filter(WATERYEAR > 1997, SITECODE %in% sites_keep) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    Q    = MEAN_Q * 0.02831683199881
  ) %>%
  arrange(SITECODE, Date)

# Functions for recession slope and RBFI
calc_recession <- function(df) {
  tmp <- df %>%
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
    mutate(dQ = Q - lag(Q), abs_dQ = abs(dQ)) %>%
    filter(!is.na(abs_dQ))
  tibble(RBFI = sum(tmp$abs_dQ) / sum(tmp$Q))
}

# Compute annual metrics per site × water year
annual_metrics <- discharge %>%
  group_by(SITECODE, WATERYEAR) %>%
  group_map(~ bind_cols(
    tibble(site = .y$SITECODE, year = .y$WATERYEAR),
    calc_recession(.x),
    calc_rbfi(.x)
  )) %>%
  bind_rows()

# Prepare full‐record recession data
recession_clean <- discharge %>%
  mutate(
    dQ           = Q - lag(Q),
    dQ_dt        = dQ / as.numeric(Date - lag(Date)),
    change_ratio = Q / lag(Q)
  ) %>%
  filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
  transmute(
    site            = SITECODE,
    Q,
    recession_slope = -dQ_dt
  )

# Compute annotation slopes & positions
slopes_df <- recession_clean %>%
  group_by(site) %>%
  summarise(
    slope     = round(coef(lm(recession_slope ~ Q))[2], 2),
    Q_pos     = quantile(Q,             0.9, na.rm = TRUE),
    dQ_dt_pos = quantile(recession_slope, 0.85, na.rm = TRUE),
    .groups   = "drop"
  )

# Site colors & manual ordering 
palette8 <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermilion
  "#CC79A7",  # reddish purple
  "#999999"   # grey
)
site_cols <- setNames(palette8, sites_keep)

# enforce factor ordering
annual_metrics <- annual_metrics %>% mutate(site = factor(site, levels = sites_keep))
recession_clean <- recession_clean   %>% mutate(site = factor(site, levels = sites_keep))
slopes_df$site <- factor(slopes_df$site, levels = sites_keep)

# 1) Annual RBFI: 4 columns × 2 rows (leave empty for last two)
p_rbfi <- ggplot(annual_metrics, aes(year, RBFI, color = site, group = site)) +
  geom_point(size = 1.1) +
  geom_line(size = 0.7) +
  facet_wrap(~ site, ncol = 4) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(x = "Water Year", y = "Richards–Baker Flashiness Index") +
  theme_classic() +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(color = "black", fill = "white", size = 0.5),
    strip.text       = element_text(color = "black")
  )
ggsave("rbfi_by_site_wy.png", p_rbfi, path = output_dir,
       width = 12, height = 6, units = "in", dpi = 300)

# 2) Recession curve: 5 columns × 2 rows
p_curve <- ggplot(recession_clean, aes(Q, recession_slope, color = site)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.8) +
  geom_text(
    data = slopes_df,
    aes(label = paste0("slope = ", slope)),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.6,
    size = 3, color = "black"
  ) +
  facet_wrap(~ site, ncol = 4, scales = "free") +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(
    x = expression(Discharge~(m^3/s)),
    y = expression(Recession~Slope~(-dQ/dt))
  ) +
  theme_classic() +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(color = "black", fill = "white", size = 0.5),
    strip.text       = element_text(color = "black")
  )
ggsave("recession_curve_by_site.png", p_curve, path = output_dir,
       width = 12, height = 6, units = "in", dpi = 300)

# 3) Annual recession slope vs WY: 5 columns × 2 rows
p_slope <- ggplot(annual_metrics, aes(year, slope, color = site, group = site)) +
  geom_point(size = 1.1) +
  geom_line(size = 0.7) +
  facet_wrap(~ site, ncol = 4) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(
    x = "Water Year",
    y = expression(Recession~Slope~(-dQ/dt~vs~Q))
  ) +
  theme_classic() +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(color = "black", fill = "white", size = 0.5),
    strip.text       = element_text(color = "black")
  )
ggsave("recession_limb_by_site_wy.png", p_slope, path = output_dir,
       width = 12, height = 6, units = "in", dpi = 300)
