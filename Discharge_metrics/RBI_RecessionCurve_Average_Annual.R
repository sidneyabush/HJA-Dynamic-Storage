library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(scales)
library(colorspace)
library(tidyr)

theme_set(theme_classic(base_size = 14))

# Clear environment
rm(list = ls())

# Directories
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Only these 10 sites
sites_keep <- c(
  "GSLOOK","GSWS01","GSWS02","GSWS03","GSWS09",
  "GSWS10","GSWSMC","GSWS06","GSWS07","GSWS08"
)

# Read & prep data
da_df <- read_csv(file.path(base_dir, "Q", "drainage_area.csv"))
discharge <- read_csv(file.path(base_dir, "Q", "HF00402_v14.csv")) %>%
  #filter(WATERYEAR > 1997, SITECODE %in% sites_keep) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    Q    = MEAN_Q * 0.02831683199881  # m³/s
  ) %>%
  arrange(SITECODE, Date)

# Need to update this function to take the log-log: 
# lm_model <- lm(log(recession_slope) ~ log(Q), data = recession_data)


# 1) Recession slope
calc_recession <- function(df) {
  tmp <- df %>%
    mutate(
      dQ           = Q - lag(Q),
      dQ_dt        = dQ / as.numeric(Date - lag(Date)),
      change_ratio = Q / lag(Q)
    ) %>%
    filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
    mutate(recession_slope = -dQ_dt)
  
  tibble(slope = coef(lm(log(recession_slope) ~ log(Q), data = tmp))[2])
}

# 2) RBI
calc_RBI <- function(df) {
  tmp     <- df %>% mutate(dQ = Q - lag(Q)) %>% filter(!is.na(dQ))
  total_Q <- sum(df$Q, na.rm = TRUE)
  props   <- abs(tmp$dQ) / total_Q
  
  tibble(RBI = sum(props, na.rm = TRUE))
}

# Compute annual metrics
annual_metrics <- discharge %>%
  group_by(SITECODE, WATERYEAR) %>%
  group_map(~ bind_cols(
    tibble(site = .y$SITECODE, year = .y$WATERYEAR),
    calc_recession(.x),
    calc_RBI(.x)
  )) %>%
  bind_rows() %>%
  mutate(site = factor(site, levels = sites_keep))

# Prepare full‐record recession data in mm/day
recession_clean <- discharge %>%
  mutate(
    dQ             = Q - lag(Q),
    dQ_dt          = dQ / as.numeric(Date - lag(Date)),        # m³/s per day
    recession_slope = -dQ_dt,                                   # m³/s per day
    # convert to mm/day: (m³/s) ÷ basin area (m²) × 86400 s/day × 1000 mm/m
    Q_mm_day       = Q / DA_M2 * 86400 * 1000,                 
    slope_mm_day   = recession_slope / DA_M2 * 86400 * 1000    
  ) %>%
  filter(!is.na(Q_mm_day), !is.na(slope_mm_day)) %>%
  transmute(
    site          = factor(SITECODE, levels = sites_keep),
    Q_mm_day,
    slope_mm_day
  )

# Annotation positions (force positive slopes)
slopes_df <- recession_clean %>%
  group_by(site) %>%
  summarise(
    slope_mm     = abs(round(coef(lm(slope_mm_day ~ Q_mm_day))[2], 2)),
    Q_pos_mm     = quantile(Q_mm_day,     0.9,  na.rm = TRUE),
    slope_pos_mm = quantile(slope_mm_day, 0.85, na.rm = TRUE),
    .groups      = "drop"
  )

# Color palette → lighten by 10%
palette10 <- c(
  "#AA4499", "#882255", "#CC6677", "#DDCC77", "#999933",
  "#117733", "#44AA99", "#88CCEE", "#6699CC", "#332288"
)
site_cols <- setNames(lighten(palette10, amount = 0.1), sites_keep)

# 1) Annual RBI by site
p_RBI <- ggplot(annual_metrics, aes(year, RBI, color = site, fill = site)) +
  geom_line(size = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ site, ncol = 2) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  scale_fill_manual(values = alpha(site_cols, 0.3), guide = FALSE) +
  labs(x = "Water Year", y = "Richards–Baker Index") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )
ggsave("RBI_by_site_wy.png", p_RBI, path = output_dir,
       width = 8, height = 9, units = "in", dpi = 300)

# 2) Annual recession‐curve slope by site
p_slope <- ggplot(annual_metrics, aes(year, slope, color = site, fill = site)) +
  geom_line(size = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ site, ncol = 2) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  scale_fill_manual(values = alpha(site_cols, 0.3), guide = FALSE) +
  labs(x = "Water Year", y = "Recession Limb Slope") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )
ggsave("recession_slope_by_site_wy.png", p_slope, path = output_dir,
       width = 8, height = 9, units = "in", dpi = 300)

# 3) Instantaneous recession‐limb curves: both axes in mm/day with slope annotation
p_curve <- ggplot(recession_clean, aes(x = Q_mm_day, y = abs(slope_mm_day), color = site)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  geom_text(
    data        = slopes_df,
    inherit.aes = FALSE,
    aes(
      x     = Q_pos_mm,
      y     = slope_pos_mm,
      label = paste0("slope = ", slope_mm)
    ),
    hjust = 0.3,
    vjust = -9,
    size  = 4,
    color = "black"
  ) +
  facet_wrap(~ site, ncol = 2, scales = "fixed") +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(
    x = "Q (mm day⁻¹)",
    y = "–dQ/dt (mm day⁻¹)"
  ) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )
ggsave("recession_curve_by_site_mmday.png", p_curve, path = output_dir,
       width = 8, height = 9, units = "in", dpi = 300)


summary_site <- annual_metrics %>%
  select(site, RBI, slope) %>%
  pivot_longer(
    cols        = c(RBI, slope),
    names_to    = "metric",
    values_to   = "value"
  ) %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value,   na.rm = TRUE),
    .groups  = "drop"
  )

# Replace your facet_wrap line with this to set custom panel titles:
p_summary <- ggplot(summary_site, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                width = 0.2) +
  facet_wrap(
    ~ metric,
    scales   = "free_y",
    ncol     = 1,     # stack panels vertically
    labeller = labeller(metric = c(
      RBI  = "Richards–Baker Index",
      slope = "Recession Limb Slope"
    ))
  ) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(x = NULL, y = "Mean ± 1 SD") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )


ggsave("MeanSD_RBI_RecessionSlope_by_site.png", p_summary, path = output_dir,
       width = 5, height = 8, units = "in", dpi = 300)
