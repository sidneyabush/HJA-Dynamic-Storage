# ---------------------------------------------
# Hydrometric Storage & Discharge‐Storage Workflow
# ---------------------------------------------

# 0) Load libraries and set theme
library(dplyr)
library(readr)
library(ggplot2)
library(colorspace)
library(tidyr)
library(multcompView)
library(lubridate)

theme_set(theme_classic(base_size = 16))

# Clear environment
rm(list = ls())

# Directories
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/Hydrometric"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Site order & palette
site_order <- c(
  "GSWS09","GSWS10","GSLOOK","GSWS01","GSWS02",
  "GSWS03","GSMACK","GSWS06","GSWS07","GSWS08"
)
palette10  <- c(
  "#AA4499", "#882255", "#CC6677", "#DDCC77", "#999933",
  "#117733", "#44AA99", "#88CCEE", "#6699CC", "#332288"
)
site_cols <- setNames(lighten(palette10, amount = 0.1), site_order)

# Axis‐label mappings
axis_labels <- c(
  fdc_slope             = "Flow Duration Curve Slope",
  mean_bf               = "Mean Baseflow",
  Q5norm                = "Normalized Q5",
  CV_Q5norm             = "Normalized CV_Q5",
  RBI                   = "Richards–Baker Index",
  recession_curve_slope = "Recession-curve Slope",
  S_annual_mm           = "DS Storage-Discharge",
  DS_sum                = "DS Drawdown"
)

# ------------------------------------------------------------------
# 1) Read Storage Metrics and reshape to long
# ------------------------------------------------------------------
hydrometric_storage <- read_csv(
  file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual.csv"),
  show_col_types = FALSE
) %>%
  select(-`...1`) %>%
  rename(RBI = rbfi) %>%
  # recode GSWSMC → GSMACK
  mutate(
    site = if_else(site == "GSWSMC", "GSMACK", site)
  ) %>%
  filter(site %in% site_order) %>%
  mutate(site = factor(site, levels = site_order))

hydrometric_storage_long <- hydrometric_storage %>%
  pivot_longer(
    cols      = -c(site, year),
    names_to  = "metric",
    values_to = "value"
  )

# ------------------------------------------------------------------
# 2) Read Dynamic Storage Discharge CSV and compute water-year
# ------------------------------------------------------------------
ds_discharge_file <- file.path(
  base_dir, "DynamicStorage",
  "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"
)

ds_discharge_long <- read_csv(ds_discharge_file, show_col_types = FALSE) %>%
  # rename wateryear → year and keep the discharge column
  select(
    site = site,
    year = wateryear,
    S_annual_mm
  ) %>%
  # recode GSLOOK_FULL → GSLOOK
  mutate(
    site = if_else(site == "GSLOOK_FULL", "GSLOOK", site)
  ) %>%
  # pivot into long form
  pivot_longer(
    cols      = S_annual_mm,
    names_to  = "metric",
    values_to = "value"
  ) %>%
  # ensure site is a factor in the proper order
  mutate(site = factor(site, levels = site_order))

# ------------------------------------------------------------------
# 3) Read DS Drawdown .csv
# ------------------------------------------------------------------
ds_draw_file  <- file.path(
  base_dir, "DynamicStorage", "DS_drawdown_annual.csv")


ds_draw_long <- read_csv(ds_draw_file, show_col_types = FALSE) %>%
  select(-...1) %>%
  # columns are: SITECODE, waterYear, DS_sum
  rename(
    site   = SITECODE,
    year   = waterYear,
    DS_sum = DS_sum
  ) %>%
  # recode GSWSMA → GSMACK
  mutate(
    site = if_else(site == "GSWSMA", "GSMACK", site)
  ) %>%
  pivot_longer(
    cols      = DS_sum,
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(site = factor(site, levels = site_order))


# ------------------------------------------------------------------
# 4) Combine all metrics into one long table
# ------------------------------------------------------------------
storage_long <- bind_rows(hydrometric_storage_long, ds_discharge_long, ds_draw_long)

# ------------------------------------------------------------------
# 5) Per‐metric time‐series & summary plots (loop over every metric)
# ------------------------------------------------------------------
for (m in unique(storage_long$metric)) {
  df         <- filter(storage_long, metric == m)
  axis_label <- axis_labels[m]
  
  # time‐series by site
  p_ts <- ggplot(df, aes(year, value, color = site, fill = site)) +
    geom_line(size = 0.6) +
    geom_point(size = 1) +
    facet_wrap(~ site, ncol = 2) +
    scale_color_manual(values = site_cols, guide = FALSE) +
    scale_fill_manual(values = alpha(site_cols, 0.3), guide = FALSE) +
    labs(x = "Water Year", y = axis_label) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, size = 1),
      axis.line        = element_blank(),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(color = "black", hjust = 0)
    )
  ggsave(
    paste0("storage_", m, "_by_site_wy.png"),
    p_ts, path = output_dir, width = 8, height = 9, units = "in", dpi = 300
  )
  
  # summary (mean ± SD + Tukey letters)
  sum_df <- df %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value,   na.rm = TRUE),
      .groups  = "drop"
    )
  aov_r <- aov(value ~ site, data = df)
  tuk   <- TukeyHSD(aov_r, "site")$site
  pvals <- setNames(tuk[, "p adj"], rownames(tuk))
  lets  <- multcompLetters(pvals)$Letters
  sum_df$group <- lets[as.character(sum_df$site)]
  
  y_max  <- max(sum_df$mean_val + sum_df$sd_val, na.rm = TRUE)
  y_min  <- min(sum_df$mean_val - sum_df$sd_val, na.rm = TRUE)
  offset <- 0.05 * (y_max - y_min)
  sum_df <- sum_df %>% mutate(label_y = mean_val + sd_val + offset)
  
  p_sum <- ggplot(sum_df, aes(site, mean_val, color = site)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  width = 0.2) +
    geom_text(aes(label = group, y = label_y), size = 4, vjust = 0) +
    scale_x_discrete(limits = site_order) +
    scale_color_manual(values = site_cols, guide = FALSE) +
    labs(x = NULL, y = paste0("Mean ±1 SD of ", axis_label)) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.border     = element_rect(color = "black", fill = NA, size = 1),
      axis.line        = element_blank()
    )
  ggsave(
    paste0("storage_summary_", m, ".png"),
    p_sum, path = output_dir, width = 6, height = 4, units = "in", dpi = 300
  )
}

# ------------------------------------------------------------------
# 6) Grid of ALL methods (mean ± SD + letters)
# ------------------------------------------------------------------
summary_all <- storage_long %>%
  group_by(metric, site) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value,   na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  group_by(metric) %>%
  do({
    dat   <- .
    dfm   <- filter(storage_long, metric == dat$metric[1])
    aov_r <- aov(value ~ site, data = dfm)
    tuk_r <- TukeyHSD(aov_r, "site")$site
    pv    <- setNames(tuk_r[, "p adj"], rownames(tuk_r))
    let   <- multcompLetters(pv)$Letters
    dat$group   <- let[as.character(dat$site)]
    y_mx        <- max(dat$mean_val + dat$sd_val, na.rm = TRUE)
    y_mn        <- min(dat$mean_val - dat$sd_val, na.rm = TRUE)
    dat$label_y <- dat$mean_val + dat$sd_val + 0.05*(y_mx - y_mn)
    dat
  }) %>%
  ungroup()

p_grid <- ggplot(summary_all, aes(site, mean_val, color = site)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                width = 0.2) +
  geom_text(aes(label = group, y = label_y), size = 3.5, vjust = 0) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y",
             labeller = as_labeller(axis_labels)) +
  scale_x_discrete(limits = site_order) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(x = NULL, y = "Mean ± 1 SD") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black", hjust = 0)
  )

ggsave(
  "grid_all_methods.png",
  plot  = p_grid,
  path  = output_dir,
  width = 9, height = 11, units = "in", dpi = 300
)
