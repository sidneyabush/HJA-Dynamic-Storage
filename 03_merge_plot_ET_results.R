################################################################################
## 03. MERGE & PLOT COMPARISONS (CALIBRATION, TIME SERIES, GRID)
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# Directories
input_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"
plot_dir    <- file.path(input_dir, "plots", "ET_methods_comparison")
scatter_dir <- file.path(plot_dir, "scatter")
dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load merged monthly‐calibrated Hamon (includes all PT columns)
merged <- read_csv(file.path(input_dir, "Hamon_ET_methods_timeseries_monthly_calibrated.csv"),
                   show_col_types = FALSE)

site_list <- unique(merged$SITECODE)

# 2. Define palette & mapping
method_colors <- c(
  "PT-Zhang (2024)"                   = "#0173B2",
  "PT-Szilagyi (2014)"                = "#DE8F05",
  "Hamon Uncalibrated"                = "#C3D7A4",
  "Hamon (Monthly Calibrated, Zhang)" = "#029E73",
  "Hamon (Monthly Calibrated, Szilagyi)" = "#D55E00"
)
method_label_to_col <- c(
  "PT-Zhang (2024)"                   = "ET_PT_zhang",
  "PT-Szilagyi (2014)"                = "ET_PT_szilagyi",
  "Hamon Uncalibrated"                = "ET_Hamon_uncalibrated",
  "Hamon (Monthly Calibrated, Zhang)" = "ET_Hamon_pt_zhang_monthly",
  "Hamon (Monthly Calibrated, Szilagyi)" = "ET_Hamon_pt_szilagyi_monthly"
)

# 3a. Calibration scatterplots (2013–2019)
for(site in site_list){
  d <- merged %>%
    filter(SITECODE == site,
           DATE >= as.Date("2013-01-01"),
           DATE <= as.Date("2019-12-31")) %>%
    transmute(
      PT_Zhang           = ET_PT_zhang,
      PT_Szilagyi        = ET_PT_szilagyi,
      Hamon_Uncalibrated = ET_Hamon_uncalibrated,
      Hamon_Cal_Zhang    = ET_Hamon_pt_zhang_monthly,
      Hamon_Cal_Szilagyi = ET_Hamon_pt_szilagyi_monthly
    )
  
  p1 <- ggplot(d, aes(x = PT_Zhang)) +
    geom_point(aes(y = Hamon_Uncalibrated, color = "Hamon Uncalibrated"),
               alpha = 0.7, size = 1) +
    geom_point(aes(y = Hamon_Cal_Zhang, color = "Hamon Calibrated"),
               alpha = 0.7, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(
      "Hamon Uncalibrated" = "#C3D7A4",
      "Hamon Calibrated"   = "#029E73"
    )) +
    labs(title = paste(site, "– vs PT-Zhang"),
         x = "PT-Zhang ET (mm/day)",
         y = "Hamon ET (mm/day)",
         color = "") +
    theme_bw() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          panel.grid           = element_blank())
  
  p2 <- ggplot(d, aes(x = PT_Szilagyi)) +
    geom_point(aes(y = Hamon_Uncalibrated, color = "Hamon Uncalibrated"),
               alpha = 0.7, size = 1) +
    geom_point(aes(y = Hamon_Cal_Szilagyi, color = "Hamon Calibrated"),
               alpha = 0.7, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(
      "Hamon Uncalibrated" = "#C3D7A4",
      "Hamon Calibrated"   = "#D55E00"
    )) +
    labs(title = paste(site, "– vs PT-Szilagyi"),
         x = "PT-Szilagyi ET (mm/day)",
         y = NULL, color = "") +
    theme_bw() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          panel.grid           = element_blank())
  
  combined <- p1 + p2 + plot_layout(ncol = 2, guides = "collect")
  
  ggsave(
    filename = file.path(scatter_dir, paste0("scatter_calib_", site, ".png")),
    plot     = combined,
    width    = 14,
    height   = 7
  )
}

# 3b. Time series plots by site (full record)
for(site in site_list){
  ts_long <- merged %>%
    filter(SITECODE == site) %>%
    select(DATE, all_of(method_label_to_col)) %>%
    pivot_longer(-DATE, names_to = "Method", values_to = "ET_mm_day") %>%
    mutate(Method = factor(Method,
                           levels = method_label_to_col,
                           labels = names(method_label_to_col)))
  
  p_ts <- ggplot(ts_long, aes(x = DATE, y = ET_mm_day, color = Method)) +
    geom_line(size = 0.6) +
    scale_color_manual(values = method_colors) +
    labs(title = paste("Daily ET at", site),
         x = "Date", y = "ET (mm/day)", color = "Method") +
    theme_bw() +
    theme(legend.position      = "bottom",
          axis.text.x          = element_text(angle = 45, hjust = 1),
          panel.grid.minor     = element_blank())
  
  ggsave(
    filename = file.path(plot_dir, paste0("TS_", site, ".png")),
    plot     = p_ts,
    width    = 14,
    height   = 8
  )
}

# 3c. Grid plots: each method across all sites
for(lbl in names(method_label_to_col)){
  col   <- method_label_to_col[[lbl]]
  grid  <- merged %>%
    filter(!is.na(.data[[col]])) %>%
    transmute(DATE, SITECODE, ET_mm_day = .data[[col]])
  
  p_grid <- ggplot(grid, aes(x = DATE, y = ET_mm_day)) +
    geom_line(color = method_colors[lbl]) +
    facet_wrap(~ SITECODE, scales = "free_y", ncol = 3) +
    labs(title = lbl, x = "Date", y = "ET (mm/day)") +
    theme_bw() +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  
  ggsave(
    filename = file.path(plot_dir, paste0("grid_", gsub(" ", "_", lbl), ".png")),
    plot     = p_grid,
    width    = 14,
    height   = 10
  )
}
