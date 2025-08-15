# --- STEP 3 & 4: MERGE + PLOT ALL ET METHODS ----------------------

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# --- 0. DIRECTORIES & PALETTES -------------------------------------
input_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"
plot_dir    <- file.path(input_dir, "plots", "ET_methods_comparison")
scatter_dir <- file.path(plot_dir, "scatter")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)

method_colors <- c(
  "PT-Zhang (2024)"                      = "#0173B2",
  "PT-Szilagyi (2014)"                   = "#DE8F05",
  "Hamon Uncalibrated"                   = "#C3D7A4",
  "Hamon (Monthly Calibrated, Zhang)"    = "#029E73",
  "Hamon (Monthly Calibrated, Szilagyi)" = "#D55E00",
  "Hamon (Interpolated, Zhang)"          = "#56B4E9",
  "Hamon (Interpolated, Szilagyi)"       = "#CC79A7"
)

method_label_to_col <- c(
  "PT-Zhang (2024)"                      = "ET_PT_zhang",
  "PT-Szilagyi (2014)"                   = "ET_PT_szilagyi",
  "Hamon Uncalibrated"                   = "ET_Hamon_uncalibrated",
  "Hamon (Monthly Calibrated, Zhang)"    = "ET_Hamon_pt_zhang_monthly",
  "Hamon (Monthly Calibrated, Szilagyi)" = "ET_Hamon_pt_szilagyi_monthly",
  "Hamon (Interpolated, Zhang)"          = "ET_Hamon_pt_zhang_interp_full",
  "Hamon (Interpolated, Szilagyi)"       = "ET_Hamon_pt_szilagyi_interp_full"
)

# --- 1. LOAD MERGED & EXPORT FOR SHARING ---------------------------
df <- read_csv(
  file.path(input_dir, "daily_water_balance_all_ET_methods_1997_present.csv"),
  show_col_types = FALSE
)
write_csv(df,
          file.path(input_dir, "daily_water_balance_all_ET_methods_1997_present.csv"))

stopifnot(all(method_label_to_col %in% names(df)))
site_list <- unique(df$SITECODE)

# helper to pivot any subset of methods
pivot_methods <- function(df, methods){
  df %>%
    select(DATE, SITECODE, all_of(methods)) %>%
    pivot_longer(-c(DATE, SITECODE), names_to="Method", values_to="ET_mm_day")
}

# --- 2. SCATTERPLOTS (2013–2019 calibration window, 4 panels) ------
for(site in site_list){
  d <- df %>%
    filter(SITECODE == site,
           between(DATE, as.Date("2013-01-01"), as.Date("2019-12-31")))
  
  # Panel 1: Monthly‐Calibrated Zhang
  p1 <- ggplot(d, aes(x = ET_PT_zhang)) +
    geom_point(aes(y = ET_Hamon_pt_zhang_monthly,
                   color = "Hamon (Monthly Calibrated, Zhang)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = method_colors["Hamon (Monthly Calibrated, Zhang)"]) +
    labs(title=paste(site, "- Zhang (monthly)"),
         x="PT–Zhang ET (mm/day)", y="Hamon ET (mm/day)") +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  # Panel 2: Interpolated Zhang
  p2 <- ggplot(d, aes(x = ET_PT_zhang)) +
    geom_point(aes(y = ET_Hamon_pt_zhang_interp_full,
                   color = "Hamon (Interpolated, Zhang)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = method_colors["Hamon (Interpolated, Zhang)"]) +
    labs(title=paste(site, "- Zhang (interpolated)"),
         x="PT–Zhang ET (mm/day)", y=NULL) +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  # Panel 3: Monthly‐Calibrated Szilagyi
  p3 <- ggplot(d, aes(x = ET_PT_szilagyi)) +
    geom_point(aes(y = ET_Hamon_pt_szilagyi_monthly,
                   color = "Hamon (Monthly Calibrated, Szilagyi)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = method_colors["Hamon (Monthly Calibrated, Szilagyi)"]) +
    labs(title=paste(site, "- Szilagyi (monthly)"),
         x="PT–Szilagyi ET (mm/day)", y="Hamon ET (mm/day)") +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  # Panel 4: Interpolated Szilagyi
  p4 <- ggplot(d, aes(x = ET_PT_szilagyi)) +
    geom_point(aes(y = ET_Hamon_pt_szilagyi_interp_full,
                   color = "Hamon (Interpolated, Szilagyi)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = method_colors["Hamon (Interpolated, Szilagyi)"]) +
    labs(title=paste(site, "- Szilagyi (interpolated)"),
         x="PT–Szilagyi ET (mm/day)", y=NULL) +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  # Combine into 2×2 with shared legend
  combo <- (p1 + p2) / (p3 + p4) +
    plot_layout(guides="collect") &
    theme(legend.position="bottom",
          legend.title=element_blank())
  
  ggsave(file.path(scatter_dir, paste0("scatter_calibration_4panel_", site, ".png")),
         combo, width=16, height=12)
}

# --- 3a. PER-SITE FULL-RECORD TIME SERIES --------------------------
for(site in site_list){
  ts_long <- pivot_methods(df, method_label_to_col) %>%
    filter(SITECODE == site)
  
  p_ts <- ggplot(ts_long, aes(DATE, ET_mm_day, color=Method)) +
    geom_line(size=0.6, alpha=0.8) +
    scale_color_manual(values=method_colors) +
    labs(title=paste("Daily ET at", site),
         x="Date", y="ET (mm/day)", color="Method") +
    theme_bw(base_size=13) +
    theme(legend.position="bottom",
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.minor=element_blank())
  
  ggsave(file.path(plot_dir, paste0("full_record_", site, ".png")),
         p_ts, width=12, height=6)
}

# --- 3b. GRID OF ALL SITES, FULL RECORD ---------------------------
all_long <- pivot_methods(df, method_label_to_col) %>%
  filter(DATE >= as.Date("1997-01-01"))

p_grid_full <- ggplot(all_long, aes(DATE, ET_mm_day, color=Method)) +
  geom_line(size=0.4, alpha=0.7) +
  facet_wrap(~SITECODE, scales="free_y", ncol=3) +
  scale_color_manual(values=method_colors) +
  labs(title="ET Methods Across All Sites (1997–Present)",
       x="Date", y="ET (mm/day)", color="Method") +
  theme_bw(base_size=12) +
  theme(strip.text       = element_text(face="bold", size=9),
        axis.text.x      = element_text(angle=45,hjust=1,size=7),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "grid_all_sites_full_record.png"),
       p_grid_full, width=16, height=10)

# --- 3c. GRID BY METHOD, FULL RECORD ------------------------------
for(label in names(method_label_to_col)){
  col <- method_label_to_col[[label]]
  d_m <- df %>%
    transmute(DATE, SITECODE, ET_mm_day = .data[[col]]) %>%
    filter(!is.na(ET_mm_day))
  
  p_m <- ggplot(d_m, aes(DATE, ET_mm_day)) +
    geom_line(color=method_colors[label], size=0.5, alpha=0.8) +
    facet_wrap(~SITECODE, scales="free_y", ncol=3) +
    labs(title=label, x="Date", y="ET (mm/day)") +
    theme_bw(base_size=12) +
    theme(strip.text       = element_text(face="bold", size=9),
          axis.text.x      = element_text(angle=45,hjust=1,size=7),
          panel.grid.minor = element_blank())
  
  fname <- gsub("[^A-Za-z0-9]+","_", label)
  ggsave(file.path(plot_dir, paste0("grid_by_method_", fname, ".png")),
         p_m, width=14, height=10)
}

# --- 4a. PER-SITE CALIBRATION-WINDOW (2013–2019) -------------------
for(site in site_list){
  cal_long <- pivot_methods(df, method_label_to_col) %>%
    filter(SITECODE == site,
           between(DATE, as.Date("2013-01-01"), as.Date("2019-12-31")))
  
  p_cal_ts <- ggplot(cal_long, aes(DATE, ET_mm_day, color=Method)) +
    geom_line(size=0.6, alpha=0.8) +
    scale_color_manual(values=method_colors) +
    labs(title=paste("Calibration Window at", site),
         x="Date", y="ET (mm/day)", color="Method") +
    theme_bw(base_size=13) +
    theme(legend.position="bottom",
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.minor=element_blank())
  
  ggsave(file.path(plot_dir, paste0("cal_window_", site, ".png")),
         p_cal_ts, width=12, height=6)
}

# --- 4b. GRID ALL SITES, CALIBRATION WINDOW -----------------------
cal_grid <- all_long %>%
  filter(between(DATE, as.Date("2013-01-01"), as.Date("2019-12-31")))

p_grid_cal <- ggplot(cal_grid, aes(DATE, ET_mm_day, color=Method)) +
  geom_line(size=0.4, alpha=0.7) +
  facet_wrap(~SITECODE, scales="free_y", ncol=3) +
  scale_color_manual(values=method_colors) +
  labs(title="ET Methods Across All Sites (2013–2019)",
       x="Date", y="ET (mm/day)", color="Method") +
  theme_bw(base_size=12) +
  theme(strip.text       = element_text(face="bold", size=9),
        axis.text.x      = element_text(angle=45,hjust=1,size=7),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "grid_all_sites_2013_2019.png"),
       p_grid_cal, width=16, height=10)
