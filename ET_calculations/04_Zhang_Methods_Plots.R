# --- STEP 5: PLOT ONLY ZHANG‐BASED METHODS ------------------------

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# --- 0. DIRECTORIES & PALETTES -------------------------------------
input_dir        <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/ET"
plot_dir         <- file.path(input_dir, "plots", "ET_methods_comparison")
zhang_dir        <- file.path(plot_dir, "zhang_methods")
zhang_scatter    <- file.path(zhang_dir, "scatter")
dir.create(zhang_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(zhang_scatter, showWarnings = FALSE, recursive = TRUE)

# Only Zhang‐related methods
zhang_colors <- c(
  "PT-Zhang (2024)"                   = "#66C2A5",  
  "Hamon Uncalibrated"                = "#FC8D62",  
  "Hamon (Monthly Calibrated, Zhang)" = "#8DA0CB",  
  "Hamon (Interpolated, Zhang)"       = "#E78AC3"   
)

zhang_methods <- c(
  "PT-Zhang (2024)"                   = "ET_PT_zhang",
  "Hamon Uncalibrated"                = "ET_Hamon_uncalibrated",
  "Hamon (Monthly Calibrated, Zhang)" = "ET_Hamon_pt_zhang_monthly",
  "Hamon (Interpolated, Zhang)"       = "ET_Hamon_pt_zhang_interp_full"
)

# --- 1. LOAD MERGED DATA -------------------------------------------
df <- read_csv(
  file.path(input_dir, "daily_water_balance_all_ET_methods_1997_present.csv"),
  show_col_types = FALSE
)
stopifnot(all(zhang_methods %in% names(df)))
site_list <- unique(df$SITECODE)

pivot_methods <- function(df, methods){
  df %>%
    select(DATE, SITECODE, all_of(methods)) %>%
    pivot_longer(-c(DATE, SITECODE), names_to="Method", values_to="ET_mm_day")
}

# --- 2. ZHANG‐ONLY SCATTERPLOTS (2013–2019) ------------------------
for(site in site_list){
  d <- df %>%
    filter(SITECODE == site,
           between(DATE, as.Date("2013-01-01"), as.Date("2019-12-31")))
  
  p_month <- ggplot(d, aes(x = ET_PT_zhang)) +
    geom_point(aes(y = ET_Hamon_pt_zhang_monthly,
                   color = "Hamon (Monthly Calibrated, Zhang)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = zhang_colors["Hamon (Monthly Calibrated, Zhang)"]) +
    labs(title=paste(site, "- Zhang (monthly)"),
         x="PT–Zhang ET", y="Hamon ET") +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  p_interp <- ggplot(d, aes(x = ET_PT_zhang)) +
    geom_point(aes(y = ET_Hamon_pt_zhang_interp_full,
                   color = "Hamon (Interpolated, Zhang)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = zhang_colors["Hamon (Interpolated, Zhang)"]) +
    labs(title=paste(site, "- Zhang (interpolated)"),
         x="PT–Zhang ET", y=NULL) +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  p_uncal <- ggplot(d, aes(x = ET_PT_zhang)) +
    geom_point(aes(y = ET_Hamon_uncalibrated,
                   color = "Hamon Uncalibrated"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = zhang_colors["Hamon Uncalibrated"]) +
    labs(title=paste(site, "- Hamon Uncalibrated"),
         x="PT–Zhang ET", y="Hamon ET") +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  p_pt <- ggplot(d, aes(x = ET_PT_zhang)) +
    geom_point(aes(y = ET_PT_zhang,
                   color = "PT-Zhang (2024)"),
               size=1.2, alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    scale_color_manual(values = zhang_colors["PT-Zhang (2024)"]) +
    labs(title=paste(site, "- PT-Zhang"),
         x="PT–Zhang ET", y=NULL) +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
  
  combo <- (p_pt + p_uncal) / (p_month + p_interp) +
    plot_layout(guides="collect") &
    theme(legend.position="bottom", legend.title=element_blank())
  
  ggsave(file.path(zhang_scatter, paste0("zhang_scatter_4panel_", site, ".png")),
         combo, width=16, height=12)
}

# --- 3a. ZHANG‐ONLY PER-SITE TIME SERIES (omit PT-Zhang) ------------------
# define only the three Hamon variants
full_methods <- zhang_methods[names(zhang_methods) != "PT-Zhang (2024)"]
full_colors  <- zhang_colors[names(zhang_colors)    != "PT-Zhang (2024)"]

for(site in site_list){
  ts_long <- pivot_methods(df, full_methods) %>%
    filter(SITECODE == site)
  
  p_ts <- ggplot(ts_long, aes(DATE, ET_mm_day, color = Method)) +
    geom_line(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = full_colors) +
    labs(
      title = paste("Zhang‐only ET at", site),
      x     = "Date",
      y     = "ET (mm/day)",
      color = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position     = "bottom",
      axis.text.x         = element_text(angle = 45, hjust = 1),
      panel.grid.minor    = element_blank()
    )
  
  ggsave(
    file.path(zhang_dir, paste0("zhang_full_record_", site, ".png")),
    p_ts, width = 12, height = 6
  )
}

# --- 3b. ZHANG‐ONLY GRID ALL SITES, FULL RECORD (1997–2019, no PT-Zhang) ---
# pick only the three Hamon‐based methods
grid_methods <- zhang_methods[names(zhang_methods) != "PT-Zhang (2024)"]

all_zhang <- pivot_methods(df, grid_methods) %>%
  filter(between(DATE, as.Date("1997-01-01"), as.Date("2019-12-31")))

p_grid_zhang <- ggplot(all_zhang, aes(DATE, ET_mm_day, color=Method)) +
  geom_line(size=0.4, alpha=0.7) +
  facet_wrap(~SITECODE, scales="free_y", ncol=3) +
  # only pass the 3 remaining colors
  scale_color_manual(values = zhang_colors[names(zhang_colors) != "PT-Zhang (2024)"]) +
  labs(
    title="Zhang‐only ET Across All Sites (1997–2019)",
    x="Date", y="ET (mm/day)", color="Method"
  ) +
  theme_bw(base_size=12) +
  theme(
    strip.text       = element_text(face="bold", size=9),
    axis.text.x      = element_text(angle=45, hjust=1, size=7),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(zhang_dir, "zhang_grid_all_sites_1997_2019_no_PTZhang.png"),
  p_grid_zhang, width=16, height=10
)

# --- 4a. ZHANG‐ONLY CALIBRATION WINDOW PER-SITE (2013–2019) ---------
for(site in site_list){
  cal_long <- pivot_methods(df, zhang_methods) %>%
    filter(SITECODE == site,
           between(DATE, as.Date("2013-01-01"), as.Date("2019-12-31")))
  
  p_cal <- ggplot(cal_long, aes(DATE, ET_mm_day, color=Method)) +
    geom_line(size=0.6, alpha=0.8) +
    scale_color_manual(values=zhang_colors) +
    labs(title=paste("Zhang‐only Calibration Window at", site),
         x="Date", y="ET (mm/day)", color="Method") +
    theme_bw(base_size=13) +
    theme(legend.position="bottom",
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.minor=element_blank())
  
  ggsave(file.path(zhang_dir, paste0("zhang_cal_window_", site, ".png")),
         p_cal, width=12, height=6)
}

# --- 4b. ZHANG‐ONLY GRID CALIBRATION WINDOW ------------------------
p_grid_cal_zhang <- ggplot(
  all_zhang %>% filter(between(DATE, as.Date("2013-01-01"), as.Date("2019-12-31"))),
  aes(DATE, ET_mm_day, color=Method)
) +
  geom_line(size=0.4, alpha=0.7) +
  facet_wrap(~SITECODE, scales="free_y", ncol=3) +
  scale_color_manual(values=zhang_colors) +
  labs(title="Zhang‐only ET Across All Sites (2013–2019)",
       x="Date", y="ET (mm/day)", color="Method") +
  theme_bw(base_size=12) +
  theme(strip.text       = element_text(face="bold", size=9),
        axis.text.x      = element_text(angle=45,hjust=1,size=7),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(zhang_dir, "zhang_grid_all_sites_2013_2019.png"),
       p_grid_cal_zhang, width=16, height=10)
