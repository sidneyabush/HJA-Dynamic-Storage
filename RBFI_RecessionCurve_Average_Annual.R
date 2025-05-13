library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(scales)
library(colorspace)

# ─── OPTION 1: bigger default text ───────────────────────────────────────────
theme_set(theme_classic(base_size = 14))

# Clear environment
rm(list = ls())

# Directories
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Only these 8 sites
sites_keep <- c("GSLOOK","GSWS01","GSWS02","GSWS03",
                "GSWS06","GSWS07","GSWS08","GSWSMC")

# Read and prep data
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

# Function: recession slope + its SE
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
  cs  <- summary(fit)$coefficients
  tibble(
    slope    = cs[2,1],
    se_slope = cs[2,2]
  )
}

# Function: RBFI + its SE (based on daily abs(dQ)/total_Q)
calc_rbfi <- function(df) {
  tmp <- df %>%
    mutate(dQ = Q - lag(Q)) %>%
    filter(!is.na(dQ))
  total_Q <- sum(df$Q, na.rm = TRUE)
  props   <- abs(tmp$dQ) / total_Q
  tibble(
    RBFI    = sum(abs(tmp$dQ), na.rm = TRUE) / total_Q,
    se_RBFI = sd(props, na.rm = TRUE) / sqrt(nrow(tmp))
  )
}

# Compute annual metrics per site × water year
annual_metrics <- discharge %>%
  group_by(SITECODE, WATERYEAR) %>%
  group_map(~ bind_cols(
    tibble(site = .y$SITECODE, year = .y$WATERYEAR),
    calc_recession(.x),
    calc_rbfi(.x)
  )) %>%
  bind_rows() %>%
  mutate(site = factor(site, levels = sites_keep))

# Prepare full‐record recession data for curves
recession_clean <- discharge %>%
  mutate(
    dQ           = Q - lag(Q),
    dQ_dt        = dQ / as.numeric(Date - lag(Date)),
    change_ratio = Q / lag(Q)
  ) %>%
  filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
  transmute(
    site            = factor(SITECODE, levels = sites_keep),
    Q,
    recession_slope = -dQ_dt
  )

# Annotation: slopes & positions
slopes_df <- recession_clean %>%
  group_by(site) %>%
  summarise(
    slope     = round(coef(lm(recession_slope ~ Q))[2], 2),
    Q_pos     = quantile(Q,              0.9, na.rm = TRUE),
    dQ_dt_pos = quantile(recession_slope, 0.85, na.rm = TRUE),
    .groups   = "drop"
  )

# Color palette → lighten by 10%
orig_palette8 <- c(
  "#B276B2", "#E87272", "#EF7832", "#EEC64B",
  "#88B372", "#5FBDA9", "#87BCE3", "#443A73"
)
site_cols <- setNames(lighten(orig_palette8, amount = 0.1), sites_keep)

# 1) Annual RBFI ± SE by site (unchanged)
p_rbfi <- ggplot(annual_metrics, aes(x = year, y = RBFI, color = site, fill = site)) +
  geom_ribbon(aes(ymin = RBFI - se_RBFI, ymax = RBFI + se_RBFI),
              alpha = 0.3, color = NA) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  facet_wrap(~ site, ncol = 4) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  scale_fill_manual(values = alpha(site_cols, 0.3), guide = FALSE) +
  labs(x = "Water Year", y = "Richards Baker Index") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )
ggsave("rbfi_by_site_wy_se.png", p_rbfi, path = output_dir,
       width = 12, height = 6, units = "in", dpi = 300)

# 2) Annual recession‐curve slope ± SE by site: label = “β (d⁻¹)”
p_slope <- ggplot(annual_metrics, aes(x = year, y = slope, color = site, fill = site)) +
  geom_ribbon(aes(ymin = slope - se_slope, ymax = slope + se_slope),
              alpha = 0.3, color = NA) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  facet_wrap(~ site, ncol = 4) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  scale_fill_manual(values = alpha(site_cols, 0.3), guide = FALSE) +
  labs(
    x = "Water Year",
    y = "Recession Limb Slope (d⁻¹)"
  ) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )
ggsave("recession_slope_by_site_wy_se.png", p_slope, path = output_dir,
       width = 12, height = 6, units = "in", dpi = 300)

# 3) Instantaneous recession‐limb curves: label = “–dQ/dt (m³ s⁻¹ d⁻¹)”
p_curve <- ggplot(recession_clean, aes(x = Q, y = recession_slope, color = site)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  geom_text(data = slopes_df,
            aes(label = paste0("slope = ", slope)),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6,
            size = 4, color = "black") +
  facet_wrap(~ site, ncol = 4, scales = "free") +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(
    x = expression(Discharge~(m^3/s)),
    y = "–dQ/dt (m³ s⁻¹ d⁻¹)"
  ) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black")
  )
ggsave("recession_curve_by_site.png", p_curve, path = output_dir,
       width = 12, height = 6, units = "in", dpi = 300)
