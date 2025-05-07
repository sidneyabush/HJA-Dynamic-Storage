# Load needed libraries
librarian::shelf(dplyr, ggplot2, data.table, ggpubr, reshape2, EflowStats, zoo, lubridate, readr)

# Clear environment
rm(list = ls())

# 1. Set your working directory to the data folder
setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/HJA_Hydrology_updatesthrough WY2024/Data_Cleaning_SBush/data")

# Set up a Figures folder one level above your current WD
fig_dir <- file.path(dirname(getwd()), "Figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# Now create a “Storage” subfolder inside Figures
storage_dir <- file.path(fig_dir, "Storage")
if (!dir.exists(storage_dir)) dir.create(storage_dir)


# ----------------------------------------------------------------------------
# 1. Read and Prepare the Data; Remove Duplicate Date Entries per Stream_ID
# ----------------------------------------------------------------------------
all_data_3hourly <- read_csv(
  "HJA_GSLOOK_3hr_WY_1995_2024.csv",
  show_col_types = FALSE
) %>%
  mutate(
    DATE_TIME = ymd_hms(DATE_TIME, tz = "UTC"),
    Date      = as_date(DATE_TIME)
  )

all_data_daily <- all_data_3hourly %>%
  group_by(SITECODE, WATERYEAR, Date) %>%    # <-- group by site & date
  summarise(
    across(                                  # mean of all your numeric vars
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


# -----------------------------------------------------------
# 2. Calculate Flashiness (RBFI) for Each SITECODE ----
# -----------------------------------------------------------
# For each stream, calculate daily discharge changes and compute RBFI.
flashiness <- all_data_daily %>%
  group_by(SITECODE) %>%
  arrange(Date) %>%                     # Ensure dates are in order for each stream
  mutate(dQ = MEAN_Q_AREA - lag(MEAN_Q_AREA),                # Daily change in discharge
         abs_dQ = abs(dQ)) %>%           # Absolute change in discharge
  filter(!is.na(abs_dQ)) %>%             # Remove NA from the first row (due to lag)
  summarise(
    total_discharge = sum(MEAN_Q_AREA, na.rm = TRUE),         # Total discharge over the period
    total_change = sum(abs_dQ, na.rm = TRUE),         # Total absolute change
    RBFI = total_change / total_discharge           # Richards-Baker Flashiness Index
  ) %>%
  ungroup()

# View the flashiness data frame with RBFI values for each SITECODE
print(flashiness)

# Keep only the SITECODE and RBFI columns
flashiness_export <- flashiness %>%
  dplyr::select(SITECODE, RBFI)

# Export the result as a CSV file
write_csv(flashiness_export, "flashiness_by_SITECODE.csv")

# -----------------------------------------------------------
# 2) Compute RBFI by SITECODE × WATERYEAR
# -----------------------------------------------------------
flashiness_wy <- all_data_daily %>%
  group_by(SITECODE, WATERYEAR) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    dQ     = MEAN_Q_AREA - lag(MEAN_Q_AREA),
    abs_dQ = abs(dQ)
  ) %>%
  filter(!is.na(abs_dQ)) %>%
  summarise(
    total_discharge = sum(MEAN_Q_AREA, na.rm = TRUE),
    total_change    = sum(abs_dQ,   na.rm = TRUE),
    RBFI            = total_change / total_discharge,
    .groups = "drop"
  )

# peek
print(flashiness_wy)

# -----------------------------------------------------------
# 3) Plot with facet_grid
# -----------------------------------------------------------
# Build the plot object
p_flashiness <- ggplot(flashiness_wy, aes(x = WATERYEAR, y = RBFI)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ SITECODE,
             ncol   = 3,        # 3 columns → ~4 rows for 10 sites
             scales = "fixed"  # each panel gets its own y-scale
  ) +  
  labs(
    x = "Water Year",
    y = "Richards–Baker Flashiness Index",
    title = "Annual Flashiness (RBFI) by Site"
  ) +
  theme_classic() +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.spacing = unit(0.3, "lines")
  )

# Save it
ggsave(
  filename = "flashiness_by_site_wy.png",
  plot     = p_flashiness,
  path     = storage_dir,
  width    = 7,      # adjust as needed
  height   = 8,      # adjust as needed
  units    = "in",
  dpi      = 300
)

