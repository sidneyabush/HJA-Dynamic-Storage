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

# ----------------------------------------------------------------------------
# 2. Calculate Daily Differences and Identify Recession Days
# ----------------------------------------------------------------------------
all_data <- all_data_daily %>%
  dplyr::arrange(SITECODE, Date) %>%
  dplyr::group_by(SITECODE) %>%
  dplyr::mutate(
    dQ = MEAN_Q_AREA - lag(MEAN_Q_AREA),
    change_dQ = MEAN_Q_AREA / lag(MEAN_Q_AREA),
    dQ_dt = dQ / as.numeric(Date - lag(Date))) %>%
  dplyr::filter(!is.na(dQ_dt)) %>% # Remove NA values (first row)
  dplyr::filter(!change_dQ < 0.7) 

# Calculate the recession slope (-dQ/dt)
recession_data <- all_data %>%
  dplyr::filter(dQ < 0) %>%  # Keep only recession periods
  dplyr::mutate(recession_slope = -dQ_dt)  # Make it positive for the slope

# ----------------------------------------------------------------------------
# 3. Compute Aggregate Recession Slope per Stream
# ----------------------------------------------------------------------------
# For each stream, if there are at least 50 recession days, fit a linear model (recession_slope ~ Q)
# and extract the slope coefficient.
recession_slopes <- recession_data %>%
  # 1) drop any rows where slope or Q isn't finite
  filter(is.finite(recession_slope), is.finite(MEAN_Q_AREA)) %>%
  
  # 2) group by stream
  group_by(SITECODE) %>%
  
  # 3) keep only streams with ≥50 valid recession days
  #filter(n() >= 50) %>%
  
  # 4) summarize once per stream
  summarise(
    n_days = n(),
    slope  = coef(lm(recession_slope ~ MEAN_Q_AREA))[2], 
    .groups = "drop"
  )

# View the result (one row per Stream_ID)
print(recession_slopes)

# -------------------------------
# Plot: dQ/dT vs Q for Multiple Sites
# -------------------------------
# Define the sites of interest by providing their Stream_ID values.
# Replace with the actual identifiers from your dataset if needed.
# make sure your interest vector is correct
sites_of_interest <- c("GSLOOK", "GSWSMC", "GSWS01", "GSWS02", "GSWS03", "GSWS06", "GSWS07",
                       "GSWS08", "GSWS09", "GSWS10")

multi_site_data <- recession_data %>%
  filter(SITECODE %in% sites_of_interest)

# drop any rows where y or x is not finite
multi_site_data_clean <- multi_site_data %>%
  filter(
    is.finite(recession_slope),
    is.finite(MEAN_Q_AREA)
  )

# now compute slopes and annotation positions
slopes <- multi_site_data %>%
  filter(is.finite(recession_slope), is.finite(MEAN_Q_AREA)) %>%
  group_by(SITECODE) %>%
  summarise(
    slope     = round(coef(lm(recession_slope ~ MEAN_Q_AREA))[2], 2),
    Q_pos     = quantile(MEAN_Q_AREA,      0.9, na.rm = TRUE),
    dQ_dt_pos = quantile(recession_slope,  0.9, na.rm = TRUE),
    .groups   = "drop"
  )


# Step 3: Create the plot with slope annotation
p_dQ_dt <- ggplot(multi_site_data, aes(x = MEAN_Q_AREA, y = recession_slope)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_smooth(method = "lm", color = "skyblue", se = FALSE) + 
  facet_wrap(~ SITECODE, scales = "free",  ncol   = 3) +  # One panel per site
  # Add slope annotations from the slopes dataframe
  geom_text(
    data = slopes, 
    aes(x = Q_pos, y = dQ_dt_pos, label = paste("slope =", slope)),
    hjust = 1, vjust = 1, size = 3, color = "skyblue"
  ) +
  labs(
    title = NULL,
    x = "Discharge (Q)",
    y = "dQ/dT"
  ) +
  theme_classic()

# Display the plot
print(p_dQ_dt)

# Save it
ggsave(
  filename = "recession_limb_by_site.png",
  plot     = p_dQ_dt,
  path     = storage_dir,
  width    = 7,      # adjust as needed
  height   = 8,      # adjust as needed
  units    = "in",
  dpi      = 300
)


# ----------------------------------------------------------------------------
# 5. Export the Results
# ----------------------------------------------------------------------------
write.csv(recession_slopes, "Recession_Slopes_by_SITECODE_Aggregate.csv")


# Now do by water year: 
# Compute one recession‐slope per (Stream, Water Year)
recession_slopes_wy <- recession_data %>%
  filter(is.finite(recession_slope), is.finite(MEAN_Q_AREA)) %>%
  group_by(SITECODE, WATERYEAR) %>%
  summarise(
    n_days = n(),
    slope  = if (n_days >= 50) {
      coef(lm(recession_slope ~ MEAN_Q_AREA))[2]
    } else {
      NA_real_
    },
    .groups = "drop"
  )


recession_limb <- ggplot(recession_slopes_wy, aes(x = WATERYEAR, y = slope)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ SITECODE,
             ncol   = 3,        # 3 columns → ~4 rows for 10 sites
             scales = "fixed"  # each panel gets its own y-scale
  ) +
  labs(
    x     = "Water Year",
    y     = "Recession Slope (–dQ/dt vs Q)",
    title = "Annual Recession Slope by Site"
  ) +
  theme_classic()

# Save it
ggsave(
  filename = "recession_limb_by_site_wy.png",
  plot     = recession_limb,
  path     = storage_dir,
  width    = 7,      # adjust as needed
  height   = 8,      # adjust as needed
  units    = "in",
  dpi      = 300
)

