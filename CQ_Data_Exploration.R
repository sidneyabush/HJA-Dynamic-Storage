# Some data exploration: 

# Clear environment
rm(list = ls())

# Specify the main folder path
main_folder_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem"

# Generate the filename with the current date
file_name <- paste0("HJA_CQ_merged_slopes_", Sys.Date(), ".csv")

# Now more post-processing: 
cq_data <- read.csv(file.path(main_folder_path, file_name)) %>%
  dplyr::select(Stream_Name, Date, Qcms, variable, value, slope) %>%  # Adjust column names as needed
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::filter(!variable %in% c("pH", "spc"))  # Replace with additional variables you want to remove

# Create a table with the start and end date of each variable from each site
date_range_table <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarize(
    start_date = min(Date, na.rm = TRUE),  # Earliest date
    end_date = max(Date, na.rm = TRUE)     # Latest date
  ) %>%
  ungroup()

# Find the common date range across all sites and solutes
common_date_range <- date_range_table %>%
  summarize(
    common_start_date = max(start_date, na.rm = TRUE),  # Latest of the earliest dates
    common_end_date = min(end_date, na.rm = TRUE)       # Earliest of the latest dates
  )

# Find the limiting factors (latest start date and earliest end date) along with the variables
limiting_factors <- date_range_table %>%
  filter(
    start_date == max(start_date, na.rm = TRUE) |
      end_date == min(end_date, na.rm = TRUE)
  ) %>%
  mutate(limiting_type = case_when(
    start_date == max(start_date, na.rm = TRUE) ~ "Limiting Start",
    end_date == min(end_date, na.rm = TRUE) ~ "Limiting End"
  ))

# Find the top 5 limiting variables by latest start date and earliest end date
top_limiting_start <- date_range_table %>%
  arrange(desc(start_date)) %>%
  head(5)

top_limiting_end <- date_range_table %>%
  arrange(end_date) %>%
  head(5)

# Print the results
print(common_date_range)
print("Limiting Factors:")
print(limiting_factors)

print("Top 5 Limiting Variables by Start Date:")
print(top_limiting_start)

print("Top 5 Limiting Variables by End Date:")
print(top_limiting_end)

# Create two additional .csv files: 
# One with data filtered by common date range that has all sites WITHIN the Lookout Creek Watershed (1,2,3,6,7,8,LOOK,MACK)
lookout_creek_sites <- c("GSWS01", "GSWS02", "GSWS03", "GSWS06", "GSWS07", "GSWS08", "GSLOOK", "GSMACK")

# Filter the date range table to include only Lookout Creek sites
lookout_creek_date_range_table <- date_range_table %>%
  filter(Stream_Name %in% lookout_creek_sites)

# Find the common start and end dates where data are available for all variables across all Lookout Creek sites
lookout_creek_common_date_range <- lookout_creek_date_range_table %>%
  group_by(variable) %>%
  summarize(
    variable_common_start_date = max(start_date, na.rm = TRUE),  # Latest start date for each variable
    variable_common_end_date = min(end_date, na.rm = TRUE)       # Earliest end date for each variable
  ) %>%
  summarize(
    common_start_date = max(variable_common_start_date, na.rm = TRUE),  # Latest of the latest start dates
    common_end_date = min(variable_common_end_date, na.rm = TRUE)       # Earliest of the earliest end dates
  )

# Filter the dataset by the common date range specific to Lookout Creek sites
cq_data_lookout_filtered <- cq_data %>%
  # Filter by the common date range specific to Lookout Creek sites
  filter(
    Date >= lookout_creek_common_date_range$common_start_date,
    Date <= lookout_creek_common_date_range$common_end_date
  ) %>%
  # Filter to include only sites within the Lookout Creek Watershed
  filter(Stream_Name %in% lookout_creek_sites)

write.csv(cq_data_lookout_filtered, file = file.path(main_folder_path, "cq_data_lookout_filtered.csv"), row.names = FALSE)


# Now perform the same operations for all sites
# Calculate the common date range for all variables across all sites
common_date_range_all_sites <- date_range_table %>%
  group_by(variable) %>%
  summarize(
    variable_common_start_date = max(start_date, na.rm = TRUE),  # Latest start date for each variable
    variable_common_end_date = min(end_date, na.rm = TRUE)       # Earliest end date for each variable
  ) %>%
  summarize(
    common_start_date = max(variable_common_start_date, na.rm = TRUE),  # Latest of the latest start dates
    common_end_date = min(variable_common_end_date, na.rm = TRUE)       # Earliest of the earliest end dates
  )

# Filter the dataset by the common date range for all sites
cq_data_all_filtered <- cq_data %>%
  # Filter by the common date range for all sites
  filter(
    Date >= common_date_range_all_sites$common_start_date,
    Date <= common_date_range_all_sites$common_end_date
  )

write.csv(cq_data_all_filtered, file = file.path(main_folder_path, "cq_data_all_filtered.csv"), row.names = FALSE)

