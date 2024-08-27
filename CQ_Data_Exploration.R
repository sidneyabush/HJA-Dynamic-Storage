# Some data exploration: Can run after running CQ_wrangling.R/ CQ_Sidney.R

# Load necessary libraries
librarian::shelf(tidyverse, dplyr, ggplot2, tidyr)

# Clear environment
rm(list = ls())

# Specify the main folder path
main_folder_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/SiSyn/"

# Generate the filename with the current date
file_name <- paste0("HJA_CQ_merged_master_", Sys.Date(), ".csv")

# Step 1: Read in and clean cq data ----
cq_data <- read.csv(file.path(main_folder_path, file_name)) %>%
  dplyr::select(Stream_Name, Date, Qcms, variable, value) %>%  # Adjust column names as needed
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::filter(!variable %in% c("pH")) %>%  # Remove variables not interested in (e.g., pH)
  dplyr::mutate(variable = dplyr::recode(variable,
                                         "dissolved org C" = "DOC",
                                         "specific conductivity" = "spc")) %>%  # Replace with actual variable names
  tidyr::drop_na()  # Remove rows with any NA values

# Verify the structure of cq_data
str(cq_data)

# Step 2: Manually calculate slopes for each combination of Stream_Name and variable ----

# Initialize an empty data frame to store the results
slope_results <- data.frame(Stream_Name = character(),
                            variable = character(),
                            slope = numeric(),
                            stringsAsFactors = FALSE)

# Get unique combinations of Stream_Name and variable
unique_combinations <- unique(cq_data %>% dplyr::select(Stream_Name, variable))

# Iterate over each combination and calculate the slope
for (i in 1:nrow(unique_combinations)) {
  subset_data <- cq_data %>%
    filter(Stream_Name == unique_combinations$Stream_Name[i],
           variable == unique_combinations$variable[i],
           Qcms > 0,
           value > 0)
  
  if (nrow(subset_data) > 1) {
    # Fit the linear model on log-transformed data
    lm_result <- tryCatch({
      lm(log10(value) ~ log10(Qcms), data = subset_data)
    }, error = function(e) NULL)
    
    # Extract the slope if the model fit was successful
    if (!is.null(lm_result)) {
      slope <- coef(lm_result)[2]
    } else {
      slope <- NA
    }
  } else {
    slope <- NA
  }
  
  # Append the result to the data frame
  slope_results <- rbind(slope_results, data.frame(
    Stream_Name = unique_combinations$Stream_Name[i],
    variable = unique_combinations$variable[i],
    slope = slope,
    stringsAsFactors = FALSE
  ))
}

# Verify the slope results
str(slope_results)

# Join the slope data back to the original dataset
cq_data <- cq_data %>%
  left_join(slope_results, by = c("Stream_Name", "variable"))

# Step 3: Calculate the average slope for each Stream_Name and variable ----
average_slope_data <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarize(average_slope = mean(slope, na.rm = TRUE)) %>%
  ungroup()

# Verify the structure of average_slope_data
str(average_slope_data)

# Join the average slope data back to the original dataset
cq_data <- cq_data %>%
  left_join(average_slope_data, by = c("Stream_Name", "variable"))

# Step 4: Categorize the average slope and plot CQ for each Stream_Name and variable ----
# CQ for All data ----
# Calculate p-values and classify significance
cq_data <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  mutate(
    p_value = tryCatch({
      if (n() > 1) {
        summary(lm(log10(value) ~ log10(Qcms)))$coefficients[2, 4]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    significant = case_when(
      p_value < 0.05 & average_slope > 0.1 ~ "Mobilizing Significant",
      p_value < 0.05 & average_slope < -0.1 ~ "Diluting Significant",
      p_value < 0.05 & average_slope <= 0.1 & average_slope >= -0.1 ~ "Chemostatic Significant",
      TRUE ~ "Not Significant"
    )
  ) %>%
  ungroup()

# Create a table with the start and end date of each variable from each site
date_range_table <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarize(
    start_date = min(Date, na.rm = TRUE),  # Find the earliest date
    end_date = max(Date, na.rm = TRUE)     # Find the latest date
  ) %>%
  ungroup()  # Ungroup to avoid any unintended grouping in subsequent operations

# Find the common start time and end time across all sites and solutes
common_date_range <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarize(
    start_date = min(Date, na.rm = TRUE),  # Find the earliest date for each site and solute
    end_date = max(Date, na.rm = TRUE)     # Find the latest date for each site and solute
  ) %>%
  ungroup() %>%
  summarize(
    common_start_date = max(start_date, na.rm = TRUE),  # Find the latest of the earliest dates
    common_end_date = min(end_date, na.rm = TRUE)       # Find the earliest of the latest dates
  )

# Print the common date range
print(common_date_range)

# Find the start and end dates for each site and solute
date_range_table <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarize(
    start_date = min(Date, na.rm = TRUE),  # Find the earliest date
    end_date = max(Date, na.rm = TRUE)     # Find the latest date
  ) %>%
  ungroup()

# Print the table
print(date_range_table)

# Find the latest start date (limiting start) and earliest end date (limiting end)
limiting_start <- date_range_table %>%
  filter(start_date == max(start_date, na.rm = TRUE))

limiting_end <- date_range_table %>%
  filter(end_date == min(end_date, na.rm = TRUE))

# Print the limiting factors
print("Limiting Start Date Factor:")
print(limiting_start)

print("Limiting End Date Factor:")
print(limiting_end)

