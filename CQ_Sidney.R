# Create some preliminary CQ plots of CQ for HJA data

# Load necessary libraries
librarian::shelf(tidyverse, googledrive, purrr, readxl, supportR)

# Clear environment
rm(list = ls())

# Specify the main folder path
main_folder_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/SiSyn/"

# Generate the filename with the current date
file_name <- paste0("HJA_CQ_merged_master_", Sys.Date(), ".csv")

# Read in the CSV file, select specific columns, convert Date to a date format, and filter by specific variables
cq_data <- read.csv(file.path(main_folder_path, file_name)) %>%
  dplyr::select(Stream_Name, Date, Qcms, variable, value) %>%  # Adjust column names as needed
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::filter(!variable %in% c("pH")) %>%  # Remove variables not interested in (e.g., pH)
  dplyr::mutate(variable = dplyr::recode(variable,
                                         "dissolved org C" = "DOC",
                                         "specific conductivity" = "spc"))  # Replace with actual variable names



