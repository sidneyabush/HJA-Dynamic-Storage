## Preliminary CQ Analysis for Long-Term HJA Dataset
## Adapted from Nick's code for SiSyn CQ project

# Load necessary libraries
librarian::shelf(tidyverse, googledrive, purrr, readxl, supportR)

# Clear environment
rm(list = ls())

# Google Drive folder URLs
folder_ids <- c(
  "0AIPkWhVuXjqFUk9PVA",  # main folder
  "1dTENIB5W2ClgW0z-8NbjqARiaGO2_A7W",  # chemistry folder
  "1hbkUsTdo4WAEUnlPReOUuXdeeXm92mg-"  # discharge folder
)

# Function to fetch and filter the most recent files
fetch_latest_files <- function(folder_id) {
  googledrive::drive_ls(googledrive::as_id(folder_id)) %>%
    dplyr::filter(
      stringr::str_detect(name, "Site_Reference_Table") |
        stringr::str_detect(name, "masterdata_discharge") |
        stringr::str_detect(name, "masterdata_chem")
    ) %>%
    dplyr::arrange(desc(name)) %>%
    dplyr::slice(1)
}

# Fetch files from all directories
wanted_files <- purrr::map_dfr(folder_ids, fetch_latest_files)

# Specify the directory path for downloads
dir_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/SiSyn/raw_data"
dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

# Download the files
purrr::walk2(wanted_files$name, wanted_files$id, ~ googledrive::drive_download(.y, path = file.path(dir_path, .x), overwrite = TRUE))

# Read and process data
ref_table <- readxl::read_excel(file.path(dir_path, "Site_Reference_Table.xlsx")) %>%
  dplyr::select(LTER, Stream_Name, Discharge_File_Name) %>%
  dplyr::distinct()

chem_v1 <- read.csv(file.path(dir_path, wanted_files$name[stringr::str_detect(wanted_files$name, "masterdata_chem")])) %>%
  dplyr::rename(Date = date) %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::filter(!is.na(Date) & !is.na(variable) & !is.na(value))

disc_v1 <- read.csv(file.path(dir_path, wanted_files$name[stringr::str_detect(wanted_files$name, "masterdata_discharge")])) %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::filter(!is.na(Date) & !is.na(Qcms))

# Filter data by reference table
disc_v3 <- disc_v1 %>%
  dplyr::filter(Discharge_File_Name %in% ref_table$Discharge_File_Name)

chem_v3 <- chem_v1 %>%
  dplyr::filter(Stream_Name %in% ref_table$Stream_Name)

# Merge datasets
cq_v1 <- disc_v3 %>%
  dplyr::full_join(chem_v3, by = c("Stream_Name", "Date"))

cq_v2 <- cq_v1 %>%
  dplyr::filter(!is.na(Qcms) & !is.na(variable) & !is.na(value)) %>%
  dplyr::distinct() %>%
  dplyr::select(LTER = LTER.x, Stream_Name, Date, Qcms, Dataset, variable, units, value)

# Filter by LTER - only want HJ Andrews sites
cq_v2_filtered <- cq_v2 %>%
  dplyr::filter(LTER == "AND")

# Generate the filename with the current date
file_name <- paste0("HJA_CQ_merged_master_", Sys.Date(), ".csv")

# Save the filtered data with the dynamic filename to the main folder
main_folder_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/SiSyn/"
write.csv(cq_v2_filtered, file = file.path(main_folder_path, file_name))


# End ----
