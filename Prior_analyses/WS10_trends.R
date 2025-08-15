# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear environment
rm(list = ls())

setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/HJA_storm_sampling2020/SBush_data_exploration/Final_Dataset_20240313")

# # Open a single PDF file for all plots (this ensures they are saved as separate pages)
# pdf_path <- pdf("/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/WS10_Prelim_Plots.pdf", width = 8, height = 10)
# 
# # Open the PDF file
# pdf(pdf_path, width = 8, height = 10)

## -----------------------------------------
#### CHEMISTRY - PROPORTIONAL ####
## -----------------------------------------

# Load proportional stream chemistry data (Oct 1968 - May 2019)
prop_chem_1 <- read.csv("./PUBLISHED/CF00201_v6.csv", header = TRUE, na.strings = c("n.a."), stringsAsFactors = FALSE)
prop_chem_1$DATE_TIME <- as.POSIXct(prop_chem_1$DATE_TIME, format="%Y-%m-%d %H:%M")

# Filter chemistry data for WS10
prop_chem_1 <- prop_chem_1 %>%
  filter(SITECODE == "GSWS10")

# Extract relevant solutes (including DOC)
solutes <- c("NO3N", "NA.", "K", "CA", "MG", "CL", "DOC")
chem_subset <- prop_chem_1 %>%
  select(DATE_TIME, SITECODE, all_of(solutes)) %>%
  rename(DATE = DATE_TIME)  # Rename for merging

## DAILY PUBLISHED DISCHARGE DATA ----
meanDailyQ_pub <- read.csv("./PUBLISHED/HF00402_v14_daily.csv", header = TRUE, na.strings = c("n.a."), stringsAsFactors = FALSE)
meanDailyQ_pub$DATE <- as.Date(meanDailyQ_pub$DATE, format="%Y-%m-%d")

# Filter for WS10 and select required columns
valid_sites <- c("GSWS10")
Q_cols <- c("SITECODE", "DATE", "MEAN_Q")

meanDailyQ_pub <- meanDailyQ_pub %>%
  filter(SITECODE %in% valid_sites) %>%
  select(all_of(Q_cols)) %>%
  mutate(SITECODE = as.factor(SITECODE))

# Merge discharge and chemistry data by date
combined_data <- left_join(chem_subset, meanDailyQ_pub, by = c("DATE", "SITECODE"))

# Convert to long format for faceting
long_data <- combined_data %>%
  pivot_longer(cols = all_of(solutes), names_to = "Solute", values_to = "Concentration")

# Ensure DATE column is in POSIXct format
long_data$DATE <- as.POSIXct(long_data$DATE)
meanDailyQ_pub$DATE <- as.POSIXct(meanDailyQ_pub$DATE)

# Create discharge data separately
discharge_data <- meanDailyQ_pub %>%
  mutate(Solute = "Discharge", Concentration = MEAN_Q)

# Combine discharge and solute data
plot_data <- bind_rows(discharge_data, long_data)

# Explicitly set "Discharge" to be the first panel manually
plot_data$Solute <- factor(plot_data$Solute, levels = c("Discharge", "DOC", "CA", "CL", "K", "MG", "NA.", "NO3N", "PO4P", "SO4S"))

# # Define the correct PDF path
# pdf_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/WS10_Prelim_Plots.pdf"
# 
# # Open a single PDF file for all plots (this ensures they are saved as separate pages)
# pdf(pdf_path, width = 8, height = 11)

## -----------------------------------------
#### Update DOC Dates to Show Only 3 Most Recent Years ####
## -----------------------------------------

# Get the most recent year in the dataset
latest_year <- max(format(long_data$DATE, "%Y"), na.rm = TRUE)

# Get the most recent year in the dataset (convert to numeric)
latest_year <- as.numeric(max(format(long_data$DATE, "%Y"), na.rm = TRUE))

# Filter for the last 3 years
doc_dates <- long_data %>%
  filter(Solute == "DOC" & !is.na(Concentration) & 
           as.numeric(format(DATE, "%Y")) >= (latest_year - 2)) %>%
  summarise(start_date = min(DATE), end_date = max(DATE))


# Filter solutes and discharge data within the updated DOC date range
solutes_data <- long_data %>%
  filter(DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Solutes")

discharge_data <- meanDailyQ_pub %>%
  filter(DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Solute = "Discharge", Concentration = MEAN_Q, Panel = "Discharge")

# Combine datasets
plot_data <- bind_rows(discharge_data, solutes_data)

# Ensure DATE is in correct format
plot_data$DATE <- as.POSIXct(plot_data$DATE)
solutes_data$DATE <- as.POSIXct(solutes_data$DATE)
discharge_data$DATE <- as.POSIXct(discharge_data$DATE)

# ## ------------------------------
# #### PLOT 1: Solutes with Discharge as Gray Line (Discharge at the Top) ####
# ## ------------------------------
# print(
#   ggplot(plot_data, aes(x = DATE, y = Concentration)) +
#     geom_line(data = discharge_data, aes(y = Concentration), color = "gray", size = 1) +
#     geom_line(data = long_data, aes(y = Concentration, color = Solute), size = 1) +
#     facet_grid(rows = vars(Solute), scales = "free_y", switch = "y") +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +  # Remove gridlines
#     labs(title = "WS10",
#          x = "Date", y = NULL, color = "Solutes") +
#     theme(legend.position = "bottom")
# )

## ------------------------------
#### PLOT 2: Discharge in Its Own Panel Above Solutes ####
## ------------------------------
print(
  ggplot(plot_data, aes(x = DATE, y = Concentration)) +
    geom_line(data = discharge_data, aes(y = Concentration), color = "darkgray", size = 1) +
    geom_line(data = solutes_data, aes(y = Concentration, color = Solute), size = 1) +
    facet_grid(rows = vars(Panel), scales = "free_y", switch = "y") +
    scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") +  # Add more x-axis tick marks
    theme_bw() +
    theme(panel.grid = element_blank()) +  # Remove gridlines
    labs(title = "WS10",
         x = "Date", y = NULL, color = "Solutes") +
    theme(legend.position = "bottom")
)

## ------------------------------
#### PLOT 3: Revised Solutes (Updated DOC Date Range) ####
## ------------------------------
# Define the new solute list for Plot 3
solutes_revised <- c("DOC", "CA", "NA.", "MG", "CL")  # Adjust as needed

# Filter solutes within the updated DOC date range
solutes_data_revised <- long_data %>%
  filter(Solute %in% solutes_revised & DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Solutes")

# Filter discharge data within the updated DOC date range
discharge_data_revised <- discharge_data %>%
  filter(DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Discharge")

# Combine datasets
plot_data_revised <- bind_rows(discharge_data_revised, solutes_data_revised)

# Ensure "Discharge" is at the top
plot_data_revised$Panel <- factor(plot_data_revised$Panel, levels = c("Discharge", "Solutes"))

print(
  ggplot(plot_data_revised, aes(x = DATE, y = Concentration)) +
    geom_line(data = discharge_data_revised, aes(y = Concentration), color = "gray", size = 1) +
    geom_line(data = solutes_data_revised, aes(y = Concentration, color = Solute), size = 1) +
    scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") +
  facet_grid(rows = vars(Panel), scales = "free_y", switch = "y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +  # Remove gridlines
    labs(title = "WS10",
         x = "Date", y = NULL, color = "Solutes") +
    theme(legend.position = "bottom")
)

## ------------------------------
#### PLOT 4: Revised Solutes (Updated DOC Date Range) ####
## ------------------------------
# Define the new solute list for Plot 3
solutes_revised <- c("DOC", "CA", "NA.", "MG")  # Adjust as needed

# Filter solutes within the updated DOC date range
solutes_data_revised <- long_data %>%
  filter(Solute %in% solutes_revised & DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Solutes")

# Filter discharge data within the updated DOC date range
discharge_data_revised <- discharge_data %>%
  filter(DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Discharge")

# Combine datasets
plot_data_revised <- bind_rows(discharge_data_revised, solutes_data_revised)

# Ensure "Discharge" is at the top
plot_data_revised$Panel <- factor(plot_data_revised$Panel, levels = c("Discharge", "Solutes"))

print(
  ggplot(plot_data_revised, aes(x = DATE, y = Concentration)) +
    geom_line(data = discharge_data_revised, aes(y = Concentration), color = "gray", size = 1) +
    geom_line(data = solutes_data_revised, aes(y = Concentration, color = Solute), size = 1) +
    scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") +
    facet_grid(rows = vars(Panel), scales = "free_y", switch = "y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +  # Remove gridlines
    labs(title = "WS10",
         x = "Date", y = NULL, color = "Solutes") +
    theme(legend.position = "bottom")
)

## ------------------------------
#### PLOT 5: Revised Solutes (Updated DOC Date Range) ####
## ------------------------------
# Define the new solute list for Plot 3
solutes_revised <- c("DOC", "MG")  # Adjust as needed

# Filter solutes within the updated DOC date range
solutes_data_revised <- long_data %>%
  filter(Solute %in% solutes_revised & DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Solutes")

# Filter discharge data within the updated DOC date range
discharge_data_revised <- discharge_data %>%
  filter(DATE >= doc_dates$start_date & DATE <= doc_dates$end_date) %>%
  mutate(Panel = "Discharge")

# Combine datasets
plot_data_revised <- bind_rows(discharge_data_revised, solutes_data_revised)

# Ensure "Discharge" is at the top
plot_data_revised$Panel <- factor(plot_data_revised$Panel, levels = c("Discharge", "Solutes"))

library(ggplot2)

# Define colorblind-friendly palette for specific solutes
colorblind_palette <- c(
  "DOC" = "#009E73",
  "MG" = "#E69F00")

# Ensure Solute column is a factor
solutes_data_revised$Solute <- as.factor(solutes_data_revised$Solute)

# Rename facet labels
facet_labels <- c("Solutes" = "Concentration (mg/L)", "Discharge" = "Discharge (cfs)")

# Define the plot
p <- ggplot(plot_data_revised, aes(x = DATE, y = Concentration)) +
  geom_line(data = discharge_data_revised, aes(y = Concentration), color = "black", size = 1) +
  geom_line(data = solutes_data_revised, aes(y = Concentration, color = Solute), size = 1) +
  scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") +
  scale_color_manual(values = colorblind_palette) +  # Apply colorblind-friendly colors
  facet_grid(rows = vars(Panel), scales = "free_y", switch = "y", labeller = labeller(Panel = facet_labels)) +
  theme_bw() +
  labs(title = "WS10",
       x = "Date", y = NULL, color = "Solutes") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),  # Increase title size
    legend.title = element_blank(),  
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.title.x = element_text(vjust = -0.5),  # Move x-axis label down using vjust
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.position = c(0.9, 0.4),  # Place legend inside the second panel
    strip.text = element_text(size = 18, face = "bold")  # Increase facet title text size
  )

# Display the plot
print(p)

# Save the plot as a TIFF file
output_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/WS10_plot.tiff"
ggsave(output_path, plot = p, width = 10, height = 8, dpi = 300, units = "in", device = "tiff")
