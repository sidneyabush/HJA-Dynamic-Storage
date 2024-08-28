# Load necessary libraries
librarian::shelf(tidyverse, dplyr, ggplot2, tidyr, lubridate)

# Clear environment
rm(list = ls())

# Set the working directory
setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/")

# Step 1: Read in the dataset for all sites ----
cq_all <- read.csv("cq_data_all_filtered.csv") %>%
  filter(!variable %in% c("TN", "TP"))  # Remove "TN" and "TP" variables

# Define the desired order of sites for all data
all_sites_order <- c("GSWS08", "GSWS07", "GSWS06", "GSMACK", "GSWS03", "GSWS02", "GSWS01", "GSWS09", "GSWS10", "GSLOOK")

# Define seasons based on months
assign_season <- function(month) {
  case_when(
    month %in% c(7, 8) ~ "Summer",        # July, August
    month %in% c(9, 10, 11, 12) ~ "Fall",  # September to December
    month %in% c(1, 2, 3) ~ "Winter",      # January to March
    month %in% c(4, 5, 6) ~ "Spring"       # April to June
  )
}

# Step 2: Add season and process data ----
process_cq_data <- function(cq_data, site_order) {
  # Ensure the Stream_Name is a factor with the specified order
  cq_data <- cq_data %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))
  
  # Add a Season column based on the month
  cq_data <- cq_data %>%
    mutate(
      Date = as.Date(Date),
      Month = month(Date),
      Season = assign_season(Month)
    )
  
  # Calculate the average slope for each Stream_Name, variable, and season
  average_slope_data <- cq_data %>%
    group_by(Stream_Name, variable, Season) %>%
    summarize(average_slope = mean(slope, na.rm = TRUE)) %>%
    ungroup()
  
  # Join the average slope data back to the original dataset
  cq_data <- cq_data %>%
    left_join(average_slope_data, by = c("Stream_Name", "variable", "Season"))
  
  return(cq_data)
}

# Step 3: Calculate p-values, R² values, and classify relationships ----
calculate_statistics <- function(cq_data) {
  stats <- cq_data %>%
    group_by(Stream_Name, variable, Season) %>%
    summarize(
      p_value = tryCatch({
        if (n() > 1) {
          summary(lm(log10(value) ~ log10(Qcms)))$coefficients[2, 4]
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_),
      r_squared = tryCatch({
        if (n() > 1) {
          summary(lm(log10(value) ~ log10(Qcms)))$r.squared
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_),
      slope = tryCatch({
        if (n() > 1) {
          coef(lm(log10(value) ~ log10(Qcms)))[2]
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_
      )
    ) %>%
    mutate(
      significance = case_when(
        p_value < 0.05 & slope > 0.1 ~ "Mobilizing Significant",
        p_value < 0.05 & slope < -0.1 ~ "Diluting Significant",
        p_value < 0.05 & slope <= 0.1 & slope >= -0.1 ~ "Chemostatic Significant",
        p_value >= 0.05 ~ "Not Significant",
        TRUE ~ "NA"
      )
    ) %>%
    ungroup()
  
  return(stats)
}

# Step 4: Plotting function without lines ----
plot_cq_data <- function(cq_data, filename, dataset_type) {
  # Get the start and end dates for the subtitle
  start_date <- min(cq_data$Date, na.rm = TRUE)
  end_date <- max(cq_data$Date, na.rm = TRUE)
  
  # Create the plot (points only, no lines)
  cq_plot <- ggplot(cq_data, aes(x = Qcms, y = value, color = Season)) +
    geom_point(alpha = 0.6, size = 1) +  # Points colored by season
    scale_x_log10() +  # Logarithmic scale for discharge
    scale_y_log10() +  # Logarithmic scale for concentration
    facet_grid(Stream_Name ~ variable, scales = "free") +  # Facet by Stream_Name and variable
    scale_color_manual(
      values = c(
        "Summer" = "#E69F00",    # Orange for Summer
        "Fall" = "#009E73",      # Green for Fall
        "Winter" = "#56B4E9",    # Blue for Winter
        "Spring" = "#F0E442"     # Yellow for Spring
      )
    ) +
    labs(
      title = paste("C-Q Relationships for", dataset_type),
      subtitle = paste("Date Range:", start_date, "to", end_date),
      x = "Discharge (Q)",
      y = "Concentration (C)",
      color = "Season"  # Legend title for color
    ) +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 14),  
      legend.key.size = unit(1, "lines"),  
      legend.key.height = unit(1, "lines"),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.background = element_rect(fill = "white", color = NA),  # White panel background
      plot.background = element_rect(fill = "white", color = NA)    # White plot background
    )
  
  # Print the plot
  print(cq_plot)
  
  # Save the plot
  ggsave(
    filename = filename,       
    plot = cq_plot,             
    width = 10,                 
    height = 8,                
    dpi = 300                  
  )
}

# Step 5: Process the data and calculate statistics ----
cq_all_processed <- process_cq_data(cq_all, all_sites_order)
stats_table <- calculate_statistics(cq_all_processed)

# Save the statistics table as a CSV
write.csv(stats_table, "cq_stats_by_season.csv", row.names = FALSE)

# Step 6: Plot the data for all sites by season ----
plot_cq_data(cq_all_processed, "cq_plot_ALL_by_Season.png", "All Data by Season")
