# Load necessary libraries
librarian::shelf(tidyverse, dplyr, ggplot2, tidyr, lubridate)

# Clear environment
rm(list = ls())

# Set the working directory
setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/")

# Step 1: Read in the dataset for all sites and remove "TN" and "TP" ----
cq_all <- read.csv("cq_data_all_filtered.csv") %>%
  filter(!variable %in% c("TN", "TP"))  # Remove "TN" and "TP" variables

# Define the desired order of sites for all data
all_sites_order <- c("GSWS08", "GSWS07", "GSWS06", "GSMACK", "GSWS03", "GSWS02", "GSWS01", "GSWS09", "GSWS10", "GSLOOK")

# Step 2: Add Year and process data ----
process_cq_data <- function(cq_data, site_order) {
  # Ensure the Stream_Name is a factor with the specified order
  cq_data <- cq_data %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))
  
  # Add a Year column based on the Date
  cq_data <- cq_data %>%
    mutate(
      Date = as.Date(Date),
      Year = year(Date)
    )
  
  # Calculate the average slope for each Stream_Name, variable, and year
  slope_data <- cq_data %>%
    group_by(Stream_Name, variable, Year) %>%
    summarize(slope = tryCatch({
      if (n() > 1) {
        coef(lm(log10(value) ~ log10(Qcms)))[2]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_), .groups = 'drop')
  
  return(slope_data)
}

# Step 3: Plotting function for slope over time ----
plot_slope_over_time <- function(slope_data, filename, dataset_type) {
  # Create the plot (slope over time)
  slope_plot <- ggplot(slope_data, aes(x = Year, y = slope, color = variable)) +
    geom_line(aes(group = variable), size = 1) +  # Line plot for slope over time
    geom_point(size = 2) +  # Points at each year
    facet_wrap(~ Stream_Name, scales = "free_y", ncol = 2) +  # Facet by Stream_Name
    labs(
      title = paste("C-Q Slope Over Time for", dataset_type),
      x = "Year",
      y = "Slope (b)",
      color = "Solute"
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
  print(slope_plot)
  
  # Save the plot
  ggsave(
    filename = filename,       
    plot = slope_plot,             
    width = 12,                 
    height = 8,                
    dpi = 300                  
  )
}

# Step 4: Process the data ----
slope_data <- process_cq_data(cq_all, all_sites_order)

# Step 5: Plot the slope over time for all sites and solutes ----
plot_slope_over_time(slope_data, "slope_over_time_ALL.png", "All Data")
