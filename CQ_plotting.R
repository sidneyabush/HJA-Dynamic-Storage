# Load necessary libraries
librarian::shelf(tidyverse, dplyr, ggplot2, tidyr)

# Clear environment
rm(list = ls())

# Set the working directory
setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/")

# Step 1: Read in the datasets ----
cq_lookout <- read.csv("cq_data_lookout_filtered.csv")
cq_all <- read.csv("cq_data_all_filtered.csv")

# Define the desired order of sites
all_sites_order <- c("GSWS08", "GSWS07", "GSWS06", "GSMACK", "GSWS03", "GSWS02", "GSWS01", "GSWS09", "GSWS10", "GSLOOK")
lookout_sites_order <- c("GSWS08", "GSWS07", "GSWS06", "GSMACK", "GSWS03", "GSWS02", "GSWS01", "GSLOOK")

# Define a function to calculate average slope and classify significance
process_cq_data <- function(cq_data, site_order) {
  # Ensure the Stream_Name is a factor with the specified order
  cq_data <- cq_data %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))
  
  # Calculate the average slope for each Stream_Name and variable
  average_slope_data <- cq_data %>%
    group_by(Stream_Name, variable) %>%
    summarize(average_slope = mean(slope, na.rm = TRUE)) %>%
    ungroup()
  
  # Join the average slope data back to the original dataset
  cq_data <- cq_data %>%
    left_join(average_slope_data, by = c("Stream_Name", "variable"))
  
  # Calculate p-values and classify significance
  cq_data <- cq_data %>%
    mutate(slope_category = case_when(
      average_slope > 0.1 ~ "Mobilizing (> 0.1)",
      average_slope <= 0.1 & average_slope >= -0.1 ~ "Chemostatic (-0.1 to 0.1)",
      average_slope < -0.1 ~ "Diluting (< -0.1)",
      TRUE ~ "NA"
    ),
    # Set the factor levels for slope_category to ensure the correct legend order
    slope_category = factor(slope_category, levels = c("Mobilizing (> 0.1)", "Chemostatic (-0.1 to 0.1)", "Diluting (< -0.1)"))
    ) %>%
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
  
  return(cq_data)
}

# Define a function to plot the data with title and subtitle, and a white background
plot_cq_data <- function(cq_data, filename, dataset_type) {
  # Get the start and end dates for the subtitle
  start_date <- min(cq_data$Date, na.rm = TRUE)
  end_date <- max(cq_data$Date, na.rm = TRUE)
  
  # Create the plot
  cq_plot <- ggplot(cq_data, aes(x = Qcms, y = value, color = slope_category)) +
    geom_point(alpha = 0.6, size = 1) +  # Points colored by slope category
    geom_smooth(
      data = cq_data, 
      method = "lm", se = FALSE, 
      aes(linetype = significant), 
      color = "black",  # Set the line color to black
      size = 0.5  # Line thickness
    ) +  # Add linear regression lines only for significant relationships
    scale_x_log10() +  # Logarithmic scale for discharge
    scale_y_log10() +  # Logarithmic scale for concentration
    facet_grid(Stream_Name ~ variable, scales = "free") +  # Facet by Stream_Name and variable
    scale_color_manual(
      values = c(
        "Mobilizing (> 0.1)" = "#66C2A5",  # Muted Green
        "Chemostatic (-0.1 to 0.1)" = "#FC8D62",  # Muted Orange
        "Diluting (< -0.1)" = "#8DA0CB",  # Muted Blue
        "NA" = "grey"  # NA values as grey
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Mobilizing Significant" = "solid",
        "Diluting Significant" = "solid",
        "Chemostatic Significant" = "solid"
      )
    ) +
    labs(
      title = paste("C-Q Relationships for", dataset_type),
      subtitle = paste("Date Range:", start_date, "to", end_date),
      x = "Discharge (Q)",
      y = "Concentration (C)",
      color = "C-Q Behavior",  # Legend title for color
      linetype = "Significance"  # Legend title for line type
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

# Step 2: Process and plot the data for both datasets ----
cq_all_processed <- process_cq_data(cq_all, all_sites_order)
plot_cq_data(cq_all_processed, "cq_plot_ALL.png", "All Data")

cq_lookout_processed <- process_cq_data(cq_lookout, lookout_sites_order)
plot_cq_data(cq_lookout_processed, "cq_plot_LOOKOUT.png", "Lookout Creek Data")
