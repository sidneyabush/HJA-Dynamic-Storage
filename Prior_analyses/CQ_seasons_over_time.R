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

# Define seasons based on months
assign_season <- function(month) {
  case_when(
    month %in% c(7, 8) ~ "Summer",        # July, August
    month %in% c(9, 10, 11, 12) ~ "Fall",  # September to December
    month %in% c(1, 2, 3) ~ "Winter",      # January to March
    month %in% c(4, 5, 6) ~ "Spring"       # April to June
  )
}

# Step 2: Add Year, Season, and process data ----
process_cq_data <- function(cq_data, site_order) {
  # Ensure the Stream_Name is a factor with the specified order
  cq_data <- cq_data %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))
  
  # Add a Year and Season column based on the Date
  cq_data <- cq_data %>%
    mutate(
      Date = as.Date(Date),
      Year = year(Date),
      Month = month(Date),
      Season = assign_season(Month)
    )
  
  # Calculate the slope for each Stream_Name, variable, year, and season
  slope_data <- cq_data %>%
    group_by(Stream_Name, variable, Year, Season) %>%
    summarize(slope = tryCatch({
      if (n() > 1) {
        coef(lm(log10(value) ~ log10(Qcms)))[2]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_), .groups = 'drop')
  
  # Filter out rows where the slope or variable is NA
  slope_data <- slope_data %>%
    filter(!is.na(slope), !is.na(variable))
  
  return(slope_data)
}

# Step 3: Plotting function for yearly slope points by site (all seasons on the same plot) ----
plot_slope_by_site <- function(slope_data, site, dataset_type) {
  # Filter the data by site
  site_data <- slope_data %>%
    filter(Stream_Name == site)
  
  # Ensure there is data to plot
  if (nrow(site_data) == 0) {
    warning(paste("No data available for site:", site))
    return(NULL)
  }
  
  # Create the plot (slope points for each year)
  slope_plot <- ggplot(site_data, aes(x = Year, y = slope, color = Season)) +
    geom_point(size = 2, alpha = 0.8) +  
    geom_line()+
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +  # Facet by variable
    scale_color_manual(
      values = c(
        "Summer" = "#E69F00",    # Orange for Summer
        "Fall" = "#009E73",      # Green for Fall
        "Winter" = "#56B4E9",    # Blue for Winter
        "Spring" = "#F0E442"     # Yellow for Spring
      )
    ) +
    labs(
      title = paste("C-Q Slope Over Time for", dataset_type, "at", site),
      x = "Year",
      y = "Slope (b)",
      color = "Season"
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
  filename <- paste0("slope_over_time_", site, "_all_seasons.png")
  ggsave(
    filename = filename,       
    plot = slope_plot,             
    width = 10,                 
    height = 8,                
    dpi = 300                  
  )
}

# Step 4: Process the data ----
slope_data <- process_cq_data(cq_all, all_sites_order)

# Step 5: Plot the slope over time for each site, combining all seasons ----
for (site in all_sites_order) {
  plot_slope_by_site(slope_data, site, "All Data")
}
