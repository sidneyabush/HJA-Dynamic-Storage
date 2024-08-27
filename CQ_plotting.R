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
  #dplyr::filter(!variable %in% c("pH")) %>%  # Remove variables not interested in (e.g., pH)
  dplyr::mutate(variable = dplyr::recode(variable,
                                         "dissolved org C" = "DOC",
                                         "specific conductivity" = "spc")) %>%  # Replace with actual variable names
  tidyr::drop_na()  # Remove rows with any NA values


# Step 2: Calculate the average slope for each Stream_Name and variable ----
average_slope_data <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarize(average_slope = mean(slope, na.rm = TRUE)) %>%
  ungroup()

# Join the average slope data back to the original dataset
cq_data <- cq_data %>%
  left_join(average_slope_data, by = c("Stream_Name", "variable"))

# Step 4: Calculate p-values and classify significance ----
cq_data <- cq_data %>%
  mutate(slope_category = case_when(
    average_slope > 0.1 ~ "Mobilizing (> 0.1)",
    average_slope <= 0.1 & average_slope >= -0.1 ~ "Chemostatic (-0.1 to 0.1)",
    average_slope < -0.1 ~ "Diluting (< -0.1)",
    TRUE ~ "NA"
  )) %>%
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

# Step 5: Plotting the data ----
# Plotting the data with black lines for significant relationships ----
cq_plot <- ggplot(cq_data, aes(x = Qcms, y = value, color = slope_category)) +
  geom_point(alpha = 0.6, size = 1) +  # Points colored by slope category
  geom_smooth(
    data = cq_data , 
    method = "lm", se = FALSE, 
    aes(linetype = significant), 
    color = "black",  # Set the line color to black
    size = .5  # Line thickness
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
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# Print the plot
print(cq_plot)

# Save the plot
ggsave(
  filename = "cq_plot_ALL.png",       
  plot = cq_plot,             
  width = 10,                 
  height = 8,                
  dpi = 300                  
)
