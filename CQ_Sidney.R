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

# Verify the structure of cq_data after joining
str(cq_data)

# Step 4: Plot CQ for each Stream_Name and variable, with lines colored by average slope ----

cq_plot <- ggplot(cq_data, aes(x = Qcms, y = value, color = average_slope)) +
  geom_point(alpha = 0.6) +  # Scatter plot with transparency
  geom_smooth(method = "lm", se = FALSE, aes(group = interaction(Stream_Name, variable)), color = "black") +  # Add linear regression lines
  scale_x_log10() +  # Logarithmic scale for discharge
  scale_y_log10() +  # Logarithmic scale for concentration
  facet_grid(Stream_Name ~ variable, scales = "free") +  # Facet by Stream_Name and variable
  scale_color_gradient2(low = "blue", mid = "yellow", high = "green", midpoint = 0, na.value = "grey") +  # Gradient color based on average slope, change later for chemostatic, mobilizing, diluting... 
  labs(
    x = "Discharge (Q)",
    y = "Concentration (C)",
    title = "Concentration-Discharge (C-Q) Relationships",
    subtitle = "Grouped by Site and Variable with Line Coloring Based on Average Slope",
    color = "Average Slope"
  ) +
  theme_minimal()

# Print the plot
print(cq_plot)
