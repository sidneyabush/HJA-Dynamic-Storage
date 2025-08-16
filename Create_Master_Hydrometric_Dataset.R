library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)

rm(list = ls())

parse_my_date <- function(d) {
  # turn everything into strings
  dd <- as.character(d)
  
  # detect which look like Excel serials (pure digits)
  is_serial <- grepl("^[0-9]+$", dd)
  
  # prepare output
  out <- rep(as.Date(NA), length(dd))
  
  # 1) Excel serials → Date
  if (any(is_serial)) {
    out[is_serial] <- as_date(
      as.numeric(dd[is_serial]),
      origin = "1899-12-30"
    )
  }
  
  # 2) everything else → lubridate
  if (any(!is_serial)) {
    out[!is_serial] <- parse_date_time(
      dd[!is_serial],
      orders = c("Ymd", "Y-m-d", "mdy", "m/d/Y", "dmy", "d/m/Y", "dbY"),
      quiet = TRUE
    ) %>% as_date()
  }
  
  return(out)
}

make_inter_long <- function(fname, var) {
  raw <- read_csv(file.path(met_dir, fname), show_col_types = FALSE) %>%
    mutate(DATE = parse_my_date(DATE))
  
  long <- raw %>%
    pivot_longer(
      cols           = ends_with("_inter"),
      names_to       = "SITECODE",
      names_pattern  = "(.*)_inter$",
      values_to      = var,
      values_drop_na = TRUE    # keep true zeros, only drop actual NAs
    )
  
  all_dates     <- seq.Date(min(raw$DATE, na.rm = TRUE),
                            max(raw$DATE, na.rm = TRUE),
                            by = "day")
  station_names <- unique(long$SITECODE)
  
  expand_grid(DATE = all_dates, SITECODE = station_names) %>%
    left_join(long, by = c("DATE", "SITECODE")) %>%
    arrange(SITECODE, DATE)
}

read_mack_precip <- function(fname) {
  read_csv(file.path(met_dir, fname), show_col_types = FALSE) %>%
    mutate(
      DATE     = parse_my_date(DATE),
      SITECODE = recode(SITECODE, "GSWSMC" = "GSMACK")
    ) %>%
    filter(SITECODE == "GSMACK") %>%
    select(DATE, SITECODE, PRECIP_TOT_DAY) %>%
    rename(P_mm_d = PRECIP_TOT_DAY)
}

calculate_vpd <- function(temp_celsius, rh_percent) {
  # Calculate saturation vapor pressure (kPa) using temperature
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  ea <- es * (rh_percent / 100)
  vpd <- es - ea
  return(vpd)
}

create_relationship_plot <- function(data, site1, site2, variable, r_squared, complete_count) {
  site1_col <- paste0(site1)
  site2_col <- paste0(site2)
  
  p <- ggplot(data, aes_string(x = paste0("`", site1_col, "`"), y = paste0("`", site2_col, "`"))) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = paste("Relationship for", variable, "between stations"),
         subtitle = paste(site1, "and", site2, "(n =", complete_count, ")"),
         x = site1, y = site2) +
    annotate("text", 
             x = min(data[[site1_col]], na.rm = TRUE) + 0.8 * (max(data[[site1_col]], na.rm = TRUE) - min(data[[site1_col]], na.rm = TRUE)),
             y = min(data[[site2_col]], na.rm = TRUE) + 0.1 * (max(data[[site2_col]], na.rm = TRUE) - min(data[[site2_col]], na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared),
             hjust = 1, 
             fontface = "bold")
  
  return(p)
}

create_multiple_regression_plot <- function(target_site, predictor_sites, variable, model_summary, complete_count) {
  title <- paste("Multiple Regression Model for", variable)
  subtitle <- paste(target_site, "predicted from", paste(predictor_sites, collapse = ", "))
  model_data <- environment(model_summary$call)$data
  model_info <- data.frame(
    Metric = c("R²", "Adjusted R²", "F-statistic", "p-value", "Sample Size"),
    Value = c(
      round(model_summary$r.squared, 4),
      round(model_summary$adj.r.squared, 4),
      round(model_summary$fstatistic[1], 2),
      format.pval(
        pf(model_summary$fstatistic[1], 
           model_summary$fstatistic[2], 
           model_summary$fstatistic[3], 
           lower.tail = FALSE),
        digits = 3
      ),
      complete_count
    )
  )
  
  # Prepare coefficient information
  coef_info <- data.frame(
    Term = rownames(model_summary$coefficients),
    Estimate = model_summary$coefficients[, "Estimate"],
    `Std. Error` = model_summary$coefficients[, "Std. Error"],
    `t value` = model_summary$coefficients[, "t value"],
    `Pr(>|t|)` = model_summary$coefficients[, "Pr(>|t|)"]
  )
  
  # Create the plot
  p <- ggplot() +
    theme_minimal(base_size = 12) +
    annotate("text", x = 0.5, y = 0.9, label = title, 
             fontface = "bold", size = 5, hjust = 0.5) +
    annotate("text", x = 0.5, y = 0.85, label = subtitle, 
             fontface = "italic", size = 4, hjust = 0.5) +
    annotate("text", x = 0.5, y = 0.75, label = "Model Metrics:", 
             fontface = "bold", size = 4, hjust = 0.5) +
    annotate("text", x = 0.3, y = 0.7, label = paste(model_info$Metric, collapse = "\n"), 
             hjust = 1, size = 3) +
    annotate("text", x = 0.7, y = 0.7, label = paste(model_info$Value, collapse = "\n"), 
             hjust = 0, size = 3) +
    annotate("text", x = 0.5, y = 0.5, label = "Coefficients:", 
             fontface = "bold", size = 4, hjust = 0.5) +
    annotate("text", x = 0.2, y = 0.45, 
             label = paste(c("Term", coef_info$Term), collapse = "\n"), 
             hjust = 0, fontface = "bold", size = 3) +
    annotate("text", x = 0.4, y = 0.45, 
             label = paste(c("Estimate", round(coef_info$Estimate, 4)), collapse = "\n"), 
             hjust = 0, size = 3) +
    annotate("text", x = 0.6, y = 0.45, 
             label = paste(c("p-value", format.pval(coef_info$`Pr(>|t|)`, digits = 3)), collapse = "\n"), 
             hjust = 0, size = 3) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = title, subtitle = subtitle)
  
  return(p)
}

plot_triplet_station_comparisons <- function(interpolated_data) {
  triplet_stations <- c("WS7MET", "VANMET", "H15MET")
  
  triplet_rh_data <- interpolated_data %>%
    filter(SITECODE %in% triplet_stations) %>%
    select(DATE, SITECODE, RH_d_pct) %>%
    pivot_wider(names_from = SITECODE, values_from = RH_d_pct)
  
  # Remove rows where we don't have all stations data
  complete_triplet_data <- triplet_rh_data %>%
    filter(!is.na(WS7MET) & !is.na(VANMET) & !is.na(H15MET))
  
  # Create WS7MET vs VANMET plot
  p1 <- ggplot(complete_triplet_data, aes(x = WS7MET, y = VANMET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = "Comparison of Relative Humidity",
         subtitle = "WS7MET vs VANMET",
         x = "WS7MET Relative Humidity (%)",
         y = "VANMET Relative Humidity (%)")
  
  # Calculate R-squared for WS7MET vs VANMET
  model_ws7_van <- lm(VANMET ~ WS7MET, data = complete_triplet_data)
  r_squared_ws7_van <- summary(model_ws7_van)$r.squared
  
  p1 <- p1 + 
    annotate("text", 
             x = min(complete_triplet_data$WS7MET, na.rm = TRUE) + 
               0.8 * (max(complete_triplet_data$WS7MET, na.rm = TRUE) - 
                        min(complete_triplet_data$WS7MET, na.rm = TRUE)),
             y = min(complete_triplet_data$VANMET, na.rm = TRUE) + 
               0.1 * (max(complete_triplet_data$VANMET, na.rm = TRUE) - 
                        min(complete_triplet_data$VANMET, na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared_ws7_van),
             hjust = 1, 
             fontface = "bold")
  
  # Create WS7MET vs H15MET plot
  p2 <- ggplot(complete_triplet_data, aes(x = WS7MET, y = H15MET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = "Comparison of Relative Humidity",
         subtitle = "WS7MET vs H15MET",
         x = "WS7MET Relative Humidity (%)",
         y = "H15MET Relative Humidity (%)")
  
  # Calculate R-squared for WS7MET vs H15MET
  model_ws7_h15 <- lm(H15MET ~ WS7MET, data = complete_triplet_data)
  r_squared_ws7_h15 <- summary(model_ws7_h15)$r.squared
  p2 <- p2 + 
    annotate("text", 
             x = min(complete_triplet_data$WS7MET, na.rm = TRUE) + 
               0.8 * (max(complete_triplet_data$WS7MET, na.rm = TRUE) - 
                        min(complete_triplet_data$WS7MET, na.rm = TRUE)),
             y = min(complete_triplet_data$H15MET, na.rm = TRUE) + 
               0.1 * (max(complete_triplet_data$H15MET, na.rm = TRUE) - 
                        min(complete_triplet_data$H15MET, na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared_ws7_h15),
             hjust = 1, 
             fontface = "bold")
  
  # Also create VANMET vs H15MET plot 
  p3 <- ggplot(complete_triplet_data, aes(x = VANMET, y = H15MET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = "Comparison of Relative Humidity",
         subtitle = "VANMET vs H15MET",
         x = "VANMET Relative Humidity (%)",
         y = "H15MET Relative Humidity (%)")
  
  # Calculate R-squared for VANMET vs H15MET
  model_van_h15 <- lm(H15MET ~ VANMET, data = complete_triplet_data)
  r_squared_van_h15 <- summary(model_van_h15)$r.squared
  p3 <- p3 + 
    annotate("text", 
             x = min(complete_triplet_data$VANMET, na.rm = TRUE) + 
               0.8 * (max(complete_triplet_data$VANMET, na.rm = TRUE) - 
                        min(complete_triplet_data$VANMET, na.rm = TRUE)),
             y = min(complete_triplet_data$H15MET, na.rm = TRUE) + 
               0.1 * (max(complete_triplet_data$H15MET, na.rm = TRUE) - 
                        min(complete_triplet_data$H15MET, na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared_van_h15),
             hjust = 1, 
             fontface = "bold")
  
  # Create time series comparison of all three stations
  triplet_long <- triplet_rh_data %>%
    pivot_longer(cols = c("WS7MET", "VANMET", "H15MET"), 
                 names_to = "Station", 
                 values_to = "RH_pct")
  
  p4 <- ggplot(triplet_long, aes(x = DATE, y = RH_pct, color = Station)) +
    geom_line(linewidth = 0.5) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom"
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Time Series Comparison of Relative Humidity",
         subtitle = "Triplet Stations: WS7MET, VANMET, and H15MET",
         x = "Date",
         y = "Relative Humidity (%)")
  
  # Save the plots
  ggsave(file.path(output_dir, "plots", "RH_comparison_WS7MET_vs_VANMET.png"), 
         plot = p1, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "plots", "RH_comparison_WS7MET_vs_H15MET.png"), 
         plot = p2, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "plots", "RH_comparison_VANMET_vs_H15MET.png"), 
         plot = p3, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "plots", "RH_timeseries_triplet_stations.png"), 
         plot = p4, width = 12, height = 6, dpi = 300)
  
  # Print model summaries
  cat("\nRegression model for WS7MET vs VANMET:\n")
  print(summary(model_ws7_van))
  
  cat("\nRegression model for WS7MET vs H15MET:\n")
  print(summary(model_ws7_h15))
  
  cat("\nRegression model for VANMET vs H15MET:\n")
  print(summary(model_van_h15))
  
  # Return the models in case they're needed
  return(list(
    model_ws7_van = model_ws7_van,
    model_ws7_h15 = model_ws7_h15,
    model_van_h15 = model_van_h15
  ))
}

# Extract all required station pairs and triplets that need interpolation
extract_station_groups <- function(site_mapping) {
  pairs <- list()
  triplets <- list()
  
  for (site in names(site_mapping)) {
    site_info <- site_mapping[[site]]
    
    # Check each variable for groups of stations
    for (var in c("temp", "precip", "rh", "netrad")) {
      stations <- site_info[[var]]
      
      if (length(stations) == 2) {
        # This is a pair that needs interpolation
        pair_name <- paste(stations[1], stations[2], sep = "_")
        if (!pair_name %in% names(pairs)) {
          pairs[[pair_name]] <- list(site1 = stations[1], site2 = stations[2])
        }
      } else if (length(stations) == 3) {
        # This is a triplet that needs multiple regression
        triplet_name <- paste(stations, collapse = "_")
        if (!triplet_name %in% names(triplets)) {
          triplets[[triplet_name]] <- list(
            site1 = stations[1], 
            site2 = stations[2], 
            site3 = stations[3]
          )
        }
      }
    }
  }
  
  return(list(
    pairs = pairs,
    triplets = triplets
  ))
}

# Extract all unique stations needed
extract_stations <- function(site_mapping) {
  stations <- c()
  
  for (site in names(site_mapping)) {
    site_info <- site_mapping[[site]]
    
    # Add stations from each variable
    for (var in c("temp", "precip", "rh", "netrad")) {
      stations <- c(stations, site_info[[var]])
    }
  }
  
  return(unique(stations))
}

# Create a function to interpolate missing values based on linear relationship between two sites
interpolate_pair <- function(data, site1, site2, variable) {
  # Filter data for the two sites and check for duplicates
  check_dupes <- data %>% 
    filter(SITECODE %in% c(site1, site2)) %>%
    group_by(DATE, SITECODE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)
  
  if(nrow(check_dupes) > 0) {
    print(paste("Warning: Found", nrow(check_dupes), "duplicate date-site combinations"))
    print(check_dupes)
    
    # Summarize duplicates by taking the mean
    data <- data %>%
      group_by(DATE, SITECODE) %>%
      summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
    
    print("Duplicates have been averaged")
  }
  
  # Now proceed with the filtered data
  pair_data <- data %>% 
    filter(SITECODE %in% c(site1, site2)) %>%
    select(DATE, SITECODE, !!sym(variable)) %>%
    pivot_wider(names_from = SITECODE, values_from = !!sym(variable), values_fn = mean)
  
  # Create column names
  site1_col <- paste0(site1)
  site2_col <- paste0(site2)
  
  # Check if there are enough data points to build a model
  complete_rows <- complete.cases(pair_data[, c(site1_col, site2_col)])
  complete_count <- sum(complete_rows)
  
  if(complete_count < 5) {
    warning(paste("Not enough complete cases for", site1, "and", site2, "for variable", variable, 
                  "(only", complete_count, "overlapping records)"))
    # Return the original data without interpolation
    return(data)
  }
  
  # Fit linear model using complete cases only
  model <- lm(formula = paste0("`", site2_col, "` ~ `", site1_col, "`"), 
              data = pair_data[complete_rows, ])
  
  # Get model stats
  r_squared <- summary(model)$r.squared
  p_value <- summary(model)$coefficients[2,4]  # p-value for the slope
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  
  print(paste("Linear model for", variable, "between", site1, "and", site2, ":"))
  print(paste("  Formula:", site2, "=", round(intercept, 4), "+", round(slope, 4), "*", site1))
  print(paste("  R-squared:", round(r_squared, 4)))
  print(paste("  p-value:", format.pval(p_value, digits = 4)))
  print(paste("  Sample size:", complete_count, "overlapping records"))
  
  # Optional: Visualize relationship
  plot_data <- pair_data[complete_rows, ]
  
  # Create and save a nice looking plot
  p <- create_relationship_plot(plot_data, site1, site2, variable, r_squared, complete_count)
  
  # Print in the console
  print(p)
  
  # Save the plot to a PNG file
  plot_filename <- file.path(output_dir, "plots", paste0("relationship_", variable, "_", site1, "_", site2, ".png"))
  ggsave(plot_filename, plot = p, width = 8, height = 6, dpi = 300)
  
  # Apply interpolation for site2 where site1 has data but site2 is missing
  s1_has_data_s2_missing <- !is.na(pair_data[[site1_col]]) & is.na(pair_data[[site2_col]])
  if(any(s1_has_data_s2_missing)) {
    predictions <- intercept + slope * pair_data[[site1_col]][s1_has_data_s2_missing]
    
    # Create a temporary dataframe for the interpolated values
    temp_df <- data.frame(
      DATE = pair_data$DATE[s1_has_data_s2_missing],
      SITECODE = site2,
      value = predictions
    )
    names(temp_df)[3] <- variable
    
    # Add interpolated values to original data
    data <- bind_rows(
      data,
      temp_df
    ) %>% distinct(DATE, SITECODE, .keep_all = TRUE)
  }
  
  # Apply interpolation for site1 where site2 has data but site1 is missing
  s2_has_data_s1_missing <- is.na(pair_data[[site1_col]]) & !is.na(pair_data[[site2_col]])
  if(any(s2_has_data_s1_missing)) {
    # Use inverse relationship
    inv_intercept <- -intercept/slope
    inv_slope <- 1/slope
    predictions <- inv_intercept + inv_slope * pair_data[[site2_col]][s2_has_data_s1_missing]
    
    # Create a temporary dataframe for the interpolated values
    temp_df <- data.frame(
      DATE = pair_data$DATE[s2_has_data_s1_missing],
      SITECODE = site1,
      value = predictions
    )
    names(temp_df)[3] <- variable
    
    # Add interpolated values to original data
    data <- bind_rows(
      data,
      temp_df
    ) %>% distinct(DATE, SITECODE, .keep_all = TRUE)
  }
  
  return(data)
}

# Function to interpolate missing values in a triplet of sites using multiple regression
interpolate_triplet <- function(data, site1, site2, site3, variable) {
  # Filter data for the three sites and check for duplicates
  check_dupes <- data %>% 
    filter(SITECODE %in% c(site1, site2, site3)) %>%
    group_by(DATE, SITECODE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)
  
  if(nrow(check_dupes) > 0) {
    print(paste("Warning: Found", nrow(check_dupes), "duplicate date-site combinations"))
    print(check_dupes)
    
    # Summarize duplicates by taking the mean
    data <- data %>%
      group_by(DATE, SITECODE) %>%
      summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
    
    print("Duplicates have been averaged")
  }
  
  # Now proceed with the filtered data
  triplet_data <- data %>% 
    filter(SITECODE %in% c(site1, site2, site3)) %>%
    select(DATE, SITECODE, !!sym(variable)) %>%
    pivot_wider(names_from = SITECODE, values_from = !!sym(variable), values_fn = mean)
  
  # Create column names
  site1_col <- paste0(site1)
  site2_col <- paste0(site2)
  site3_col <- paste0(site3)
  
  # For each site, create a multiple regression model using the other two sites as predictors
  # 1. Model for site1 using site2 and site3
  model1_rows <- complete.cases(triplet_data[, c(site1_col, site2_col, site3_col)])
  model1_count <- sum(model1_rows)
  
  if(model1_count < 10) {
    warning(paste("Not enough complete cases for triplet", site1, site2, site3, "for variable", variable, 
                  "(only", model1_count, "overlapping records)"))
    # Return the original data without interpolation
    return(data)
  }
  
  # Create the three models
  model1 <- lm(formula = paste0("`", site1_col, "` ~ `", site2_col, "` + `", site3_col, "`"), 
               data = triplet_data[model1_rows, ])
  
  model2 <- lm(formula = paste0("`", site2_col, "` ~ `", site1_col, "` + `", site3_col, "`"), 
               data = triplet_data[model1_rows, ])
  
  model3 <- lm(formula = paste0("`", site3_col, "` ~ `", site1_col, "` + `", site2_col, "`"), 
               data = triplet_data[model1_rows, ])
  
  # Get model summaries
  model1_summary <- summary(model1)
  model2_summary <- summary(model2)
  model3_summary <- summary(model3)
  
  # Print model statistics
  print(paste("Multiple regression model for", variable, "- predicting", site1, "from", site2, "and", site3, ":"))
  print(paste("  R-squared:", round(model1_summary$r.squared, 4)))
  print(paste("  Adjusted R-squared:", round(model1_summary$adj.r.squared, 4)))
  print(paste("  Sample size:", model1_count, "complete records"))
  
  print(paste("Multiple regression model for", variable, "- predicting", site2, "from", site1, "and", site3, ":"))
  print(paste("  R-squared:", round(model2_summary$r.squared, 4)))
  print(paste("  Adjusted R-squared:", round(model2_summary$adj.r.squared, 4)))
  
  print(paste("Multiple regression model for", variable, "- predicting", site3, "from", site1, "and", site2, ":"))
  print(paste("  R-squared:", round(model3_summary$r.squared, 4)))
  print(paste("  Adjusted R-squared:", round(model3_summary$adj.r.squared, 4)))
  
  # Create and save plots for the models
  p1 <- create_multiple_regression_plot(site1, c(site2, site3), variable, model1_summary, model1_count)
  p2 <- create_multiple_regression_plot(site2, c(site1, site3), variable, model2_summary, model1_count)
  p3 <- create_multiple_regression_plot(site3, c(site1, site2), variable, model3_summary, model1_count)
  
  # Save the plots
  plot_filename1 <- file.path(output_dir, "plots", paste0("multireg_", variable, "_", site1, "_from_", site2, "_", site3, ".png"))
  plot_filename2 <- file.path(output_dir, "plots", paste0("multireg_", variable, "_", site2, "_from_", site1, "_", site3, ".png"))
  plot_filename3 <- file.path(output_dir, "plots", paste0("multireg_", variable, "_", site3, "_from_", site1, "_", site2, ".png"))
  
  ggsave(plot_filename1, plot = p1, width = 8, height = 6, dpi = 300)
  ggsave(plot_filename2, plot = p2, width = 8, height = 6, dpi = 300)
  ggsave(plot_filename3, plot = p3, width = 8, height = 6, dpi = 300)
  
  # Now apply the interpolation for missing values
  
  # Case 1: Only site1 is missing, use model1 to predict it
  case1 <- is.na(triplet_data[[site1_col]]) & !is.na(triplet_data[[site2_col]]) & !is.na(triplet_data[[site3_col]])
  if(any(case1)) {
    # Create a temp dataframe for prediction
    temp_df <- triplet_data[case1, c("DATE", site2_col, site3_col)]
    names(temp_df) <- c("DATE", site2_col, site3_col)
    
    # Predict site1 using model1
    predictions <- predict(model1, newdata = temp_df)
    
    # Create dataframe of interpolated values
    interp_df <- data.frame(
      DATE = triplet_data$DATE[case1],
      SITECODE = site1,
      value = predictions
    )
    names(interp_df)[3] <- variable
    
    # Add to the data
    data <- bind_rows(
      data,
      interp_df
    ) %>% distinct(DATE, SITECODE, .keep_all = TRUE)
  }
  
  # Case 2: Only site2 is missing, use model2 to predict it
  case2 <- !is.na(triplet_data[[site1_col]]) & is.na(triplet_data[[site2_col]]) & !is.na(triplet_data[[site3_col]])
  if(any(case2)) {
    # Create a temp dataframe for prediction
    temp_df <- triplet_data[case2, c("DATE", site1_col, site3_col)]
    names(temp_df) <- c("DATE", site1_col, site3_col)
    
    # Predict site2 using model2
    predictions <- predict(model2, newdata = temp_df)
    
    # Create dataframe of interpolated values
    interp_df <- data.frame(
      DATE = triplet_data$DATE[case2],
      SITECODE = site2,
      value = predictions
    )
    names(interp_df)[3] <- variable
    
    # Add to the data
    data <- bind_rows(
      data,
      interp_df
    ) %>% distinct(DATE, SITECODE, .keep_all = TRUE)
  }
  
  # Case 3: Only site3 is missing, use model3 to predict it
  case3 <- !is.na(triplet_data[[site1_col]]) & !is.na(triplet_data[[site2_col]]) & is.na(triplet_data[[site3_col]])
  if(any(case3)) {
    # Create a temp dataframe for prediction
    temp_df <- triplet_data[case3, c("DATE", site1_col, site2_col)]
    names(temp_df) <- c("DATE", site1_col, site2_col)
    
    # Predict site3 using model3
    predictions <- predict(model3, newdata = temp_df)
    
    # Create dataframe of interpolated values
    interp_df <- data.frame(
      DATE = triplet_data$DATE[case3],
      SITECODE = site3,
      value = predictions
    )
    names(interp_df)[3] <- variable
    
    # Add to the data
    data <- bind_rows(
      data,
      interp_df
    ) %>% distinct(DATE, SITECODE, .keep_all = TRUE)
  }
  
  return(data)
}

# Function to constrain interpolated values (cap RH at 100%, precipitation at 0)
constrain_interpolated_values <- function(data) {
  # For RH: Cap values > 100% to exactly 100%
  if ("RH_d_pct" %in% names(data)) {
    over_100_count <- sum(data$RH_d_pct > 100, na.rm = TRUE)
    if (over_100_count > 0) {
      cat(paste0("Capping ", over_100_count, " RH values that were > 100% to exactly 100%\n"))
      data <- data %>%
        mutate(RH_d_pct = ifelse(RH_d_pct > 100, 100, RH_d_pct))
    }
  }
  
  # For precipitation: Cap negative values to 0
  if ("P_mm_d" %in% names(data)) {
    neg_precip_count <- sum(data$P_mm_d < 0, na.rm = TRUE)
    if (neg_precip_count > 0) {
      cat(paste0("Capping ", neg_precip_count, " precipitation values that were < 0 to exactly 0\n"))
      data <- data %>%
        mutate(P_mm_d = ifelse(P_mm_d < 0, 0, P_mm_d))
    }
  }
  
  return(data)
}

# Function to process all pairs, triplets and interpolate data
process_station_groups <- function(data, station_groups, variables) {
  # First process all the pairs
  interpolated_data <- data
  
  # Keep track of which site pairs were actually interpolated
  interpolated_pairs <- list()
  interpolated_triplets <- list()
  
  # Process pairs
  for (pair_name in names(station_groups$pairs)) {
    pair <- station_groups$pairs[[pair_name]]
    site1 <- pair$site1
    site2 <- pair$site2
    
    # Check if we can interpolate this pair by looking at data availability
    pair_can_be_interpolated <- FALSE
    
    for (var in variables) {
      # Extract data for this pair
      pair_data <- data %>% 
        filter(SITECODE %in% c(site1, site2)) %>%
        select(DATE, SITECODE, !!sym(var)) %>%
        pivot_wider(names_from = SITECODE, values_from = !!sym(var), values_fn = mean)
      
      # Count complete cases
      complete_count <- sum(complete.cases(pair_data[, c(site1, site2)]))
      
      if (complete_count >= 5) {
        pair_can_be_interpolated <- TRUE
        break
      }
    }
    
    if (pair_can_be_interpolated) {
      # Perform interpolation for all variables
      for (var in variables) {
        interpolated_data <- interpolate_pair(interpolated_data, site1, site2, var)
      }
      # Add to the list of successful interpolations
      interpolated_pairs[[pair_name]] <- pair
    } else {
      message(paste("Warning: Pair", site1, "and", site2, 
                    "doesn't have enough overlapping data for interpolation."))
      message("Both sites will be included as standalone without interpolation.")
    }
  }
  
  # Process triplets
  for (triplet_name in names(station_groups$triplets)) {
    triplet <- station_groups$triplets[[triplet_name]]
    site1 <- triplet$site1
    site2 <- triplet$site2
    site3 <- triplet$site3
    
    # Check if we can interpolate this triplet by looking at data availability
    triplet_can_be_interpolated <- FALSE
    
    for (var in variables) {
      # Extract data for this triplet
      triplet_data <- data %>% 
        filter(SITECODE %in% c(site1, site2, site3)) %>%
        select(DATE, SITECODE, !!sym(var)) %>%
        pivot_wider(names_from = SITECODE, values_from = !!sym(var), values_fn = mean)
      
      # Count complete cases where all three sites have data
      complete_count <- sum(complete.cases(triplet_data[, c(site1, site2, site3)]))
      
      if (complete_count >= 10) {  # Need more data for multiple regression
        triplet_can_be_interpolated <- TRUE
        break
      }
    }
    
    if (triplet_can_be_interpolated) {
      # Perform interpolation for all variables
      for (var in variables) {
        interpolated_data <- interpolate_triplet(interpolated_data, site1, site2, site3, var)
      }
      # Add to the list of successful interpolations
      interpolated_triplets[[triplet_name]] <- triplet
    } else {
      message(paste("Warning: Triplet", site1, site2, "and", site3, 
                    "doesn't have enough overlapping data for interpolation."))
      message("Sites will be interpolated using pairwise relationships instead.")
      
      # Fall back to pairwise interpolation for the triplet
      pair1_name <- paste(site1, site2, sep = "_")
      pair2_name <- paste(site1, site3, sep = "_")
      pair3_name <- paste(site2, site3, sep = "_")
      
      fallback_pairs <- list(
        list(site1 = site1, site2 = site2),
        list(site1 = site1, site2 = site3),
        list(site1 = site2, site2 = site3)
      )
      
      for (pair in fallback_pairs) {
        for (var in variables) {
          interpolated_data <- interpolate_pair(interpolated_data, pair$site1, pair$site2, var)
        }
      }
    }
  }
  
  # Apply constraints to the interpolated data (cap RH at 100%, precip at 0)
  interpolated_data <- constrain_interpolated_values(interpolated_data)
  
  # After all interpolation, calculate VPD from interpolated T and RH
  # This ensures VPD is calculated using the interpolated temperature and humidity values
  interpolated_data <- interpolated_data %>%
    rowwise() %>%
    mutate(VPD_kPa = if_else(!is.na(T_C) & !is.na(RH_d_pct), 
                             calculate_vpd(T_C, RH_d_pct), 
                             as.numeric(NA))) %>%
    ungroup()
  
  return(list(
    data = interpolated_data,
    interpolated_pairs = interpolated_pairs,
    interpolated_triplets = interpolated_triplets
  ))
}


# Modified function to create ONLY the watershed datasets (no met station data)
create_watershed_datasets <- function(interpolated_data, site_mapping, variables) {
  site_datasets <- list()
  
  for (site_name in names(site_mapping)) {
    site_info <- site_mapping[[site_name]]
    site_data <- data.frame()
    
    # Process each variable
    for (var_idx in seq_along(variables)) {
      var <- variables[var_idx]
      
      # Skip if we're beyond the available mappings or handling VPD differently
      if (var_idx > length(names(site_info)) && var != "VPD_kPa") {
        next
      }
      
      # For VPD, use the same station mapping as temperature
      if (var == "VPD_kPa") {
        var_name <- "temp"  # Use temperature stations for VPD
      } else {
        var_name <- names(site_info)[var_idx]  # temp, precip, rh, netrad
      }
      
      stations <- site_info[[var_name]]
      
      if (length(stations) == 1) {
        # Single station - just extract the data
        single_data <- interpolated_data %>%
          filter(SITECODE == stations[1]) %>%
          select(DATE, !!sym(var))
        
        if (nrow(site_data) == 0) {
          # First variable, create the dataframe
          site_data <- single_data %>% mutate(SITECODE = site_name)
        } else {
          # Add variable to existing dataframe
          site_data <- site_data %>% 
            left_join(single_data, by = "DATE")
        }
      } else if (length(stations) >= 2) {
        # Multiple stations - average the interpolated values
        multi_data <- interpolated_data %>%
          filter(SITECODE %in% stations) %>%
          select(DATE, SITECODE, !!sym(var)) %>%
          pivot_wider(names_from = SITECODE, values_from = !!sym(var)) %>%
          mutate(avg_value = rowMeans(select(., stations), na.rm = TRUE)) %>%
          select(DATE, avg_value)
        
        # Rename the averaged column
        names(multi_data)[2] <- var
        
        if (nrow(site_data) == 0) {
          # First variable, create the dataframe
          site_data <- multi_data %>% mutate(SITECODE = site_name)
        } else {
          # Add variable to existing dataframe
          site_data <- site_data %>% 
            left_join(multi_data, by = "DATE")
        }
      }
    }
    
    # Store the site dataset
    site_datasets[[site_name]] <- site_data
  }
  
  return(site_datasets)
}

# Add this function to create comparison plots for specific triplet stations
plot_triplet_station_comparisons <- function(interpolated_data) {
  # Extract RH data for the specific stations
  triplet_stations <- c("WS7MET", "VANMET", "H15MET")
  
  triplet_rh_data <- interpolated_data %>%
    filter(SITECODE %in% triplet_stations) %>%
    select(DATE, SITECODE, RH_d_pct) %>%
    pivot_wider(names_from = SITECODE, values_from = RH_d_pct)
  
  # Remove rows where we don't have all stations data
  complete_triplet_data <- triplet_rh_data %>%
    filter(!is.na(WS7MET) & !is.na(VANMET) & !is.na(H15MET))
  
  # Create WS7MET vs VANMET plot
  p1 <- ggplot(complete_triplet_data, aes(x = WS7MET, y = VANMET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = "Comparison of Relative Humidity",
         subtitle = "WS7MET vs VANMET",
         x = "WS7MET Relative Humidity (%)",
         y = "VANMET Relative Humidity (%)")
  
  # Calculate R-squared for WS7MET vs VANMET
  model_ws7_van <- lm(VANMET ~ WS7MET, data = complete_triplet_data)
  r_squared_ws7_van <- summary(model_ws7_van)$r.squared
  
  # Add R-squared annotation to plot
  p1 <- p1 + 
    annotate("text", 
             x = min(complete_triplet_data$WS7MET, na.rm = TRUE) + 
               0.8 * (max(complete_triplet_data$WS7MET, na.rm = TRUE) - 
                        min(complete_triplet_data$WS7MET, na.rm = TRUE)),
             y = min(complete_triplet_data$VANMET, na.rm = TRUE) + 
               0.1 * (max(complete_triplet_data$VANMET, na.rm = TRUE) - 
                        min(complete_triplet_data$VANMET, na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared_ws7_van),
             hjust = 1, 
             fontface = "bold")
  
  # Create WS7MET vs H15MET plot
  p2 <- ggplot(complete_triplet_data, aes(x = WS7MET, y = H15MET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = "Comparison of Relative Humidity",
         subtitle = "WS7MET vs H15MET",
         x = "WS7MET Relative Humidity (%)",
         y = "H15MET Relative Humidity (%)")
  
  # Calculate R-squared for WS7MET vs H15MET
  model_ws7_h15 <- lm(H15MET ~ WS7MET, data = complete_triplet_data)
  r_squared_ws7_h15 <- summary(model_ws7_h15)$r.squared
  
  # Add R-squared annotation to plot
  p2 <- p2 + 
    annotate("text", 
             x = min(complete_triplet_data$WS7MET, na.rm = TRUE) + 
               0.8 * (max(complete_triplet_data$WS7MET, na.rm = TRUE) - 
                        min(complete_triplet_data$WS7MET, na.rm = TRUE)),
             y = min(complete_triplet_data$H15MET, na.rm = TRUE) + 
               0.1 * (max(complete_triplet_data$H15MET, na.rm = TRUE) - 
                        min(complete_triplet_data$H15MET, na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared_ws7_h15),
             hjust = 1, 
             fontface = "bold")
  
  # Also create VANMET vs H15MET plot for completeness
  p3 <- ggplot(complete_triplet_data, aes(x = VANMET, y = H15MET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = "Comparison of Relative Humidity",
         subtitle = "VANMET vs H15MET",
         x = "VANMET Relative Humidity (%)",
         y = "H15MET Relative Humidity (%)")
  
  # Calculate R-squared for VANMET vs H15MET
  model_van_h15 <- lm(H15MET ~ VANMET, data = complete_triplet_data)
  r_squared_van_h15 <- summary(model_van_h15)$r.squared
  
  # Add R-squared annotation to plot
  p3 <- p3 + 
    annotate("text", 
             x = min(complete_triplet_data$VANMET, na.rm = TRUE) + 
               0.8 * (max(complete_triplet_data$VANMET, na.rm = TRUE) - 
                        min(complete_triplet_data$VANMET, na.rm = TRUE)),
             y = min(complete_triplet_data$H15MET, na.rm = TRUE) + 
               0.1 * (max(complete_triplet_data$H15MET, na.rm = TRUE) - 
                        min(complete_triplet_data$H15MET, na.rm = TRUE)),
             label = sprintf("R² = %.3f", r_squared_van_h15),
             hjust = 1, 
             fontface = "bold")
  
  # Create time series comparison of all three stations
  triplet_long <- triplet_rh_data %>%
    pivot_longer(cols = c("WS7MET", "VANMET", "H15MET"), 
                 names_to = "Station", 
                 values_to = "RH_pct")
  
  p4 <- ggplot(triplet_long, aes(x = DATE, y = RH_pct, color = Station)) +
    geom_line(linewidth = 0.5) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom"
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Time Series Comparison of Relative Humidity",
         subtitle = "Triplet Stations: WS7MET, VANMET, and H15MET",
         x = "Date",
         y = "Relative Humidity (%)")
  
  # Save the plots
  ggsave(file.path(output_dir, "plots", "RH_comparison_WS7MET_vs_VANMET.png"), 
         plot = p1, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "plots", "RH_comparison_WS7MET_vs_H15MET.png"), 
         plot = p2, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "plots", "RH_comparison_VANMET_vs_H15MET.png"), 
         plot = p3, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "plots", "RH_timeseries_triplet_stations.png"), 
         plot = p4, width = 12, height = 6, dpi = 300)
  
  # Print model summaries
  cat("\nRegression model for WS7MET vs VANMET:\n")
  print(summary(model_ws7_van))
  
  cat("\nRegression model for WS7MET vs H15MET:\n")
  print(summary(model_ws7_h15))
  
  cat("\nRegression model for VANMET vs H15MET:\n")
  print(summary(model_van_h15))
  
  # Return the models in case they're needed
  return(list(
    model_ws7_van = model_ws7_van,
    model_ws7_h15 = model_ws7_h15,
    model_van_h15 = model_van_h15
  ))
}

# After processing all sites, create a summary plot for each site
create_site_summary_plots <- function(site_datasets, site_mapping, variables) {
  for (site_name in names(site_datasets)) {
    site_data <- site_datasets[[site_name]]
    site_info <- site_mapping[[site_name]]
    
    # Create a time series plot for each variable
    for (var_idx in seq_along(variables)) {
      var <- variables[var_idx]
      
      # Handle the case where var_idx exceeds the number of variable mappings
      if (var_idx <= length(names(site_info)) || var == "VPD_kPa" || var == "Q_mm_d") {
        # For VPD, use the temperature stations mapping
        if (var == "VPD_kPa") {
          var_name <- "temp"
          stations <- site_info[["temp"]]
        } else if (var == "Q_mm_d") {
          var_name <- "discharge"
          stations <- "N/A"
        } else {
          var_name <- names(site_info)[var_idx]  # temp, precip, rh, netrad
          stations <- site_info[[var_name]]
        }
      } else {
        stations <- "N/A"
      }
      
      # Get descriptive names for the variables
      var_labels <- c(
        "T_C" = "Temperature (°C)",
        "P_mm_d" = "Precipitation (mm/day)",
        "RH_d_pct" = "Relative Humidity (%)",
        "NR_Wm2_d" = "Net Radiation (W/m²)",
        "VPD_kPa" = "Vapor Pressure Deficit (kPa)",
        "Q_mm_d" = "Discharge (mm/day)"
      )
      
      # Check if variable exists in the dataset
      if (!(var %in% names(site_data))) {
        next  # Skip this variable if it's not in the dataset
      }
      
      # Create a better title with relevant info
      title <- paste0(site_name, " ", var_labels[var])
      
      if (var == "Q_mm_d") {
        subtitle <- "Discharge data from gauging station"
      } else if (var == "VPD_kPa") {
        if (length(stations) == 1) {
          subtitle <- paste0("Calculated from temperature and RH at ", stations[1])
        } else if (length(stations) == 2) {
          subtitle <- paste0("Calculated from temperature and RH (average of ", paste(stations, collapse = ", "), ")")
        } else {
          subtitle <- paste0("Calculated from temperature and RH (average of ", paste(stations, collapse = ", "), ")")
        }
      } else if (length(stations) == 1) {
        subtitle <- paste0("Data from ", stations[1])
      } else if (length(stations) == 2) {
        subtitle <- paste0("Average of ", stations[1], " and ", stations[2])
      } else {
        subtitle <- paste0("Average of ", paste(stations, collapse = ", "))
      }
      
      # Create the time series plot
      p <- ggplot(site_data, aes(x = DATE, y = !!sym(var))) +
        geom_line(color = "blue", linewidth = 0.5) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black")
        ) +
        labs(title = title,
             subtitle = subtitle,
             x = "Date",
             y = var_labels[var])
      
      # Save the plot
      plot_filename <- file.path(output_dir, "plots", paste0("timeseries_", site_name, "_", var, ".png"))
      ggsave(plot_filename, plot = p, width = 10, height = 6, dpi = 300)
    }
  }
}

# Modified main workflow function to only generate watershed data and include VPD
process_meteorological_data <- function(combined_met, site_mapping, variables) {
  # Clean any duplicates
  combined_met_clean <- combined_met %>%
    group_by(DATE, SITECODE) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
  
  # Extract pairs and triplets that need interpolation
  station_groups <- extract_station_groups(site_mapping)
  
  cat("Station pairs that need interpolation:\n")
  for (pair_name in names(station_groups$pairs)) {
    pair <- station_groups$pairs[[pair_name]]
    cat(sprintf("- %s: %s and %s\n", pair_name, pair$site1, pair$site2))
  }
  
  cat("\nStation triplets that need interpolation:\n")
  for (triplet_name in names(station_groups$triplets)) {
    triplet <- station_groups$triplets[[triplet_name]]
    cat(sprintf("- %s: %s, %s, and %s\n", triplet_name, triplet$site1, triplet$site2, triplet$site3))
  }
  
  # Process pairs and triplets and interpolate data
  # This now also calculates VPD from interpolated T and RH
  results <- process_station_groups(combined_met_clean, station_groups, variables)
  interpolated_data <- results$data
  interpolated_pairs <- results$interpolated_pairs
  interpolated_triplets <- results$interpolated_triplets
  
  # Print which pairs were successfully interpolated
  cat("\nSuccessfully interpolated pairs:\n")
  for (pair_name in names(interpolated_pairs)) {
    pair <- interpolated_pairs[[pair_name]]
    cat(sprintf("- %s and %s\n", pair$site1, pair$site2))
  }
  
  # Print which triplets were successfully interpolated
  cat("\nSuccessfully interpolated triplets:\n")
  for (triplet_name in names(interpolated_triplets)) {
    triplet <- interpolated_triplets[[triplet_name]]
    cat(sprintf("- %s, %s, and %s\n", triplet$site1, triplet$site2, triplet$site3))
  }
  
  # Update variables to include VPD for watershed datasets
  watershed_variables <- c(variables, "VPD_kPa")
  
  # Create site-specific datasets (only watershed data, not met station data)
  watershed_datasets <- create_watershed_datasets(interpolated_data, site_mapping, watershed_variables)
  
  # Combine all watershed datasets into one
  all_watersheds_data <- bind_rows(watershed_datasets)
  
  # Return results
  return(list(
    interpolated_data = interpolated_data,
    watershed_datasets = watershed_datasets,
    all_watersheds_data = all_watersheds_data,
    interpolated_pairs = interpolated_pairs,
    interpolated_triplets = interpolated_triplets
  ))
}

# Function to add discharge data to watershed datasets directly
add_discharge_to_watersheds <- function(watershed_datasets, discharge) {
  # Read drainage area data
  da_df <- read_csv(file.path(met_dir, "drainage_area.csv"))
  
  discharge_processed <- discharge %>%
    left_join(da_df, by = "SITECODE") %>%
    filter(!is.na(DA_M2)) %>%
    mutate(
      DATE   = parse_my_date(DATE),    # handles both “YYYY-MM-DD” and Excel serials
      Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000
    ) %>%
    dplyr::select(DATE, SITECODE, Q_mm_d)
  
  # Add discharge data to each watershed dataset
  for (site_name in names(watershed_datasets)) {
    # Extract discharge data for this site
    site_discharge <- discharge_processed %>%
      filter(SITECODE == site_name) %>%
      select(DATE, Q_mm_d)
    
    # Add discharge data to the watershed dataset
    if (nrow(site_discharge) > 0) {
      watershed_datasets[[site_name]] <- watershed_datasets[[site_name]] %>%
        left_join(site_discharge, by = "DATE")
    }
  }
  
  return(watershed_datasets)
}

###########################
# PART 3: DATA PROCESSING #----
###########################

# Set the data directory
met_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data/all_hydromet"

# Define output directory for results
output_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/MET"

# Create output directories if they don't exist
dir.create(output_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "data"),  showWarnings = FALSE)

# --- helper: add Q_mm_d to each watershed dataset, without re-parsing DATE ---
add_discharge_to_watersheds <- function(watershed_datasets, discharge) {
  da_df <- read_csv(file.path(met_dir, "drainage_area.csv"), show_col_types = FALSE) %>%
    mutate(SITECODE = recode(SITECODE, "GSWSMC" = "GSMACK"))
  
  discharge_processed <- discharge %>%
    left_join(da_df, by = "SITECODE") %>%
    filter(!is.na(DA_M2)) %>%
    mutate(Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000) %>%
    select(DATE, SITECODE, Q_mm_d)
  
  for (site_name in names(watershed_datasets)) {
    site_dis <- discharge_processed %>% filter(SITECODE == site_name)
    if (nrow(site_dis) > 0) {
      watershed_datasets[[site_name]] <-
        watershed_datasets[[site_name]] %>%
        left_join(site_dis, by = "DATE")
    }
  }
  
  watershed_datasets
}

# Load temperature data
Temp <- make_inter_long("Temperature_original_&_filled_1979_2023_v2.csv", "Temp") %>%
  select(DATE, SITECODE, Temp) %>%
  rename(T_C = Temp)

# Load precipitation data
Precip <- make_inter_long("Precipitation_original_&_filled_1979_2023.csv", "Precip") %>%
  select(DATE, SITECODE, Precip) %>%
  rename(P_mm_d = Precip)

# Add GSMACK precipitation data
MACK_Precip <- read_mack_precip("MS00403_v2.csv")
Precip <- bind_rows(Precip, MACK_Precip)

# Load relative humidity data
RH <- read_csv(file.path(met_dir, "MS00102_v9.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  select(SITECODE, DATE, RELHUM_MEAN_DAY) %>%
  rename(RH_d_pct = RELHUM_MEAN_DAY)

# Load net radiation data
NetRad <- read_csv(file.path(met_dir, "MS05025_v3.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  select(SITECODE, DATE, NR_TOT_MEAN_DAY) %>%
  rename(NR_Wm2_d = NR_TOT_MEAN_DAY)

# Combine all meteorological datasets
combined_met <- Temp %>%
  full_join(Precip, by = c("DATE","SITECODE")) %>%
  full_join(RH,     by = c("DATE","SITECODE")) %>%
  full_join(NetRad, by = c("DATE","SITECODE")) %>%
  filter(DATE >= ymd("1997-10-01")) %>%
  arrange(SITECODE, DATE)

# Define the site mapping
site_mapping <- list(
  "GSWS09" = list(temp=c("PRIMET"),     precip=c("PRIMET"),     rh=c("PRIMET","CS2MET"), netrad=c("PRIMET")),
  "GSWS10" = list(temp=c("PRIMET"),     precip=c("PRIMET"),     rh=c("PRIMET","CS2MET"), netrad=c("PRIMET")),
  "GSWS01" = list(temp=c("PRIMET"),     precip=c("PRIMET"),     rh=c("PRIMET","CS2MET"), netrad=c("PRIMET")),
  "GSWS02" = list(temp=c("PRIMET"),     precip=c("PRIMET"),     rh=c("PRIMET","CS2MET"), netrad=c("PRIMET")),
  "GSWS03" = list(temp=c("PRIMET"),     precip=c("PRIMET"),     rh=c("PRIMET","CS2MET"), netrad=c("PRIMET")),
  "GSMACK" = list(temp=c("CENMET","UPLMET"), precip=c("GSMACK","UPLMET"), rh=c("CENMET","UPLMET"), netrad=c("VANMET")),
  "GSWS06" = list(temp=c("H15MET","VANMET"), precip=c("H15MET"),       rh=c("H15MET","VANMET","WS7MET"), netrad=c("VANMET")),
  "GSWS07" = list(temp=c("H15MET","VANMET"), precip=c("H15MET"),       rh=c("H15MET","VANMET","WS7MET"), netrad=c("VANMET")),
  "GSWS08" = list(temp=c("H15MET","VANMET"), precip=c("H15MET"),       rh=c("H15MET","VANMET","WS7MET"), netrad=c("VANMET")),
  "LONGER" = list(temp=c("CENMET"),     precip=c("CENMET"),     rh=c("CENMET"),         netrad=c("VANMET")),
  "COLD"   = list(temp=c("CENMET","UPLMET"), precip=c("CENMET","UPLMET"), rh=c("CENMET","UPLMET"), netrad=c("VANMET"))
)

# Variables to process (only the met vars here—VPD will be computed internally)
variables <- c("T_C","P_mm_d","RH_d_pct","NR_Wm2_d")

# Run interpolation & VPD calculation
results <- process_meteorological_data(combined_met, site_mapping, variables)

# **CORRECT** assignments from the results list:
interpolated_data     <- results$interpolated_data
watershed_datasets    <- results$watershed_datasets
all_watersheds_data   <- results$all_watersheds_data
interpolated_pairs    <- results$interpolated_pairs
interpolated_triplets <- results$interpolated_triplets

# Now that interpolated_data exists, this works:
triplet_models <- plot_triplet_station_comparisons(interpolated_data)

# read + collapse GSWSMC → GSMACK (dropping the halves)
discharge <- read_csv(file.path(met_dir, "HF00402_v14.csv"), show_col_types = FALSE) %>%
  mutate(
    DATE     = parse_my_date(DATE),
    SITECODE = recode(SITECODE,
                      "GSWSMC" = "GSMACK")) %>%
  filter(!SITECODE %in% c("GSWSMA", "GSWSMF")) %>%
  group_by(DATE, SITECODE) %>%
  summarise(
    MEAN_Q = sum(MEAN_Q, na.rm = TRUE),
    .groups = "drop"
  )

# stitch Q_mm_d into each watershed
watershed_datasets <- add_discharge_to_watersheds(watershed_datasets, discharge)

# **and then** immediately rebuild your master table so Q_mm_d shows up:
all_watersheds_data <- dplyr::bind_rows(watershed_datasets)

all_watersheds_data <- bind_rows(watershed_datasets) %>%
  # remove the stray SITECODE column
  select(-SITECODE, -SITECODE.y) %>%
  # rename the joined‐in SITECODE.x back to SITECODE
  rename(SITECODE = SITECODE.x) %>%
  # optional: reorder so DATE, SITECODE come first
  select(DATE, SITECODE, everything())


# Update variable list for the final NA‐counts and summaries
variables <- c("T_C", "P_mm_d", "RH_d_pct", "NR_Wm2_d", "VPD_kPa", "Q_mm_d")

# 1) Read drainage areas and recode GSLOOK → GSLOOK_FULL
da_df <- read_csv(file.path(met_dir, "drainage_area.csv"), show_col_types = FALSE) %>%
  mutate(
    SITECODE = recode(
      SITECODE,
      "GSWSMC" = "GSMACK",
      "GSLOOK" = "GSLOOK_FULL"
    )
  )

# 2) Define which sites go into the composite
gslook_components <- c("GSWS01", "GSWS06", "LONGER", "COLD")

# 3) Build GSLOOK_FULL composite **without** Q_mm_d
gslook_full_df <- all_watersheds_data %>%
  filter(SITECODE %in% gslook_components) %>%
  group_by(DATE) %>%
  summarise(
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(SITECODE = "GSLOOK_FULL") %>%
  select(DATE, SITECODE, everything())

# 4) Compute GSLOOK_FULL discharge from the raw "GSLOOK" records:
gslook_q_full <- discharge %>%
  filter(SITECODE == "GSLOOK") %>%
  left_join(da_df, by = "SITECODE") %>%
  mutate(
    # convert m3/s to mm/day: 1 m3/s = 0.0283168 m3/s = mm/day over DA_M2
    Q_mm_d   = MEAN_Q * 0.0283168 * 86400 / DA_M2,
    SITECODE = "GSLOOK_FULL"
  ) %>%
  select(DATE, SITECODE, Q_mm_d)

# 5) Attach discharge onto the composite
gslook_full_df <- gslook_full_df %>%
  left_join(gslook_q_full, by = c("DATE", "SITECODE"))

# 6) Re-assemble your master table
all_watersheds_data <- bind_rows(
  all_watersheds_data      %>% filter(SITECODE != "GSLOOK_FULL"),
  gslook_full_df
)

# 7) And keep your list of per-site tibbles in sync
watershed_datasets[["GSLOOK_FULL"]] <- gslook_full_df

all_watersheds_data <- bind_rows(
  all_watersheds_data      %>% filter(SITECODE != "GSLOOK_FULL"),
  gslook_full_df
) %>%
  # merge the two Q_mm_d variants into one
  mutate(
    Q_mm_d = coalesce(Q_mm_d.x, Q_mm_d)      
  ) %>%
  # drop the leftover suffixed column(s)
  select(-Q_mm_d.x, -Q_mm_d.y)

# Write the watersheds data to files (NOT the met station data)
write_csv(all_watersheds_data, file.path(output_dir, "data", "watersheds_met_data_q.csv"))
