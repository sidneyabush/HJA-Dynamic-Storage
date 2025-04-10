# Clear environment
rm(list = ls())

# Set working directory (adjust as needed)
setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem/raw_data")

# Load needed libraries
librarian::shelf(dplyr, googledrive, ggplot2, data.table, lubridate, tidyr, stringr, readr)

# Read in your data and filter streams
RBFI <- read.csv("flashiness_by_stream_id.csv") %>%
  filter(startsWith(Stream_ID, "AND"))
recession <- read.csv("Recession_Slopes_by_StreamID_Aggregate.csv") %>%
  filter(startsWith(Stream_ID, "AND"))

# Join the datasets
FI_recession <- left_join(RBFI, recession, by = "Stream_ID")

# Reshape the data to long format (creates "Metric" column, not "Parameter")
FI_long <- FI_recession %>%
  pivot_longer(
    cols = c(RBFI, slope),
    names_to = "Metric",
    values_to = "Value"
  )

# Define unique Stream_ID and lookup table for colors
unique_streams <- c("AND__GSLOOK",
                    "AND__GSMACK",
                    "AND__GSWS01",
                    "AND__GSWS02",
                    "AND__GSWS06",
                    "AND__GSWS07",
                    "AND__GSWS08",
                    "AND__GSWS09",
                    "AND__GSWS10")

color_lookup <- data.frame(
  Stream_ID = unique_streams,
  RBFI_color = c("lightgray", "#d0abcc", "#dfb9bd", "#f1beab", "#fce8b6", "#a2bea8", "#aac9c3", "#bdd5e5", "#a3a4c3"),
  slope_color = c("gray30", "#963f89", "#bb616d", "#e16b43", "#f7c955", "#2d6b37", "#509e8e", "#88c2e8", "#2d2d75"),
  stringsAsFactors = FALSE
)

# Merge by Stream_ID; ensure FI_long has a matching column "Stream_ID"
FI_long <- merge(FI_long, color_lookup, by = "Stream_ID", all.x = TRUE)

# Create fill_color column based on Metric
FI_long$fill_color <- ifelse(FI_long$Metric == "RBFI", FI_long$RBFI_color, FI_long$slope_color)

ggplot(FI_long, aes(x = Stream_ID, y = Value, fill = fill_color)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_identity() +
  labs(title = "Recession Slope = Dark; Flashiness Index = Light",
       x = NULL,
       y = "average value (2001-2019)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

