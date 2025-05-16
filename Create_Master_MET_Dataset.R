library(readr)
library(dplyr)
library(lubridate)

# ── Clear environment ─────────────────────────────────────────────────────────
rm(list = ls())

# ── Directories & site list ─────────────────────────────────────────────────
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"

# Need to import: 
# Precip
# Net Radiation
# Air Temp
# RH

# The databases and dates will vary by watershed
# GSLOOK (the entirety of the watershed) will need to be averaged or interpolated

# Need to decide a method for calculating ET
# Latitude method?
# Calculate Ground Heat Flux