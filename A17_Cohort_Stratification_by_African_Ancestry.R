# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. LOAD DATA & FIX TYPES
# -----------------------------------------------------------------------------
# Import the baseline metadata
meta <- read_excel("Data/Metadata/BaselineSamples_metadata_ancestry.xlsx") %>%
  # Ensure all variables representing DNA and RNA ancestry proportions are strictly numeric
  mutate(across(ends_with(c("_DNA", "_RNA")), as.numeric))

# -----------------------------------------------------------------------------
# 2. CALCULATE MEDIAN AND STRATIFY COHORT
# -----------------------------------------------------------------------------
# Compute the statistical median of the inherited African (AFR) DNA ancestry across the cohort
median_afr <- median(meta$AFR_DNA, na.rm = TRUE)

meta_grouped <- meta %>%
  # Assign a binary categorical variable: samples at or above the median are designated "High_AFR", 
  # while those below the median are designated "Low_AFR"
  mutate(AFR_Baseline_Group = if_else(AFR_DNA >= median_afr, "High_AFR", "Low_AFR")) %>%
  # Retain only the variables required for downstream comparative analysis
  select(Patient_ID, Response, AFR_DNA, AFR_Baseline_Group) 

# Export the stratified dataset using European Excel formatting (semicolon delimited)
write_excel_csv2(meta_grouped, "Results/Tables/Metadata_AFR_Stratified.csv")
print("Grouped metadata saved successfully!")
