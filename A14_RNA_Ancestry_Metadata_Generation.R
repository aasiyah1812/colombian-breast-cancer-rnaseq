# -----------------------------------------------------------------------------
# RNA INFERRED ANCESTRY METADATA GENERATION SCRIPT (K=3, UNIQUE MAPPING)
# -----------------------------------------------------------------------------
library(tidyverse)
library(stringr)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# 1. Load data
fam_data <- read.table("Data/Admixture/FINAL_Ancestry_Data.fam", header = FALSE, stringsAsFactors = FALSE)
colnames(fam_data) <- c("FID", "Sample_ID", "PID", "MID", "Sex", "Pheno")

# Load Q data for K=3 with generic column names (K1, K2, K3)
q_data <- read.table("Data/Admixture/FINAL_Ancestry_Data.3.Q", header = FALSE)
colnames(q_data) <- paste0("K", 1:ncol(q_data))

metadata <- read.delim("Data/Reference_1000G/igsr_samples.tsv", stringsAsFactors = FALSE)

# 2. Combine and merge with IGSR metadata to identify reference samples
admixture_raw <- bind_cols(fam_data %>% select(Sample_ID), q_data)
admixture_labeled <- merge(admixture_raw, metadata, by.x = "Sample_ID", by.y = "Sample.name", all.x = TRUE)

# -----------------------------------------------------------------------------
# 3. AUTOMATICALLY DETERMINE UNIQUE ANCESTRY MAPPING USING 1000G REFERENCES
# -----------------------------------------------------------------------------
# Calculate the average proportion of each K column across 1000G superpopulations
ref_means <- admixture_labeled %>%
  filter(Superpopulation.code %in% c("AFR", "EUR", "AMR")) %>%
  group_by(Superpopulation.code) %>%
  summarise(across(starts_with("K"), mean))

cat("--- 1000G Reference Mean Proportions (Diagnostic Check) ---\n")
print(ref_means)

# Ensure strict 1-to-1 unique mapping for K=3
# Step A: AFR is the column peaking in the AFR reference
afr_col <- ref_means %>% filter(Superpopulation.code == "AFR") %>% select(starts_with("K")) %>% pivot_longer(everything()) %>% slice_max(value, n = 1) %>% pull(name)

# Step B: EUR is the column peaking in the EUR reference (excluding AFR column)
eur_col <- ref_means %>% filter(Superpopulation.code == "EUR") %>% select(starts_with("K")) %>% select(-all_of(afr_col)) %>% pivot_longer(everything()) %>% slice_max(value, n = 1) %>% pull(name)
if(length(eur_col) == 0) {
  eur_col <- ref_means %>% filter(Superpopulation.code == "EUR") %>% select(starts_with("K")) %>% pivot_longer(everything()) %>% slice_max(value, n = 1) %>% pull(name)
}

# Step C: NAM (Native American proxy via AMR reference) is the remaining unassigned column
nam_col <- setdiff(paste0("K", 1:3), c(afr_col, eur_col))

detected_mapping <- c(AFR = afr_col, EUR = eur_col, NAM = nam_col)
cat("\n--- Automated Unique Ancestry Column Mapping ---\n")
print(detected_mapping)

# Extract column names safely as plain strings
afr_col <- unname(detected_mapping["AFR"])
eur_col <- unname(detected_mapping["EUR"])
nam_col <- unname(detected_mapping["NAM"])

# -----------------------------------------------------------------------------
# 4. APPLY MAPPING AND CLEAN COHORT IDs
# -----------------------------------------------------------------------------
clean_sample_id <- function(raw_id) {
  if (str_detect(raw_id, "B138")) return("138B")
  if (str_detect(raw_id, "B205")) return("205B")
  temp_id <- str_remove_all(raw_id, "[[:punct:]]")
  core_id <- str_extract(temp_id, "^[F0-9]+B")
  return(core_id)
}

rna_metadata <- admixture_labeled %>%
  mutate(
    Patient_ID = map_chr(Sample_ID, clean_sample_id),
    Dataset_Source = if_else(is.na(Superpopulation.code), "Colombian_Cohort", "1000G_Reference"),
    # Map dynamically using unique mapped columns
    AFR = .data[[afr_col]],
    EUR = .data[[eur_col]],
    NAM = .data[[nam_col]]
  ) %>%
  select(Sample_ID, Patient_ID, Dataset_Source, Superpopulation.code, Population.code, AFR, EUR, NAM)

# 5. Save output table
output_metadata_path <- "Data/Metadata/RNA_Inferred_Ancestry_Metadata.csv"
write_csv(rna_metadata, output_metadata_path)
cat("\nRNA Ancestry metadata table successfully saved to:", output_metadata_path, "\n")
