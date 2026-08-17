# -----------------------------------------------------------------------------
# SUPPLEMENTARY ADMIXTURE PLOTS (K = 4, 5, 6)
# COHORT & GLOBAL VISUALIZATIONS
# -----------------------------------------------------------------------------
library(tidyverse)
library(stringr)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# Define K values for supplementary analysis
k_values <- c(4, 5, 6)

# Extended colorblind-friendly palette for up to 6 components
extended_palette <- c(
  "K1" = "#0072B2", # Blue
  "K2" = "#009E73", # Bluish green
  "K3" = "#E69F00", # Orange
  "K4" = "#CC79A7", # Reddish purple
  "K5" = "#56B4E9", # Sky blue
  "K6" = "#F0E442"  # Yellow
)

# =============================================================================
# 1. COHORT ADMIXTURE PLOTS (K = 4, 5, 6)
# =============================================================================
for (k in k_values) {
  cat("Processing Cohort Admixture for K =", k, "\n")
  
  fam_file     <- "Data/Admixture/FINAL_Ancestry_Data.fam"
  q_file       <- paste0("Data/Admixture/FINAL_Ancestry_Data.", k, ".Q")
  output_table <- paste0("Data/Metadata/RNA_Ancestry_Proportions_K", k, ".csv")
  output_plot  <- paste0("Results/Figures/RNA_Admixture_K", k, ".png")
  
  # Load data
  fam_data <- read_table(fam_file, col_names = c("FamilyID", "SampleID", "PaternalID", "MaternalID", "Sex", "Phenotype"), show_col_types = FALSE)
  q_data <- read_table(q_file, col_names = paste0("K", 1:k), show_col_types = FALSE)
  
  admixture_raw <- bind_cols(fam_data %>% select(SampleID), q_data)
  
  # Clean sample ID function
  clean_sample_id <- function(raw_id) {
    if (str_detect(raw_id, "B138")) return("138B")
    if (str_detect(raw_id, "B205")) return("205B")
    temp_id <- str_remove_all(raw_id, "[[:punct:]]")
    core_id <- str_extract(temp_id, "^[F0-9]+B")
    return(core_id)
  }
  
  cohort_data <- admixture_raw %>%
    mutate(Patient_ID = map_chr(SampleID, clean_sample_id)) %>%
    filter(!is.na(Patient_ID)) %>%
    select(Patient_ID, all_of(paste0("K", 1:k)))
  
  # Export table for metadata integration
  write_csv(cohort_data, output_table)
  
  # Prepare plot data
  plot_data <- cohort_data %>%
    pivot_longer(cols = all_of(paste0("K", 1:k)), names_to = "Ancestry", values_to = "Proportion")
  
  plot_data$Ancestry <- factor(plot_data$Ancestry, levels = paste0("K", 1:k))
  
  current_palette <- extended_palette[1:k]
  names(current_palette) <- paste0("K", 1:k)
  
  admixture_plot <- ggplot(plot_data, aes(x = reorder(Patient_ID, Proportion), y = Proportion, fill = Ancestry)) +
    geom_col(width = 1, color = "black", linewidth = 0.1) + 
    scale_fill_manual(
      values = current_palette,
      labels = paste0("Population ", 1:k)
    ) +
    scale_y_continuous(expand = c(0, 0)) + 
    labs(
      title = paste0("Cohort Ancestry Proportions (K = ", k, ")"),
      x = "Patient ID", 
      y = "Proportion of Inferred Ancestry", 
      fill = "Ancestral Component"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  ggsave(output_plot, plot = admixture_plot, width = 12, height = 5, dpi = 300)
}

# =============================================================================
# 2. GLOBAL ADMIXTURE PLOTS (K = 4, 5, 6)
# =============================================================================
metadata <- read.delim("Data/Reference_1000G/igsr_samples.tsv", stringsAsFactors = FALSE)

for (k in k_values) {
  cat("Processing Global Admixture for K =", k, "\n")
  
  fam_file    <- "Data/Admixture/FINAL_Ancestry_Data.fam"
  q_file      <- paste0("Data/Admixture/FINAL_Ancestry_Data.", k, ".Q")
  output_plot <- paste0("Results/Figures/Global_Admixture_K", k, "_Clean.png")
  
  fam_data <- read.table(fam_file, header = FALSE)
  colnames(fam_data) <- c("FID", "Sample_ID", "PID", "MID", "Sex", "Pheno")
  
  q_data <- read.table(q_file, header = FALSE)
  colnames(q_data) <- paste0("K", 1:k)
  
  admixture_raw <- bind_cols(fam_data %>% select(Sample_ID), q_data)
  admixture_labeled <- merge(admixture_raw, metadata, by.x = "Sample_ID", by.y = "Sample.name", all.x = TRUE)
  
  admixture_clean <- admixture_labeled %>%
    mutate(Plot_Group = case_when(
      is.na(Superpopulation.code) ~ "Cohort",
      TRUE ~ Superpopulation.code 
    )) %>%
    pivot_longer(cols = all_of(paste0("K", 1:k)), names_to = "Ancestry", values_to = "Proportion") %>%
    filter(!is.na(Plot_Group))
  
  admixture_clean$Plot_Group <- factor(admixture_clean$Plot_Group, 
                                       levels = c("AFR", "AMR", "Cohort", "EAS", "EUR", "SAS"))
  
  current_palette <- extended_palette[1:k]
  names(current_palette) <- paste0("K", 1:k)
  
  global_admixture_plot <- ggplot(admixture_clean, aes(x = Sample_ID, y = Proportion, fill = Ancestry)) +
    geom_col(width = 1, color = NA) +
    scale_fill_manual(values = current_palette) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(~Plot_Group, scales = "free_x", space = "free_x") +
    theme_classic() +
    labs(
      title = paste0("Global Admixture Reference Comparison (K = ", k, ")"),
      y = "Ancestry Proportion",
      x = NULL
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 11, angle = 90, hjust = 0),
      panel.spacing = unit(0.2, "lines"),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  
  ggsave(output_plot, plot = global_admixture_plot, width = 12, height = 5, dpi = 300)
}
