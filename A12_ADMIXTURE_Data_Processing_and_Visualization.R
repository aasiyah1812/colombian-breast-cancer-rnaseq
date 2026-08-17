# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(stringr)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. DEFINE EXACT PATHS 
# -----------------------------------------------------------------------------
fam_file     <- "Data/Admixture/FINAL_Ancestry_Data.fam"
q_file       <- "Data/Admixture/FINAL_Ancestry_Data.3.Q"
output_table <- "Data/Metadata/RNA_Ancestry_Proportions.csv"
output_plot  <- "Results/Figures/RNA_Admixture_K3.png"

# -----------------------------------------------------------------------------
# 2. LOAD THE DATA
# -----------------------------------------------------------------------------
fam_data <- read_table(fam_file, col_names = c("FamilyID", "SampleID", "PaternalID", "MaternalID", "Sex", "Phenotype"))
q_data <- read_table(q_file, col_names = c("K1", "K2", "K3"))

admixture_raw <- bind_cols(fam_data %>% select(SampleID), q_data)

# -----------------------------------------------------------------------------
# 3. MAP COLUMNS (K1=AFR, K2=EUR, K3=NAM)
# -----------------------------------------------------------------------------
admixture_clean <- admixture_raw %>%
  rename(
    AFR = K1,
    EUR = K2,
    NAM = K3
  )

# -----------------------------------------------------------------------------
# 4. ISOLATE, CLEAN, AND FORMAT THE COHORT
# -----------------------------------------------------------------------------
clean_sample_id <- function(raw_id) {
  if (str_detect(raw_id, "B138")) return("138B")
  if (str_detect(raw_id, "B205")) return("205B")
  temp_id <- str_remove_all(raw_id, "[[:punct:]]")
  core_id <- str_extract(temp_id, "^[F0-9]+B")
  return(core_id)
}

cohort_data <- admixture_clean %>%
  mutate(Patient_ID = map_chr(SampleID, clean_sample_id)) %>%
  filter(!is.na(Patient_ID)) %>%
  select(Patient_ID, AFR, EUR, NAM)

print(paste("Total patients isolated:", nrow(cohort_data)))

# -----------------------------------------------------------------------------
# 5. EXPORT THE RNA ANCESTRY TABLE FOR METADATA INTEGRATION
# -----------------------------------------------------------------------------
write_csv(cohort_data, output_table)
print(paste("RNA Ancestry table successfully saved to:", output_table))

# -----------------------------------------------------------------------------
# 6. GENERATE ADMIXTURE PLOT
# -----------------------------------------------------------------------------
plot_data <- cohort_data %>%
  pivot_longer(cols = c(AFR, EUR, NAM), names_to = "Ancestry", values_to = "Proportion")

# Explicitly set factor levels so ggplot builds them in the correct order
plot_data$Ancestry <- factor(plot_data$Ancestry, levels = c("EUR", "NAM", "AFR"))

# Standardized color palette
thesis_palette <- c("EUR" = "#0072B2", 
                    "NAM" = "#009E73", 
                    "AFR" = "#E69F00") 

admixture_plot <- ggplot(plot_data, aes(x = reorder(Patient_ID, Proportion), y = Proportion, fill = Ancestry)) +
  geom_col(width = 1, color = "black", linewidth = 0.1) + 
  scale_fill_manual(
    values = thesis_palette,
    breaks = c("EUR", "NAM", "AFR"),
    labels = c("European", "Native American", "African")
  ) +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(
    title = "Cohort Ancestry Proportions (K = 3)",
    x = "Patient ID", 
    y = "Proportion of Inferred Ancestry", 
    fill = "Ancestral Component"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave(output_plot, plot = admixture_plot, width = 12, height = 5, dpi = 300)
print(paste("Plot successfully saved to:", output_plot))
