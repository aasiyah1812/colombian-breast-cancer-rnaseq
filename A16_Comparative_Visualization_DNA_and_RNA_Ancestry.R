# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(stringr)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. LOAD DATA & FIX TYPES
# -----------------------------------------------------------------------------
# Define path and load the baseline metadata containing both DNA and RNA ancestry estimates
data_path <- "Data/Metadata/BaselineSamples_metadata_ancestry.xlsx"
data <- read_excel(data_path)

# Ensure columns ending in "_DNA" or "_RNA" are formatted strictly as numeric variables
data <- data %>%
  mutate(across(ends_with(c("_DNA", "_RNA")), as.numeric))

# -----------------------------------------------------------------------------
# 2. RESHAPE DATA WITH NEUTRAL LABELS
# -----------------------------------------------------------------------------
# Isolate and reshape the DNA baseline dataset
dna_long <- data %>% 
  select(Patient_ID, Response, ends_with("_DNA")) %>%
  pivot_longer(cols = ends_with("_DNA"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Ancestry = str_remove(Ancestry, "_DNA"), 
         Type = "Germline (DNA)") 

# Isolate and reshape the RNA transcriptomic dataset
rna_long <- data %>% 
  select(Patient_ID, Response, ends_with("_RNA")) %>%
  pivot_longer(cols = ends_with("_RNA"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Ancestry = str_remove(Ancestry, "_RNA"), 
         Type = "Transcriptomic (RNA)") 

# Extract a strict vector of Patient_IDs sorted by their baseline European DNA ancestry
patient_order <- data %>% arrange(desc(EUR_DNA)) %>% pull(Patient_ID)

# Bind the reshaped DNA and RNA datasets vertically
plot_data <- bind_rows(dna_long, rna_long) %>%
  mutate(
    Patient_ID = factor(Patient_ID, levels = patient_order),
    Type = factor(Type, levels = c("Transcriptomic (RNA)", "Germline (DNA)"))
  )

# -----------------------------------------------------------------------------
# 3. GENERATE COMPARISON PLOT
# -----------------------------------------------------------------------------
# Standardized color palette
thesis_palette <- c(
  "AFR" = "#0072B2", 
  "EUR" = "#009E73", 
  "NAM" = "#E69F00"
)

# Initialize ggplot object mapping the ordered patient IDs to the respective ancestry proportions
comparison_plot <- ggplot(plot_data, aes(x = Patient_ID, y = Proportion, fill = Ancestry)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = thesis_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  
  # Stratify the plot horizontally by biological source (Type) and vertically by clinical outcome (Response)
  facet_grid(Type ~ Response, scales = "free_x", space = "free_x") +
  
  # Apply axis labels and title with consistent styling
  labs(
    title = "Comparison of Germline (DNA) and Transcriptomic (RNA) Ancestry",
    x = "Patient ID", 
    y = "Proportion of Genetic Ancestry", 
    fill = "Ancestral Lineage"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"), 
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines")
  )

# Export the plot to the results directory at 300 dpi
ggsave("Results/Figures/Ancestry_DNA_vs_RNA_Comparison.png", plot = comparison_plot, width = 12, height = 7, dpi = 300)
print("Comparison plot successfully generated with updated colors and matching styling!")
