# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. LOAD & PREPARE DATA
# -----------------------------------------------------------------------------
# Import the consolidated baseline metadata (incorporating corrected RNA ancestry data)
data <- read_excel("Data/Metadata/BaselineSamples_metadata_ancestry.xlsx")

# Ensure ancestry fraction columns are strictly numeric
data <- data %>%
  mutate(across(ends_with(c("_DNA", "_RNA")), as.numeric))

# Calculate the difference (Delta) between the transcribed RNA and germline DNA
data_delta <- data %>%
  mutate(
    EUR_Delta = EUR_RNA - EUR_DNA,
    NAM_Delta = NAM_RNA - NAM_DNA,
    AFR_Delta = AFR_RNA - AFR_DNA
  )

# -----------------------------------------------------------------------------
# 2. STATISTICAL ANALYSIS (PAIRED WILCOXON SIGNED-RANK TEST)
# -----------------------------------------------------------------------------
run_ancestry_stats <- function(ancestry_name) {
  col_dna <- paste0(ancestry_name, "_DNA")
  col_rna <- paste0(ancestry_name, "_RNA")
  
  wt <- wilcox.test(data[[col_rna]], data[[col_dna]], paired = TRUE)
  
  tibble(
    `Ancestral Lineage` = ancestry_name,
    `Mean Germline (DNA)` = mean(data[[col_dna]], na.rm = TRUE),
    `Mean Transcriptomic (RNA)` = mean(data[[col_rna]], na.rm = TRUE),
    `Mean Shift` = mean(data_delta[[paste0(ancestry_name, "_Delta")]], na.rm = TRUE),
    `p-value` = signif(wt$p.value, digits = 3)
  )
}

stats_results <- bind_rows(
  run_ancestry_stats("EUR"),
  run_ancestry_stats("NAM"),
  run_ancestry_stats("AFR")
)

# Export the final statistical table as a CSV
write_csv(stats_results, "Results/Tables/Ancestry_Divergence_Statistics.csv")
print("--- Statistical Proof of Ancestry Shifts ---")
print(stats_results)

# -----------------------------------------------------------------------------
# 3. VISUALIZATION OF THE SHIFTS
# -----------------------------------------------------------------------------
plot_data <- data_delta %>%
  select(Patient_ID, ends_with("_Delta")) %>%
  pivot_longer(cols = ends_with("_Delta"), names_to = "Lineage", values_to = "Delta") %>%
  mutate(Lineage = str_remove(Lineage, "_Delta"))

# Enforce strict factor levels matching previous visualizations
plot_data$Lineage <- factor(plot_data$Lineage, levels = c("EUR", "NAM", "AFR"))

# Color palette mapping
thesis_palette <- c(
  "AFR" = "#0072B2", 
  "EUR" = "#009E73", 
  "NAM" = "#E69F00"
)

delta_plot <- ggplot(plot_data, aes(x = Lineage, y = Delta, fill = Lineage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  scale_fill_manual(values = thesis_palette) +
  labs(
    title = "Divergence of Transcriptomic Ancestry from Germline Baselines",
    x = "Ancestral Lineage",
    y = "Shift in Proportion (RNA - DNA)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "none"
  )

ggsave("Results/Figures/Ancestry_Delta_Boxplots.png", plot = delta_plot, width = 7, height = 6, dpi = 300)
print("Divergence boxplot successfully generated!")
