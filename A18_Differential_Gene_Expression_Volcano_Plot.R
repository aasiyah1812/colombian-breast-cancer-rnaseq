# -----------------------------------------------------------------------------
# VOLCANO PLOT: HIGH VS LOW AFRICAN ANCESTRY
# -----------------------------------------------------------------------------
library(tidyverse)
library(ggrepel) 

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# Load the annotated results utilizing the European CSV format delimiter
res_df <- read_csv2("Results/Tables/DESeq2_High_vs_Low_AFR_Annotated.csv")

# Set strict biological and statistical thresholds
p_threshold <- 0.05
lfc_threshold <- 1

# Categorize genes based on established thresholds for differential expression
volcano_data <- res_df %>%
  mutate(
    Significance = case_when(
      padj < p_threshold & log2FoldChange > lfc_threshold ~ "Upregulated in High AFR",
      padj < p_threshold & log2FoldChange < -lfc_threshold ~ "Downregulated in High AFR",
      TRUE ~ "Not Significant"
    )
  )

# Extract the top 10 significant genes based on adjusted p-value for plot labeling
top_genes <- volcano_data %>%
  filter(Significance != "Not Significant", !is.na(Gene_Symbol)) %>%
  arrange(padj) %>%
  slice_head(n = 10)

# Define color mapping for the specific regulatory states
thesis_palette <- c("Upregulated in High AFR" = "#E69F00", 
                    "Downregulated in High AFR" = "#0072B2", 
                    "Not Significant" = "gray80")

# Generate the Volcano Plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = thesis_palette) +
  
  # Superimpose reference lines for the applied LFC and p-value thresholds
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Apply text annotations for the top significant genes while preventing visual overlap
  geom_text_repel(data = top_genes, aes(label = Gene_Symbol), color = "black", size = 3.5, fontface = "italic", max.overlaps = 15) +
  
  labs(
    title = "Differential Gene Expression: High vs. Low Baseline African Ancestry",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Export the plot to the target directory
ggsave("Results/Figures/Volcano_AFR_Annotated.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
print("Annotated Volcano Plot successfully saved!")
