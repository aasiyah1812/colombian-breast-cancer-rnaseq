# -----------------------------------------------------------------------------
# HEATMAP: TOP 50 DIFFERENTIALLY EXPRESSED GENES
# -----------------------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(pheatmap)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# Load the annotated DESeq2 output
res_df <- read_csv2("Results/Tables/DESeq2_High_vs_Low_AFR_Annotated.csv")

# Isolate the top 50 significant genes strictly by adjusted p-value
top_genes_df <- res_df %>%
  filter(padj < 0.05, !is.na(Gene_Symbol)) %>%
  arrange(padj) %>%
  slice_head(n = 50)

# Extract raw Ensembl IDs to query the DESeq2 object, and symbols for plotting
top_ensembl <- top_genes_df$Gene_ID
top_symbols <- top_genes_df$Gene_Symbol

# Apply a Variance Stabilizing Transformation (VST) to the DESeq2 counts
vsd <- vst(dds, blind = FALSE)

# Subset the transformed matrix to include only the top 50 significant identifiers
mat <- assay(vsd)[top_ensembl, ]

# Replace Ensembl IDs with human-readable Gene Symbols for the final visualization
rownames(mat) <- top_symbols

# Standardize the rows (z-score scaling) to allow relative expression comparison
mat_scaled <- t(scale(t(mat)))

# Construct the clinical annotation framework for the plot header
df_anno <- as.data.frame(colData(dds)[, c("AFR_Baseline_Group", "Response")])
colnames(df_anno) <- c("Baseline Ancestry", "Therapy Response")

# Define visual variables for the clinical annotations
anno_colors <- list(
  `Baseline Ancestry` = c(High_AFR = "#E69F00", Low_AFR = "gray80"),
  `Therapy Response` = c(R = "#009E73", NR = "#0072B2")
)

# Render and export the Heatmap
png("Results/Figures/Heatmap_Top50_AFR_Annotated.png", width = 8, height = 10, units = "in", res = 300)
pheatmap(mat_scaled, 
         annotation_col = df_anno, 
         annotation_colors = anno_colors,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_colnames = FALSE, 
         fontsize_row = 8,
         fontface_row = "italic",
         main = "Expression Profile of Top 50 DGEs (High vs Low Baseline AFR)")
dev.off()
print("Annotated Heatmap successfully saved!")
