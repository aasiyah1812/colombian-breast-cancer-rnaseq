# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. LOAD & PREPARE RANKED GENE LIST
# -----------------------------------------------------------------------------
res_df <- read_csv2("Results/Tables/DESeq2_High_vs_Low_AFR_Annotated.csv")

ranked_data <- res_df %>%
  filter(!is.na(Gene_Symbol) & !is.na(log2FoldChange)) %>%
  distinct(Gene_Symbol, .keep_all = TRUE) %>%
  arrange(desc(log2FoldChange))

gene_list <- ranked_data$log2FoldChange
names(gene_list) <- ranked_data$Gene_Symbol

# -----------------------------------------------------------------------------
# 2. RUN GSEA (Gene Set Enrichment Analysis)
# -----------------------------------------------------------------------------
print("Running GSEA across the entire ranked genome...")

gsea_results <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP", 
  minGSSize    = 10,  
  maxGSSize    = 500,  
  pvalueCutoff = 0.05, 
  verbose      = FALSE
)

# -----------------------------------------------------------------------------
# 3. GENERATE THE GSEA DOTPLOT
# -----------------------------------------------------------------------------
if (nrow(as.data.frame(gsea_results)) > 0) {
  gsea_plot <- dotplot(gsea_results, showCategory = 10, split = ".sign", title = "GSEA: Systemic Pathway Shifts in High AFR Tumors") +
    facet_grid(.~.sign) + 
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 12) 
    )
  
  ggsave("Results/Figures/GSEA_Pathway_Shifts.png", plot = gsea_plot, width = 12, height = 7, dpi = 300)
  print("GSEA Plot successfully saved!")
  
} else {
  print("Diagnostic Note: GSEA identified no statistically significant pathway level shifts.")
}
