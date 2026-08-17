# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. LOAD & PREPARE DATA
# -----------------------------------------------------------------------------
# Import the annotated DESeq2 results matrix
res_df <- read_csv2("Results/Tables/DESeq2_High_vs_Low_AFR_Annotated.csv")

# Filter for statistically significant genes containing valid HUGO symbols
sig_genes <- res_df %>% 
  filter(padj < 0.05, !is.na(Gene_Symbol))

# Stratify into distinct upregulated and downregulated vectors
upregulated <- sig_genes %>% 
  filter(log2FoldChange > 0) %>% 
  pull(Gene_Symbol)

downregulated <- sig_genes %>% 
  filter(log2FoldChange < 0) %>% 
  pull(Gene_Symbol)

# Define the background reference universe (all reliably detected transcripts in the assay)
universe_genes <- res_df %>% 
  filter(!is.na(Gene_Symbol)) %>% 
  pull(Gene_Symbol)

# -----------------------------------------------------------------------------
# 2. RUN GENE ONTOLOGY (GO) ENRICHMENT 
# -----------------------------------------------------------------------------
print("Running GO Enrichment for Upregulated Genes...")
go_up <- enrichGO(
  gene          = upregulated,
  universe      = universe_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,  
  qvalueCutoff  = 0.2   
)

print("Running GO Enrichment for Downregulated Genes...")
go_down <- enrichGO(
  gene          = downregulated,
  universe      = universe_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, 
  qvalueCutoff  = 0.05
)

# --- Diagnostic Check ---
cat("\nFound", nrow(as.data.frame(go_up)), "upregulated pathways.\n")
cat("Found", nrow(as.data.frame(go_down)), "downregulated pathways.\n\n")

# -----------------------------------------------------------------------------
# 3. GENERATE DOT PLOTS 
# -----------------------------------------------------------------------------
if (nrow(as.data.frame(go_up)) > 0) {
  plot_up <- dotplot(go_up, showCategory = 15, title = "Enriched Biological Processes: Upregulated in High AFR") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  ggsave("Results/Figures/GO_Enrichment_Upregulated.png", plot = plot_up, width = 9, height = 7, dpi = 300)
  print("Upregulated plot saved!")
} else {
  print("Diagnostic Note: No statistically enriched upregulated pathways met the threshold criteria.")
}

plot_down <- dotplot(go_down, showCategory = 15, title = "Enriched Biological Processes: Downregulated in High AFR") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )
ggsave("Results/Figures/GO_Enrichment_Downregulated.png", plot = plot_down, width = 9, height = 7, dpi = 300)
print("Downregulated plot saved!")
