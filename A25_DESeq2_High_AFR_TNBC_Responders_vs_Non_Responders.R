# -----------------------------------------------------------------------------
# DESEQ2: HIGH-AFR TNBC RESPONDERS VS NON-RESPONDERS
# -----------------------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)
library(ggrepel)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# 1. Load the curated 6-patient subcohort table
meta_6_patients <- read_csv2("Results/Tables/High_AFR_TNBC_Patients.csv") %>%
  mutate(
    Response = factor(Response)
  )

# 2. Load count matrix and clean headers to match patient IDs
counts_raw <- read.delim("Data/RNAseq/Master_Expression_Matrix.txt", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
colnames(counts_raw) <- sapply(colnames(counts_raw), function(x) {
  if (str_detect(x, "B138")) return("138B")
  if (str_detect(x, "B205")) return("205B")
  return(str_extract(x, "^[F0-9]+B"))
})

# 3. Match counts strictly to the 6 patients in the subcohort table
common_samples <- intersect(meta_6_patients$Patient_ID, colnames(counts_raw))
counts_filtered <- counts_raw[, common_samples]
meta_filtered <- meta_6_patients %>% 
  filter(Patient_ID %in% common_samples) %>% 
  arrange(match(Patient_ID, common_samples))

cat("\nRunning DESeq2 on exactly", nrow(meta_filtered), "High-AFR TNBC patients...\n")
print("Patient breakdown:")
print(table(meta_filtered$Response))

# 4. Run DESeq2 Pipeline
dds_sub <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta_filtered, design = ~ Response)

keep <- rowSums(counts(dds_sub) >= 10) >= 3
dds_sub <- dds_sub[keep, ]
dds_sub <- DESeq(dds_sub)

# 5. Extract and Annotate Results
res_sub <- as.data.frame(results(dds_sub)) %>%
  rownames_to_column(var = "Gene_ID") %>%
  mutate(Clean_Ensembl = str_remove(Gene_ID, "\\..*"))

res_sub$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = res_sub$Clean_Ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_sub$Full_Name <- mapIds(org.Hs.eg.db, keys = res_sub$Clean_Ensembl, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")

# 6. Save results table
final_res_sub <- res_sub %>%
  dplyr::select(Gene_Symbol, Full_Name, Gene_ID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  arrange(pvalue)

write_excel_csv2(final_res_sub, "Results/Tables/DESeq2_HighAFR_TNBC_Responders_vs_NR.csv")

# -----------------------------------------------------------------------------
# 7. GENERATE THE VOLCANO PLOT
# -----------------------------------------------------------------------------
plot_data <- final_res_sub %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(
    Significance = case_when(
      pvalue < 0.05 & log2FoldChange > 1 ~ "Upregulated in NR",
      pvalue < 0.05 & log2FoldChange < -1 ~ "Upregulated in R",
      TRUE ~ "Not Significant"
    )
  )

top_genes <- plot_data %>%
  filter(Significance != "Not Significant", !is.na(Gene_Symbol)) %>%
  arrange(pvalue) %>%
  head(10)

volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("gray80", "#0072B2", "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(data = top_genes, aes(label = Gene_Symbol), color = "black", size = 4, fontface = "bold", max.overlaps = Inf) +
  labs(
    title = "High-AFR TNBC Therapy Response Profiling",
    subtitle = paste("N =", nrow(meta_filtered), "Patients (Responders vs. Non-Responders)"),
    x = "Log2 Fold Change",
    y = "-Log10(P-Value)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave("Results/Figures/Volcano_HighAFR_TNBC_Response.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
print("Analysis complete! Data and Volcano Plot saved.")
