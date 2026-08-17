# -----------------------------------------------------------------------------
# DESEQ2: TNBC VS NON-TNBC SUBTYPE SIGNATURE
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# Load full metadata (all 58 patients)
meta <- read_excel("Data/Metadata/BaselineSamples_metadata_ancestry.xlsx") %>%
  mutate(
    TNBC_Status = if_else(Subtype == "TNBC", "TNBC", "Non_TNBC"),
    TNBC_Status = factor(TNBC_Status, levels = c("Non_TNBC", "TNBC"))
  )

# Load and clean count matrix headers
counts_raw <- read.delim("Data/RNAseq/Master_Expression_Matrix.txt", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
colnames(counts_raw) <- sapply(colnames(counts_raw), function(x) {
  if (str_detect(x, "B138")) return("138B")
  if (str_detect(x, "B205")) return("205B")
  return(str_extract(x, "^[F0-9]+B"))
})

# Align and check samples
common_samples <- intersect(meta$Patient_ID, colnames(counts_raw))
counts_filtered <- counts_raw[, common_samples]
meta_filtered <- meta %>% filter(Patient_ID %in% common_samples) %>% arrange(match(Patient_ID, common_samples))

# Run DESeq2 Pipeline
dds_tnbc <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta_filtered, design = ~ TNBC_Status)
keep <- rowSums(counts(dds_tnbc) >= 10) >= 3
dds_tnbc <- dds_tnbc[keep, ]
dds_tnbc <- DESeq(dds_tnbc)

# Extract and annotate results
res_tnbc <- as.data.frame(results(dds_tnbc, contrast = c("TNBC_Status", "TNBC", "Non_TNBC"))) %>%
  rownames_to_column(var = "Gene_ID") %>%
  mutate(Clean_Ensembl = str_remove(Gene_ID, "\\..*"))

res_tnbc$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = res_tnbc$Clean_Ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Save the full annotated TNBC results table
write_excel_csv2(res_tnbc, "Results/Tables/DESeq2_TNBC_vs_NonTNBC_Annotated.csv")
print("TNBC DESeq2 complete and table saved successfully!")
