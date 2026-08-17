# -----------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
# Load the tidyverse suite for data manipulation and visualization
library(tidyverse)

# Set the working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 2. LOAD DATA
# -----------------------------------------------------------------------------
# Read the PLINK .eigenvec file containing the sample coordinates 
pca_data <- read.table("Data/PCA_Corrected/PCA_Results.eigenvec", header = FALSE)

# Assign column names to the first six columns to identify IDs and principal components
colnames(pca_data)[1:6] <- c("FID", "Sample_ID", "PC1", "PC2", "PC3", "PC4")

# Read the PLINK .eigenval file containing the eigenvalues 
eigenval <- read.table("Data/PCA_Corrected/PCA_Results.eigenval", header = FALSE)
colnames(eigenval) <- c("variance")

# Calculate the percentage of variance explained (PVE)
pve <- (eigenval$variance / sum(eigenval$variance)) * 100

# Read the metadata file from the 1000 Genomes Project (IGSR portal)
metadata <- read.delim("Data/Reference_1000G/igsr_samples.tsv", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 3. MERGE, CLEAN LABELS, AND SORT
# -----------------------------------------------------------------------------
# Merge the calculated PCA coordinates with the 1000 Genomes metadata
pca_labeled <- merge(pca_data, metadata, by.x = "Sample_ID", by.y = "Sample.name", all.x = TRUE)

# Create a categorical column to define plotting groups and filter out missing data
pca_labeled <- pca_labeled %>%
  mutate(Plot_Group = case_when(
    is.na(Superpopulation.code) ~ "Cohort",
    TRUE ~ Superpopulation.code 
  )) %>%
  filter(!is.na(Plot_Group)) # Removes technical artifacts generating 'NA' legends

# Convert to factor to strictly define the order, ensuring the Cohort is drawn last/on top
pca_labeled$Plot_Group <- factor(pca_labeled$Plot_Group, 
                                 levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "Cohort"))
pca_labeled <- pca_labeled %>% arrange(Plot_Group)

# Invert PC2 to align the visual orientation with standard global reference literature
pca_labeled$PC2 <- pca_labeled$PC2 * -1

# Define universal color and shape palettes
my_colors <- c("AFR" = "#E69F00", "AMR" = "#009E73", "EAS" = "#CC79A7", 
               "EUR" = "#56B4E9", "SAS" = "#D55E00", "Cohort" = "black")

my_shapes <- c("AFR" = 16, "AMR" = 16, "EAS" = 16, "EUR" = 16, "SAS" = 16, 
               "Cohort" = 17)

# -----------------------------------------------------------------------------
# 4. PLOT 1: PC1 vs PC2
# -----------------------------------------------------------------------------
plot_1v2 <- ggplot(pca_labeled, aes(x = PC1, y = PC2, color = Plot_Group, shape = Plot_Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = my_colors, na.translate = FALSE) +
  scale_shape_manual(values = my_shapes, na.translate = FALSE) +
  theme_bw() +
  labs(title = "Global PCA: PC1 vs PC2", 
       x = paste0("PC1 (", round(pve[1], 1), "%)"), 
       y = paste0("PC2 (", round(pve[2], 1), "%)"), 
       color = "Population", shape = "Population") +
  theme(legend.position = "bottom")

ggsave("Results/Figures/PCA_Annex_PC1_vs_PC2.png", plot = plot_1v2, width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 5. PLOT 2: PC1 vs PC3
# -----------------------------------------------------------------------------
plot_1v3 <- ggplot(pca_labeled, aes(x = PC1, y = PC3, color = Plot_Group, shape = Plot_Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = my_colors, na.translate = FALSE) +
  scale_shape_manual(values = my_shapes, na.translate = FALSE) +
  theme_bw() +
  labs(title = "Global PCA: PC1 vs PC3", 
       x = paste0("PC1 (", round(pve[1], 1), "%)"), 
       y = paste0("PC3 (", round(pve[3], 1), "%)"), 
       color = "Population", shape = "Population") +
  theme(legend.position = "bottom")

ggsave("Results/Figures/PCA_Annex_PC1_vs_PC3.png", plot = plot_1v3, width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 6. PLOT 3: PC1 vs PC4
# -----------------------------------------------------------------------------
plot_1v4 <- ggplot(pca_labeled, aes(x = PC1, y = PC4, color = Plot_Group, shape = Plot_Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = my_colors, na.translate = FALSE) +
  scale_shape_manual(values = my_shapes, na.translate = FALSE) +
  theme_bw() +
  labs(title = "Global PCA: PC1 vs PC4", 
       x = paste0("PC1 (", round(pve[1], 1), "%)"), 
       y = paste0("PC4 (", round(pve[4], 1), "%)"), 
       color = "Population", shape = "Population") +
  theme(legend.position = "bottom")

ggsave("Results/Figures/PCA_Annex_PC1_vs_PC4.png", plot = plot_1v4, width = 8, height = 6, dpi = 300)
