# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# -----------------------------------------------------------------------------
# 1. LOAD DATA & RE-APPLY QUARTILES
# -----------------------------------------------------------------------------
meta <- read_excel("Data/Metadata/BaselineSamples_metadata_ancestry.xlsx") %>%
  mutate(across(ends_with(c("_DNA", "_RNA")), as.numeric))

q1_afr <- quantile(meta$AFR_DNA, 0.25, na.rm = TRUE)
q3_afr <- quantile(meta$AFR_DNA, 0.75, na.rm = TRUE)

meta_clinical <- meta %>%
  mutate(
    AFR_Baseline_Group = case_when(
      AFR_DNA <= q1_afr ~ "Low Baseline AFR",
      AFR_DNA >= q3_afr ~ "High Baseline AFR",
      TRUE ~ "Intermediate" 
    )
  ) %>%
  filter(AFR_Baseline_Group != "Intermediate") %>%
  mutate(AFR_Baseline_Group = factor(AFR_Baseline_Group, levels = c("Low Baseline AFR", "High Baseline AFR")))

# -----------------------------------------------------------------------------
# 2. STATISTICAL TESTS (Fisher's Exact Test)
# -----------------------------------------------------------------------------
response_table <- table(meta_clinical$AFR_Baseline_Group, meta_clinical$Response)
fisher_response <- fisher.test(response_table)

subtype_table <- table(meta_clinical$AFR_Baseline_Group, meta_clinical$Subtype)
fisher_subtype <- fisher.test(subtype_table, simulate.p.value = TRUE)

cat("\n--- CLINICAL STATISTICS ---\n")
cat("P-value for Ancestry vs. Therapy Response:", signif(fisher_response$p.value, 3), "\n")
cat("P-value for Ancestry vs. Clinical Subtype:", signif(fisher_subtype$p.value, 3), "\n\n")

# -----------------------------------------------------------------------------
# 3. GENERATE CLINICAL PLOTS
# -----------------------------------------------------------------------------
response_colors <- c("R" = "#009E73", "NR" = "#0072B2") 

plot_response <- ggplot(meta_clinical, aes(x = AFR_Baseline_Group, fill = Response)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.3) + 
  scale_fill_manual(values = response_colors, labels = c("NR" = "Non-Responder", "R" = "Responder")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Therapy Response by Ancestry Cohort",
    x = "Genomic Cohort",
    y = "Percentage of Patients",
    fill = "Clinical Outcome"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "right"
  )

ggsave("Results/Figures/Clinical_Response_by_Ancestry.png", plot = plot_response, width = 6, height = 5, dpi = 300)

subtype_colors <- c("TNBC" = "#D55E00", "LumA" = "#E69F00", 
                    "LumBHer2-" = "#56B4E9", "LumBHer2+" = "#0072B2", "HER2" = "#CC79A7")

plot_subtype <- ggplot(meta_clinical, aes(x = AFR_Baseline_Group, fill = Subtype)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = subtype_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Breast Cancer Subtype Distribution",
    x = "Genomic Cohort",
    y = "Percentage of Patients",
    fill = "Clinical Subtype"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "right"
  )

ggsave("Results/Figures/Clinical_Subtype_by_Ancestry.png", plot = plot_subtype, width = 6.5, height = 5, dpi = 300)

print("Clinical plots successfully generated!")
