# -----------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(scales) 

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
# 2. CALCULATE EXACT PERCENTAGES FOR THE LABELS
# -----------------------------------------------------------------------------
# Summarize subset variables corresponding to neoadjuvant outcome
response_summary <- meta_clinical %>%
  group_by(AFR_Baseline_Group, Response) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(AFR_Baseline_Group) %>%
  mutate(
    Percent = Count / sum(Count),
    # Convert standard float values to formatted percentage strings utilizing the scales package
    Label = percent(Percent, accuracy = 1) 
  )

# Summarize subset variables corresponding to molecular subtype
subtype_summary <- meta_clinical %>%
  group_by(AFR_Baseline_Group, Subtype) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(AFR_Baseline_Group) %>%
  mutate(
    Percent = Count / sum(Count),
    Label = percent(Percent, accuracy = 1)
  )

# Retain standardized color profiles
response_colors <- c("R" = "#009E73", "NR" = "#0072B2") 
subtype_colors <- c("TNBC" = "#D55E00", "LumA" = "#E69F00", 
                    "LumBHer2-" = "#56B4E9", "LumBHer2+" = "#0072B2", "HER2" = "#CC79A7")

# -----------------------------------------------------------------------------
# 3. LABELED STACKED BAR CHARTS
# -----------------------------------------------------------------------------
plot_response_bar <- ggplot(response_summary, aes(x = AFR_Baseline_Group, y = Percent, fill = Response)) +
  geom_col(color = "black", linewidth = 0.3) +
  # Position text labels at the vertical center of the stacked segments
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 5) +
  scale_fill_manual(values = response_colors, labels = c("NR" = "Non-Responder", "R" = "Responder")) +
  scale_y_continuous(labels = percent) +
  labs(title = "Therapy Response by Ancestry Cohort", x = "Genomic Cohort", y = "Percentage of Patients", fill = "Clinical Outcome") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 11), axis.title = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

ggsave("Results/Figures/Clinical_Response_LabeledBar.png", plot = plot_response_bar, width = 6, height = 5, dpi = 300)

plot_subtype_bar <- ggplot(subtype_summary, aes(x = AFR_Baseline_Group, y = Percent, fill = Subtype)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "black", fontface = "bold", size = 4) +
  scale_fill_manual(values = subtype_colors) +
  scale_y_continuous(labels = percent) +
  labs(title = "Breast Cancer Subtype Distribution", x = "Genomic Cohort", y = "Percentage of Patients", fill = "Clinical Subtype") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 11), axis.title = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

ggsave("Results/Figures/Clinical_Subtype_LabeledBar.png", plot = plot_subtype_bar, width = 6.5, height = 5, dpi = 300)

print("Labeled Bar Charts generated successfully!")
