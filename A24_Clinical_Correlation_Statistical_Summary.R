# -----------------------------------------------------------------------------
# SAVE CLINICAL STATISTICS TABLE
# -----------------------------------------------------------------------------
library(tidyverse)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# Structure the findings into a standardized matrix for export
clinical_stats_table <- tibble(
  `Clinical Variable` = c("Neoadjuvant Therapy Response", "Clinical Breast Cancer Subtype"),
  `Statistical Test Used` = c("Fisher's Exact Test", "Fisher's Exact Test (Monte Carlo simulated)"),
  `p-value` = c(signif(fisher_response$p.value, 3), signif(fisher_subtype$p.value, 3)),
  `Statistical Significance (alpha = 0.05)` = c("Not Significant", "Not Significant (Trend)")
)

# Export the matrix structure utilizing the localized format command
write_excel_csv2(clinical_stats_table, "Results/Tables/Clinical_Correlations_Statistics.csv")

print("Clinical statistics table successfully saved.")
