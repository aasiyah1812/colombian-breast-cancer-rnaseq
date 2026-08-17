# -----------------------------------------------------------------------------
# ADMIXTURE CROSS-VALIDATION (CV) ERROR PLOT
# -----------------------------------------------------------------------------
library(tidyverse)

# Set working directory to the primary project folder
setwd("/path/to/project_directory")

# 1. Construct the dataset from the cross-validation log outputs
cv_data <- tibble(
  K = 2:10,
  CV_Error = c(
    0.60158, # K=2
    0.58034, # K=3
    0.56286, # K=4
    0.55738, # K=5
    0.55215, # K=6
    0.55166, # K=7
    0.55134, # K=8
    0.55085, # K=9
    0.55068  # K=10
  )
)

# 2. Generate the plot matching institutional figure style guidelines
cv_plot <- ggplot(cv_data, aes(x = K, y = CV_Error)) +
  # Draw the line (Steel Blue) and points (Burnt Orange) for visual distinction
  geom_line(color = "#286090", linewidth = 1.2) +
  geom_point(color = "#C55A28", size = 3.5) +
  
  # Ensure the x-axis explicitly shows every K value analyzed
  scale_x_continuous(breaks = 2:10) +
  
  # Apply axis and title labels
  labs(
    title = "ADMIXTURE Cross-Validation Error",
    x = "Number of Ancestral Populations (K)",
    y = "CV Error Score"
  ) +
  
  # Use theme_bw for a clean, academic white background with border
  theme_bw() +
  theme(
    # Format gridlines to highlight major intervals
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92"),
    
    # Adjust typography to match thesis figure standards
    plot.title = element_text(size = 18, hjust = 0, margin = margin(b = 10)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "gray30")
  )

# 3. Export the plot to the results directory
output_file <- "Results/Figures/Supplementary_Figure_S2_CV_Error.png"
ggsave(output_file, plot = cv_plot, width = 8, height = 5.5, dpi = 300)

print(paste("Cross-validation plot successfully saved to:", output_file))
