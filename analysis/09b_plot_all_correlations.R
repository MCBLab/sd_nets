library(dplyr)
library(vroom)
library(readr)
library(ggplot2)

# Correlation results
results <- read_csv("results/repurposing/sjs_dig_pathway_correlations.csv", show_col_types = FALSE)

# Filter for significant results (padj < 0.05)
filtered <- results %>%
  filter(padj < 0.05)

message("Generating ", nrow(filtered), " plots...")

# Network and Pathway scores
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
ssgsea_scores <- readRDS("results/ssgsea/ssgsea_hallmark_scores.rds")

# Metadata
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)

dir.create("results/plots/correlations", showWarnings = FALSE, recursive = TRUE)

common_samples <- intersect(colnames(deg_mat), colnames(ssgsea_scores)) %>%
  intersect(meta_sjs_ctrl$ID)

# Align matrices
deg_mat_subset <- deg_mat[, common_samples]
ssgsea_scores_subset <- ssgsea_scores[, common_samples]

for (i in 1:nrow(filtered)) {
  row <- filtered[i, ]
  message("  - Plotting ", i, "/", nrow(filtered), ": ", row$symbol, " vs ", row$pathway)
  
  cor_data <- data.frame(
    Degree = as.numeric(deg_mat_subset[row$gene, ]),
    PathwayScore = as.numeric(ssgsea_scores_subset[row$pathway, ]),
    ID = common_samples
  ) %>%
  left_join(meta_sjs_ctrl, by = "ID")
  
  p <- ggplot(cor_data, aes(x = Degree, y = PathwayScore, color = Condition)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, aes(group = 1), color = "black", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste(row$symbol, "vs", row$pathway),
         subtitle = paste("Spearman Rho:", round(row$rho, 3), "| padj:", format.pval(row$padj)),
         x = paste("Gene Degree (", row$symbol, ")"),
         y = "Pathway Activity (ssGSEA)",
         color = "Condition")
  
  plot_filename <- paste0("results/plots/correlations/cor_", row$symbol, "_", row$pathway, ".png")
  ggsave(plot_filename, p, width = 8, height = 7, bg = "white")
}
