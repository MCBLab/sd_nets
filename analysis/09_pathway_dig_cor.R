library(dplyr)
library(vroom)
library(readr)
library(tidyr)
library(clusterProfiler)
library(purrr)
library(ggplot2)

# --- 1. Load Data ---
message("Loading data...")
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
ssgsea_scores <- readRDS("results/ssgsea/ssgsea_hallmark_scores.rds")
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)
digs_sjs <- read_csv("results/digs/sjs_vs_ctrl_significant_padj05.csv", show_col_types = FALSE)
hallmark_gmt <- read.gmt("data/h.all.v2026.1.Hs.symbols.gmt")

# --- 2. Prepare Combinations ---
# Filter hallmark pathways for genes that are significant DIGs in sjs_vs_ctrl
combinations <- hallmark_gmt %>% 
  filter(gene %in% digs_sjs$symbol) %>%
  rename(symbol = gene, pathway = term) %>%
  left_join(digs_sjs %>% select(gene, symbol), by = "symbol")

message("Total combinations to analyze: ", nrow(combinations))

# --- 3. Perform Correlation Analysis ---
common_samples <- intersect(colnames(deg_mat), colnames(ssgsea_scores)) %>%
  intersect(meta_sjs_ctrl$ID)

message("Number of samples for correlation: ", length(common_samples))

# Align matrices
deg_mat_subset <- deg_mat[, common_samples]
ssgsea_scores_subset <- ssgsea_scores[, common_samples]

results <- combinations %>%
  mutate(correlation_results = map2(gene, pathway, function(g_id, p_name) {
    # Check if ID exists in matrices
    if (!g_id %in% rownames(deg_mat_subset) || !p_name %in% rownames(ssgsea_scores_subset)) {
      return(NULL)
    }
    
    deg_vec <- as.numeric(deg_mat_subset[g_id, ])
    path_vec <- as.numeric(ssgsea_scores_subset[p_name, ])
    
    # Spearman Correlation
    res <- cor.test(deg_vec, path_vec, method = "spearman", exact = FALSE)
    
    data.frame(
      rho = res$estimate,
      p_val = res$p.value,
      stringsAsFactors = FALSE
    )
  })) %>%
  unnest(correlation_results) %>%
  mutate(padj = p.adjust(p_val, method = "BH")) %>%
  arrange(p_val)

# --- 4. Save Results ---
dir.create("results/repurposing", showWarnings = FALSE, recursive = TRUE)
output_path <- "results/repurposing/sjs_dig_pathway_correlations.csv"
write_csv(results, output_path)

message("Correlation analysis complete. Results saved to: ", output_path)
message("Significant correlations (padj < 0.05): ", sum(results$padj < 0.05))

# --- 5. Generate Plots for Significant Correlations ---
filtered <- results %>%
  filter(padj < 0.05)

if (nrow(filtered) > 0) {
  message("Generating ", nrow(filtered), " plots...")
  dir.create("results/plots/correlations", showWarnings = FALSE, recursive = TRUE)
  
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
    
    # Clean pathway name for filename
    clean_pathway <- gsub("HALLMARK_", "", row$pathway)
    plot_filename <- paste0("results/plots/correlations/cor_", row$symbol, "_", clean_pathway, ".png")
    ggsave(plot_filename, p, width = 8, height = 7, bg = "white")
  }
} else {
  message("No significant correlations to plot.")
}
