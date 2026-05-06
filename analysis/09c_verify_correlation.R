library(dplyr)
library(vroom)
library(readr)
library(ggplot2)

# ==========================================
# USER CONFIGURATION
# ==========================================
target_gene <- "BLK"
target_pathway <- "GOBP_B_CELL_ACTIVATION"
output_dir <- "results/plots/correlations/verify"
# ==========================================

message("Loading data...")
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
ssgsea_scores <- readRDS("results/ssgsea/ssgsea_ontology_scores.rds")
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)

library(org.Hs.eg.db)
gene_id <- mapIds(org.Hs.eg.db, keys = target_gene, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")

if (is.na(gene_id)) {
  stop("Gene symbol '", target_gene, "' could not be mapped to an Ensembl ID.")
}

if (!gene_id %in% rownames(deg_mat)) {
  # Some matrices might have versions (ENSG... .1), check for matches
  matches <- grep(gene_id, rownames(deg_mat), value = TRUE)
  if (length(matches) > 0) {
    gene_id <- matches[1]
  } else {
    stop("Gene ID '", gene_id, "' (", target_gene, ") not found in degree matrix.")
  }
}

if (!target_pathway %in% rownames(ssgsea_scores)) {
  # List close matches if not found
  matches <- grep(target_pathway, rownames(ssgsea_scores), value = TRUE, ignore.case = TRUE)
  message("Pathway '", target_pathway, "' not found.")
  if (length(matches) > 0) {
    message("Did you mean one of these?\n  ", paste(matches, collapse = "\n  "))
  }
  stop("Invalid pathway name.")
}

common_samples <- intersect(colnames(deg_mat), colnames(ssgsea_scores)) %>%
  intersect(meta_sjs_ctrl$ID)

deg_vec <- as.numeric(deg_mat[gene_id, common_samples])
path_vec <- as.numeric(ssgsea_scores[target_pathway, common_samples])

res <- cor.test(deg_vec, path_vec, method = "spearman", exact = FALSE)
rho <- res$estimate
p_val <- res$p.value

message("\nResults for ", target_gene, " vs ", target_pathway, ":")
message("  Spearman Rho: ", round(rho, 4))
message("  P-value: ", format.pval(p_val))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cor_data <- data.frame(
  Degree = deg_vec,
  PathwayScore = path_vec,
  ID = common_samples
) %>%
left_join(meta_sjs_ctrl, by = "ID")

p <- ggplot(cor_data, aes(x = Degree, y = PathwayScore, color = Condition)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, aes(group = 1), color = "black", linetype = "dashed") +
  theme_minimal() +
  labs(title = paste(target_gene, "vs", target_pathway),
       subtitle = paste("Spearman Rho:", round(rho, 3), "| p-val:", format.pval(p_val)),
       x = paste("Gene Degree (", target_gene, ")"),
       y = "Pathway Activity (ssGSEA)",
       color = "Condition")

print(p)

clean_pathway <- gsub("GOBP_", "", target_pathway)
plot_filename <- file.path(output_dir, paste0(target_gene, "_", clean_pathway, ".svg"))
ggsave(plot_filename, p, width = 10, height = 7, bg = "white")
