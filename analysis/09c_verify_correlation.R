library(dplyr)
library(vroom)
library(readr)
library(ggplot2)

# ==========================================
# USER CONFIGURATION
# ==========================================
target_genes <- c("OASL", "LGALS3BP", "ISG15", "UBE2L6")
target_pathway <- "GOBP_DEFENSE_RESPONSE_TO_VIRUS"
output_dir <- "results/plots/correlations/verify"
# ==========================================

message("Loading data...")
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
ssgsea_scores <- readRDS("results/ssgsea/ssgsea_ontology_scores.rds")
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)

library(org.Hs.eg.db)

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

sjs_samples <- meta_sjs_ctrl$ID[meta_sjs_ctrl$Condition == "Sjogrens"] %>% intersect(common_samples)
ctrl_samples <- meta_sjs_ctrl$ID[meta_sjs_ctrl$Condition == "Control"] %>% intersect(common_samples)

path_vec <- as.numeric(ssgsea_scores[target_pathway, common_samples])
path_vec_sjs <- as.numeric(ssgsea_scores[target_pathway, sjs_samples])
path_vec_ctrl <- as.numeric(ssgsea_scores[target_pathway, ctrl_samples])

all_cor_data <- list()

for (target_gene in target_genes) {
  gene_id <- mapIds(org.Hs.eg.db, keys = target_gene, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")

  if (is.na(gene_id)) {
    warning("Gene symbol '", target_gene, "' could not be mapped to an Ensembl ID.")
    next
  }

  if (!gene_id %in% rownames(deg_mat)) {
    # Some matrices might have versions (ENSG... .1), check for matches
    matches <- grep(gene_id, rownames(deg_mat), value = TRUE)
    if (length(matches) > 0) {
      gene_id <- matches[1]
    } else {
      warning("Gene ID '", gene_id, "' (", target_gene, ") not found in degree matrix.")
      next
    }
  }

  deg_vec <- as.numeric(deg_mat[gene_id, common_samples])
  
  deg_vec_sjs <- as.numeric(deg_mat[gene_id, sjs_samples])
  res_sjs <- tryCatch(cor.test(deg_vec_sjs, path_vec_sjs, method = "spearman", exact = FALSE), error = function(e) list(estimate = NA, p.value = NA))
  
  deg_vec_ctrl <- as.numeric(deg_mat[gene_id, ctrl_samples])
  res_ctrl <- tryCatch(cor.test(deg_vec_ctrl, path_vec_ctrl, method = "spearman", exact = FALSE), error = function(e) list(estimate = NA, p.value = NA))

  res <- cor.test(deg_vec, path_vec, method = "spearman", exact = FALSE)
  rho <- res$estimate
  p_val <- res$p.value

  message("\nResults for ", target_gene, " vs ", target_pathway, ":")
  message("  Spearman Rho: ", round(rho, 4))
  message("  P-value: ", format.pval(p_val))
  message("  SD Rho: ", round(res_sjs$estimate, 4), " (p=", format.pval(res_sjs$p.value, digits=2), ")")
  message("  Ctrl Rho: ", round(res_ctrl$estimate, 4), " (p=", format.pval(res_ctrl$p.value, digits=2), ")")

  facet_label <- sprintf(
    "%s\nSD Rho: %.3f (p=%s)\nCtrl Rho: %.3f (p=%s)",
    target_gene,
    res_sjs$estimate, format.pval(res_sjs$p.value, digits = 2),
    res_ctrl$estimate, format.pval(res_ctrl$p.value, digits = 2)
  )

  cor_data <- data.frame(
    Degree = deg_vec,
    PathwayScore = path_vec,
    ID = common_samples,
    Gene = target_gene,
    FacetLabel = facet_label
  ) %>%
  left_join(meta_sjs_ctrl, by = "ID")

  all_cor_data[[target_gene]] <- cor_data
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

combined_data <- do.call(rbind, all_cor_data)

# Ensure FacetLabel is a factor so it plots in the order specified
facet_levels <- unique(combined_data$FacetLabel[order(match(combined_data$Gene, target_genes))])
combined_data$FacetLabel <- factor(combined_data$FacetLabel, levels = facet_levels)

p <- ggplot(combined_data, aes(x = Degree, y = PathwayScore, color = Condition, fill = Condition)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, linetype = "dashed") +
  facet_wrap(~ FacetLabel, scales = "free_x") +
  scale_color_manual(values = c("Control" = "#3498db", "Sjogrens" = "#e74c3c")) +
  scale_fill_manual(values = c("Control" = "#3498db", "Sjogrens" = "#e74c3c")) +
  theme_minimal() +
  labs(
       x = "Gene Degree",
       y = "Pathway Activity (ssGSEA)",
       color = "Condition", fill = "Condition") +
  theme(strip.text = element_text(size = 10, face = "bold"))

print(p)

clean_pathway <- gsub("GOBP_", "", target_pathway)
plot_filename <- file.path(output_dir, paste0("Grid_Genes_", clean_pathway, ".svg"))
ggsave(plot_filename, p, width = 10, height = 7, bg = "white")
plot_filename <- file.path(output_dir, paste0("Grid_Genes_", clean_pathway, ".pdf"))
ggsave(plot_filename, p, width = 10, height = 7, bg = "white")