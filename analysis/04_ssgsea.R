library(GSVA)
library(GSEABase)
library(dplyr)
library(tibble)
library(vroom)
library(limma)

normalized_counts <- readRDS("results/vst_normalized_counts_symbols.rds")
gene_sets <- getGmt("data/h.all.v2026.1.Hs.symbols.gmt")

params <- ssgseaParam(
  exprData = as.matrix(normalized_counts),
  geneSets = gene_sets,
  minSize = 10,
  maxSize = 500
)

ssgsea_res_all <- gsva(params)
dir.create("results/ssgsea", showWarnings = FALSE)
saveRDS(ssgsea_res_all, "results/ssgsea/ssgsea_hallmark_scores.rds")
write.csv(as.data.frame(ssgsea_res_all), "results/ssgsea/ssgsea_scores_matrix.csv")

metadata_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)
clusters <- vroom("data/precisesads/sjogren_clusters.csv", show_col_types = FALSE)

# 1. Cluster comparisons (subset of samples)
metadata_clusters <- metadata_sjs_ctrl %>% 
  inner_join(clusters, by = "ID") %>%
  filter(ID %in% colnames(ssgsea_res_all))

if (nrow(metadata_clusters) > 0) {
  ssgsea_res_clusters <- ssgsea_res_all[, metadata_clusters$ID]
  group <- factor(paste0("Cluster_", metadata_clusters$PREDICTION))
  
  if (length(levels(group)) >= 2) {
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    fit <- lmFit(ssgsea_res_clusters, design)
    
    for (cl in levels(group)) {
      others <- levels(group)[levels(group) != cl]
      contrast_formula <- paste0(cl, " - (", paste(others, collapse = " + "), ") / ", length(others))
      contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
      
      fit_contrast <- contrasts.fit(fit, contrast_matrix)
      fit_contrast <- eBayes(fit_contrast)
      
      res <- topTable(fit_contrast, coef = 1, number = Inf) %>%
        rownames_to_column("pathway") %>%
        rename(diff_score = logFC, p_val = P.Value, padj = adj.P.Val, t_stat = t) %>%
        arrange(p_val)
      
      clean_name <- gsub("Cluster_", "cluster_", cl)
      write.csv(res, paste0("results/ssgsea/limma_pathway_", clean_name, "_vs_others.csv"), row.names = FALSE)
    }
  }
}

# 2. Sjogrens vs Control comparison (all available samples)
metadata_general <- metadata_sjs_ctrl %>%
  filter(ID %in% colnames(ssgsea_res_all))

if (nrow(metadata_general) > 0) {
  ssgsea_res_general <- ssgsea_res_all[, metadata_general$ID]
  condition <- factor(metadata_general$Condition)
  
  if (length(levels(condition)) >= 2) {
    design_sjs <- model.matrix(~ 0 + condition)
    colnames(design_sjs) <- levels(condition)
    
    fit_sjs <- lmFit(ssgsea_res_general, design_sjs)
    # Ensure Control and Sjogrens levels exist
    if ("Sjogrens" %in% colnames(design_sjs) && "Control" %in% colnames(design_sjs)) {
        contrast_matrix_sjs <- makeContrasts(Sjogrens - Control, levels = design_sjs)
        fit_contrast_sjs <- contrasts.fit(fit_sjs, contrast_matrix_sjs)
        fit_contrast_sjs <- eBayes(fit_contrast_sjs)
        
        res_sjs <- topTable(fit_contrast_sjs, coef = 1, number = Inf) %>%
          rownames_to_column("pathway") %>%
          rename(diff_score = logFC, p_val = P.Value, padj = adj.P.Val, t_stat = t) %>%
          arrange(p_val)
        
        write.csv(res_sjs, "results/ssgsea/limma_pathway_sjs_vs_ctrl.csv", row.names = FALSE)
    }
  }
}
