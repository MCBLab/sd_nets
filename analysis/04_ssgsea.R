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

ssgsea_res <- gsva(params)
dir.create("results/ssgsea", showWarnings = FALSE)
saveRDS(ssgsea_res, "results/ssgsea/ssgsea_hallmark_scores.rds")
write.csv(as.data.frame(ssgsea_res), "results/ssgsea/ssgsea_scores_matrix.csv")

metadata_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv")
clusters <- vroom("data/precisesads/sjogren_clusters.csv")
metadata_clusters <- metadata_sjs_ctrl %>% 
  inner_join(clusters, by = "ID") %>%
  filter(ID %in% colnames(ssgsea_res))

ssgsea_res <- ssgsea_res[, metadata_clusters$ID]

group <- factor(paste0("Cluster_", metadata_clusters$PREDICTION))
design <- model.matrix(~ 0 + group)

# Limpar os nomes das colunas da matriz de desenho para remover o prefixo "group" que o R coloca
colnames(design) <- levels(group)

fit <- lmFit(ssgsea_res, design)

cluster_levels <- levels(group)

for (cl in cluster_levels) {
  
  others <- cluster_levels[cluster_levels != cl]
  
  contrast_formula <- paste0(cl, " - (", paste(others, collapse = " + "), ") / ", length(others))
  
  # Agora o makeContrasts funcionará corretamente
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  fit_contrast <- contrasts.fit(fit, contrast_matrix)
  fit_contrast <- eBayes(fit_contrast)
  
  res <- topTable(fit_contrast, coef = 1, number = Inf) %>%
    rownames_to_column("pathway") %>%
    rename(
      diff_score = logFC,
      p_val = P.Value,
      padj = adj.P.Val,
      t_stat = t
    ) %>%
    arrange(p_val)
  
  clean_name <- gsub("Cluster_", "cluster_", cl)
  write.csv(res, paste0("results/ssgsea/limma_pathway_", clean_name, "_vs_others.csv"), row.names = FALSE)
}

