library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(tibble)
library(vroom)

# --- 1. Load Data ---
ssgsea_res <- readRDS("results/ssgsea/ssgsea_hallmark_scores.rds")
clusters_meta <- vroom("data/precisesads/sjogren_clusters.csv") %>%
  rename(Cluster = PREDICTION)

# Metadata for Sjögren patients (if needed to filter or order samples)
meta_sjs <- vroom("data/precisesads/metadata_sjs_ctrl.csv") %>%
  filter(Condition == "Sjögren")

# Sync samples: only Sjögren samples in both matrix and cluster metadata
common_samples <- intersect(colnames(ssgsea_res), clusters_meta$ID)
ssgsea_res <- ssgsea_res[, common_samples]
clusters_meta <- clusters_meta %>%
  filter(ID %in% common_samples) %>%
  arrange(match(ID, common_samples))

# Ensure output directory exists
dir.create("results/plots/ssgsea", showWarnings = FALSE, recursive = TRUE)

# --- 2. Heatmap: Top Cluster-Specific Pathways ---

# Load Limma results for each cluster and find the most significant pathways
cluster_ids <- 1:4
all_pathways_list <- list()

for (cl in cluster_ids) {
  res_file <- paste0("results/ssgsea/limma_pathway_cluster_", cl, "_vs_others.csv")
  if (file.exists(res_file)) {
    # We'll take top 10 pathways per cluster (by significance)
    cl_res <- read_csv(res_file) %>%
      filter(padj < 0.05) %>%
      arrange(p_val) %>%
      head(10) %>%
      pull(pathway)
    all_pathways_list[[as.character(cl)]] <- cl_res
  }
}

# Unique union of top pathways
signature_pathways <- unique(unlist(all_pathways_list))

# Filter ssgsea_res for these pathways
heatmap_mat <- ssgsea_res[signature_pathways, ]

# Row annotations: which cluster(s) a pathway is "top" for
pathway_anno <- data.frame(pathway = signature_pathways) %>%
  mutate(
    Cluster_Association = sapply(pathway, function(p) {
      cl_found <- names(all_pathways_list)[sapply(all_pathways_list, function(l) p %in% l)]
      paste("Cluster", cl_found, collapse = " & ")
    })
  ) %>%
  column_to_rownames("pathway")

# Column annotations: sample cluster assignments
col_anno <- clusters_meta %>%
  mutate(Cluster = paste0("Cluster ", Cluster)) %>%
  column_to_rownames("ID")

clusters_meta <- clusters_meta %>%
  mutate(Cluster_Label = paste0("Cluster ", Cluster)) %>%
  arrange(Cluster) # Ordena amostras: 1, 2, 3, 4

heatmap_mat <- ssgsea_res[signature_pathways, clusters_meta$ID]

# 2. Z-score CORRETO (Garante que a matriz mantenha a forma original)
# t(apply) corrige a transposição automática do scale
heatmap_mat_scaled <- t(apply(heatmap_mat, 1, scale))
colnames(heatmap_mat_scaled) <- colnames(heatmap_mat)
rownames(heatmap_mat_scaled) <- rownames(heatmap_mat)

# 3. Preparar Anotações
col_anno <- clusters_meta %>%
  select(ID, Cluster_Label) %>%
  column_to_rownames("ID")

colnames(col_anno) <- "Cluster"

ann_colors <- list(
  Cluster = c("Cluster 1" = "#E41A1C", 
              "Cluster 2" = "#377EB8", 
              "Cluster 3" = "#4DAF4A", 
              "Cluster 4" = "#984EA3")
)

# 4. Plot Heatmap com Gaps entre Clusters
# Calculamos onde os clusters mudam para colocar um espaço visual
gaps <- cumsum(table(clusters_meta$Cluster_Label)[unique(clusters_meta$Cluster_Label)])

png("results/plots/ssgsea/heatmap_pathway_signature.png", width = 1200, height = 1000, res = 150)
pheatmap(
  heatmap_mat_scaled,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  cluster_cols = FALSE, # Mantemos FALSE pois já ordenamos manualmente
  cluster_rows = TRUE,  # Agrupa vias similares (ex: IFN-a e IFN-g juntas)
  gaps_col = gaps,      # Cria a separação visual entre os grupos
  main = "Atividade de Vias por Cluster (Z-score)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(-2, 2, length.out = 101), # Limita a escala para aumentar o contraste
  fontsize_row = 8
)
dev.off()

# --- 3. PCA of Pathway Activities ---

# Transpose matrix: samples in rows, pathways in columns
pca_input <- t(ssgsea_res)

# PCA calculation
pca_res <- prcomp(pca_input, scale. = TRUE)

# Prepare dataframe for plotting
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("ID") %>%
  left_join(clusters_meta, by = "ID") %>%
  mutate(Cluster = factor(paste0("Cluster ", Cluster)))

# Variance explained
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Plot PCA
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = ann_colors$Cluster) +
  theme_minimal() +
  labs(
    title = "PCA of Pathway Activities (ssGSEA Hallmark Scores)",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme(legend.position = "right")

ggsave("results/plots/ssgsea/pca_pathways_by_cluster.png", p_pca, width = 8, height = 6, bg = "white")