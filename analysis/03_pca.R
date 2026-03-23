library(dplyr)
library(tibble)
library(vroom)
library(ggplot2)
library(uwot)
library(RColorBrewer)

# 1. Load Data
# Expression Matrix (VST Normalized)
expr_mat <- readRDS("results/vst_normalized_counts_symbols.rds")

# Metadata
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv")
meta_detailed <- vroom("data/precisesads/metadata_sjs_detailed.csv")
clusters <- vroom("data/precisesads/sjogren_clusters.csv")

# Merge detailed clinical info and clusters onto the basic SJS vs Ctrl metadata
metadata <- meta_sjs_ctrl %>%
  left_join(meta_detailed %>% dplyr::select(-Sex, -Cohort), by = "ID") %>%
  left_join(clusters, by = "ID") %>%
  mutate(
    # AcquisitionGroup as defined in 02_digs.R (Late vs Early Acquisition)
    # This is only for SJS patients with age and duration data
    AcquisitionAge = Age - DiseaseDuration_Years,
    Menopause = ifelse(Condition == "Sjogrens" & Sex == "Female",
                       ifelse(AcquisitionAge > 48, "Post-menopause (Late)", "Pre-menopause (Early)"),
                       NA)
  )

# 3. Align Expression Matrix and Metadata
common_samples <- intersect(colnames(expr_mat), metadata$ID)
expr_mat <- expr_mat[, common_samples]
metadata <- metadata %>% 
  filter(ID %in% common_samples) %>% 
  arrange(match(ID, common_samples))

# Double check alignment
if(!all(colnames(expr_mat) == metadata$ID)) stop("Samples are not aligned!")

# 4. Feature Selection: Top 1000 most variable genes
gene_vars <- apply(expr_mat, 1, var)
keep_genes <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
expr_subset <- expr_mat[keep_genes, ]

# 5. Perform PCA
X <- t(expr_subset)
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

# Variance explained
var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Prepare dataframe for plotting
pca_df <- as.data.frame(pca_res$x[, 1:10]) %>%
  rownames_to_column("sample") %>%
  left_join(metadata, by = c("sample" = "ID"))

# 6. Generate PCA Plots
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# Function to save PCA plot
save_pca_plot <- function(df, color_var, title_suffix, filename) {
  # Remove NA for the specific coloring variable to have cleaner plots if needed
  # or keep them to see where they fall. For Menopause/Cluster we keep them as "Unknown/Control"
  
  plot_df <- df
  if (color_var == "PREDICTION") {
    plot_df$PREDICTION <- as.factor(plot_df$PREDICTION)
  }
  
  p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = !!sym(color_var))) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() + 
    labs(color = title_suffix,
         title = paste0("PCA of Expression Matrix (Top 1000 genes) - by ", title_suffix),
         x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
    scale_color_brewer(palette = "Set1", na.value = "grey80")
  
  ggsave(filename, p, width = 8, height = 6, dpi = 300)
}

# PCA by Condition (Group)
save_pca_plot(pca_df, "Condition", "Condition", "results/plots/PCA_expression_by_group.png")

# PCA by Cluster
save_pca_plot(pca_df, "PREDICTION", "Cluster", "results/plots/PCA_expression_by_cluster.png")

# PCA by Menopause (AcquisitionGroup)
save_pca_plot(pca_df, "Menopause", "Menopause", "results/plots/PCA_expression_by_menopause.png")

set.seed(123)
npc_use <- 30 # Using more PCs for UMAP
umap_res <- uwot::umap(pca_res$x[, 1:npc_use], 
                       n_neighbors = 15, 
                       min_dist = 0.1, 
                       metric = "euclidean",
                       verbose = FALSE)

umap_df <- as.data.frame(umap_res) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(sample = rownames(X)) %>%
  left_join(metadata, by = c("sample" = "ID"))

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(color = "Condition",
       title = "UMAP of Expression Matrix (Top 30 PCs)") +
  scale_color_brewer(palette = "Set1")

ggsave("results/plots/UMAP_expression_by_group.png", umap_plot, width = 8, height = 6, dpi = 300)