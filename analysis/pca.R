library(dplyr)
library(tibble)
library(vroom)
library(ggplot2)
library(uwot)
library(RColorBrewer)

deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
metadata <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

common_samples <- intersect(colnames(deg_mat), metadata$ID)
deg_mat <- deg_mat[, common_samples]
metadata <- metadata %>% 
  filter(ID %in% common_samples) %>% 
  arrange(match(ID, common_samples))

total_degrees <- colSums(deg_mat)
scaling_factors <- total_degrees / median(total_degrees)
deg_norm <- sweep(deg_mat, 2, scaling_factors, "/")
deg_norm <- log1p(deg_norm)

# Filter by Variance (Matching Notebook Logic: sd > 1)
gene_sds <- apply(deg_norm, 1, sd)
keep_genes <- names(gene_sds[gene_sds > 1])
deg_subset <- deg_norm[keep_genes, ]

cat("Genes remaining after SD > 1 filter:", nrow(deg_subset), "\n")

X <- t(deg_subset)

pca_res <- prcomp(X, center = TRUE, scale. = TRUE)
# Use top 40 PCs as defined in the original notebook function
npc_use <- min(40, ncol(pca_res$x))
var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

pca_df <- as.data.frame(pca_res$x[, 1:npc_use]) %>%
  rownames_to_column("sample") %>%
  mutate(group = metadata$Condition[match(sample, metadata$ID)])

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() + 
  labs(color = "Condition",
       title = "PCA of Degree Matrix (All genes with SD > 1)",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)"))

ggsave("results/plots/PCA_degree_SD_filtered.png", pca_plot, width = 8, height = 6, dpi = 300)

set.seed(123)
umap_res <- uwot::umap(pca_res$x[, 1:npc_use], 
                       n_neighbors = 15, 
                       min_dist = 0.1, 
                       metric = "euclidean",
                       verbose = FALSE)

umap_df <- as.data.frame(umap_res) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(sample = rownames(X),
         group = metadata$Condition[match(sample, metadata$ID)])

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(color = "Condition",
       title = paste0("UMAP of Degree Matrix (", npc_use, " PCs)"))

ggsave("results/plots/UMAP_degree_SD_filtered.png", umap_plot, width = 8, height = 6, dpi = 300)