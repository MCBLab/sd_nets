library(sva)
library(vroom)
library(ggplot2)
library(dplyr)
library(tibble)
library(Rtsne)

message("--- Loading Data ---")

file1 <- 'data/ERP129369_STAR/salmon.merged.gene_counts.tsv'
file2 <- 'data/SRP287550_STAR/salmon.merged.gene_counts.tsv'
file3 <- 'data/SRP279407/salmon.merged.gene_counts.tsv'

data1 <- vroom(file1) |> column_to_rownames("gene_id") |> select(-gene_name)
data2 <- vroom(file2) |> column_to_rownames("gene_id") |> select(-gene_name)
data3 <- vroom(file3) |> column_to_rownames("gene_id") |> select(-gene_name)

metadata1 <- vroom("data/ERP129369_kallisto/metadata_ERP129369.csv", delim = ",") |>
  select(sample = Run, condition = `Experimental_Factor:_disease (exp)`, sex) |>
  mutate(condition = ifelse(condition == "Sjogren's syndrome", "SD", "Sicca")) |> 
  # filter(condition == "SD") |> 
  filter(sex == "female") |> 
  select(-sex)

metadata2 <- vroom("data/SRP287550_kallisto/metadata_SRP287550.csv", delim = ",") |>
  select(sample = Run, condition = diagnosis) |>
  mutate(condition = ifelse(condition == "pSS", "SD", "non-SD"))

metadata3 <- vroom("data/SRP279407/samplesheet_rnaseq.csv", delim = ",") |>
  select(sample, condition) |>
  mutate(condition = ifelse(condition == "sjogren", "SD", "non-SD"))

data1 <- data1[colnames(data1) %in% metadata1$sample]
data3 <- data3[colnames(data3) %in% metadata3$sample]

metadata <- rbind(metadata1, metadata2, metadata3)

message("\n--- Generating Individual PCA and t-SNE Plots ---")

plot_individual_analysis <- function(expression_data, meta, dataset_name) {
  message(paste("Running PCA for", dataset_name, "..."))
  
  variances <- apply(expression_data, 1, var)
  ordered_variances <- order(variances, decreasing = TRUE)
  top_genes <- head(ordered_variances, 2000)
  expression_subset <- expression_data[top_genes, ]
  
  pca_result <- prcomp(t(expression_subset), scale. = TRUE)
  pca_summary <- summary(pca_result)
  pc1_var <- round(pca_summary$importance[2, 1] * 100, 1)
  pc2_var <- round(pca_summary$importance[2, 2] * 100, 1)
  
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample = rownames(pca_result$x)
  ) |>
    left_join(meta, by = "sample")

  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = paste("PCA of Top 2000 Variable Genes in", dataset_name),
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p_pca)

  message(paste("Running t-SNE for", dataset_name, "..."))
  num_pcs <- min(20, ncol(pca_result$x))
  perplexity_val <- max(1, floor((nrow(pca_result$x) - 1) / 3) -1)

  set.seed(42) # for reproducibility
  tsne_result <- Rtsne(pca_result$x[, 1:num_pcs], pca = FALSE, perplexity = perplexity_val)
  
  tsne_data <- data.frame(
    TSNE1 = tsne_result$Y[, 1],
    TSNE2 = tsne_result$Y[, 2],
    sample = rownames(pca_result$x)
  ) |>
    left_join(meta, by = "sample")

  p_tsne <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = condition)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = paste("t-SNE of Top 20 PCs in", dataset_name)
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  print(p_tsne)
}

plot_individual_analysis(data1, metadata1, "ERP129369")
plot_individual_analysis(data2, metadata2, "SRP287550")
plot_individual_analysis(data3, metadata3, "SRP279407")

message("\n--- Generating Combined PCA Plot (Before Correction) ---")

common_genes <- Reduce(intersect, list(rownames(data1), rownames(data2), rownames(data3)))
message(paste("Found", length(common_genes), "common genes across all datasets."))

expression1_common <- data1[common_genes, ]
expression2_common <- data2[common_genes, ]
expression3_common <- data3[common_genes, ]
merged_expression <- cbind(expression1_common, expression2_common, expression3_common)

variances_merged <- apply(merged_expression, 1, var)
ordered_variances_merged <- order(variances_merged, decreasing = TRUE)
top_2000_genes_merged <- head(ordered_variances_merged, 2000)
# merged_expression_top2000 <- merged_expression[variances_merged > 0, ]
merged_expression_top2000 <- merged_expression[top_2000_genes_merged, ]

message("Filtered combined data to the top 2000 most variable genes.")

batch_info <- data.frame(
  sample = colnames(merged_expression_top2000),
  batch = c(
    rep('ERP129369', ncol(expression1_common)),
    rep('SRP287550', ncol(expression2_common)),
    rep('SRP279407', ncol(expression3_common))
  ),
  row.names = colnames(merged_expression_top2000)
) |>
  left_join(metadata, by = "sample")

pca_result_before <- prcomp(t(merged_expression_top2000), scale. = TRUE)
pca_summary_before <- summary(pca_result_before)
pc1_var_before <- round(pca_summary_before$importance[2, 1] * 100, 1)
pc2_var_before <- round(pca_summary_before$importance[2, 2] * 100, 1)

pca_data_before <- data.frame(
  PC1 = pca_result_before$x[, 1],
  PC2 = pca_result_before$x[, 2],
  Batch = batch_info$batch,
  condition = batch_info$condition
)

plot_before <- ggplot(pca_data_before, aes(x = PC1, y = PC2, color = Batch, shape = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA Before Batch Correction",
    x = paste0("PC1 (", pc1_var_before, "%)"),
    y = paste0("PC2 (", pc2_var_before, "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))
print(plot_before)

message("\n--- Generating Combined t-SNE Plot (Before Correction) ---")
set.seed(42)
safe_perplexity_before <- max(1, floor((nrow(pca_result_before$x) - 1) / 3) - 1)
tsne_result_before <- Rtsne(pca_result_before$x[, 1:20], pca = FALSE, perplexity = safe_perplexity_before)
tsne_data_before <- data.frame(
  TSNE1 = tsne_result_before$Y[, 1],
  TSNE2 = tsne_result_before$Y[, 2],
  Batch = batch_info$batch,
  condition = batch_info$condition
)
plot_tsne_before <- ggplot(tsne_data_before, aes(x = TSNE1, y = TSNE2, color = Batch, shape = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE Before Batch Correction") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))
print(plot_tsne_before)


message("\n--- Running ComBat for Batch Correction ---")
batch_vector <- batch_info$batch
modcombat <- model.matrix(~1, data = batch_info) 

corrected_expression <- ComBat_seq(
  as.matrix(merged_expression),
  batch_vector,
  covar_mod = modcombat,
  group = batch_info$condition
)
message("ComBat finished.")

message("\n--- Generating Combined PCA Plot (After Correction) ---")

variances_merged <- apply(corrected_expression, 1, var)
ordered_variances_merged <- order(variances_merged, decreasing = TRUE)
top_2000_genes_merged <- head(ordered_variances_merged, 2000)
# merged_expression_top2000 <- merged_expression[variances_merged > 0, ]
corrected_expression_top2000 <- corrected_expression[top_2000_genes_merged, ]

pca_result_after <- prcomp(t(corrected_expression_top2000), scale. = TRUE)
pca_summary_after <- summary(pca_result_after)
pc1_var_after <- round(pca_summary_after$importance[2, 1] * 100, 1)
pc2_var_after <- round(pca_summary_after$importance[2, 2] * 100, 1)

pca_data_after <- data.frame(
  PC1 = pca_result_after$x[, 1],
  PC2 = pca_result_after$x[, 2],
  Batch = batch_info$batch,
  condition = batch_info$condition
)

plot_after <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch, shape = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA After Batch Correction (ComBat)",
    x = paste0("PC1 (", pc1_var_after, "%)"),
    y = paste0("PC2 (", pc2_var_after, "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))
print(plot_after)

message("\n--- Generating Combined t-SNE Plot (After Correction) ---")
set.seed(42)
safe_perplexity_after <- max(1, floor((nrow(pca_result_after$x) - 1) / 3) - 1)
tsne_result_after <- Rtsne(pca_result_after$x[, 1:20], pca = FALSE, perplexity = safe_perplexity_after)
tsne_data_after <- data.frame(
  TSNE1 = tsne_result_after$Y[, 1],
  TSNE2 = tsne_result_after$Y[, 2],
  Batch = batch_info$batch,
  condition = batch_info$condition
)
plot_tsne_after <- ggplot(tsne_data_after, aes(x = TSNE1, y = TSNE2, color = Batch, shape = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE After Batch Correction (ComBat)") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))
print(plot_tsne_after)

