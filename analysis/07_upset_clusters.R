library(UpSetR)
library(dplyr)

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

cluster_files <- list.files("results/digs", pattern = "cluster_[0-9]_vs_others_significant_padj05.csv$", full.names = TRUE)
message("Found ", length(cluster_files), " cluster files.")

gene_list <- list()
for (f in cluster_files) {
  cluster_name <- gsub("_significant_padj05.csv", "", basename(f))
  # Map "cluster_1_vs_others" to "Cluster 1" for better plot labels
  clean_name <- gsub("_vs_others", "", cluster_name)
  clean_name <- gsub("cluster_", "Cluster ", clean_name)
  
  df <- read.csv(f)
  if (nrow(df) > 0) {
    gene_list[[clean_name]] <- as.character(df$gene)
    message("Read ", length(gene_list[[clean_name]]), " genes for ", clean_name)
  } else {
    message("No significant genes for ", clean_name)
  }
}

if (length(gene_list) > 1) {
  all_unique_genes <- unique(unlist(gene_list))
  input_mat <- matrix(0, nrow = length(all_unique_genes), ncol = length(gene_list))
  colnames(input_mat) <- names(gene_list)
  rownames(input_mat) <- all_unique_genes
  
  for (name in names(gene_list)) {
    input_mat[gene_list[[name]], name] <- 1
  }
  input_df <- as.data.frame(input_mat)
  message("Final input matrix: ", nrow(input_df), " rows, ", ncol(input_df), " columns.")
  
  png("results/plots/upset_clusters.png", width = 1600, height = 1200, res = 200)
  
  print(upset(input_df, 
              nsets = ncol(input_df), 
              order.by = "freq", 
              main.bar.color = "steelblue",
              sets.bar.color = "darkred",
              nintersects = 40,
              text.scale = c(1.5, 1.5, 1.2, 1.2, 1.8, 1.5)))
  
  dev.off()
  message("UpSet plot saved to results/plots/upset_clusters.png")
  
} else {
  message("Not enough clusters with significant genes to create an UpSet plot.")
}
