library(dplyr)
library(ggplot2)
library(vroom)
library(readr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)

# --- 1. Load Data ---
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

# Ensure output directory exists
dir.create("results/plots/digs", showWarnings = FALSE, recursive = TRUE)

# Helper function to get symbols for all genes in a result set
# (Since the 'all_genes' file doesn't have symbols, but 'significant' does)
add_symbols <- function(df) {
  clean_ids <- gsub("\\..*$", "", df$gene)
  df$symbol <- mapIds(org.Hs.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  return(df)
}

# --- 2. Sjögren vs Control Analysis ---

# Load results
sjs_vs_ctrl_res <- read_csv("results/digs/sjs_vs_ctrl_all_genes.csv") %>%
  add_symbols()

# Volcano Plot
# We'll use log2(Ratio) for the x-axis and highlight top genes
pseudocount <- 0.01
volcano_data <- sjs_vs_ctrl_res %>%
  mutate(
    log2FC = log2((mean_target + pseudocount) / (mean_ref + pseudocount)),
    logP = -log10(p_val),
    Significant = ifelse(padj < 0.05 & abs(log2FC) > 0.5, "Significant", "Not Significant")
  )

# Top 15 genes to label (up and down)
top_up <- volcano_data %>%
  filter(padj < 0.05 & log2FC > 0) %>%
  arrange(p_val) %>%
  head(15)

top_down <- volcano_data %>%
  filter(padj < 0.05 & log2FC < 0) %>%
  arrange(p_val) %>%
  head(15)

top_genes <- bind_rows(top_up, top_down)

p_volcano <- ggplot(volcano_data, aes(x = log2FC, y = logP, color = Significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_text(data = top_genes, aes(label = symbol), vjust = 1.5, color = "black", check_overlap = TRUE) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Sjögren vs Control (Differential Interactivity)",
    x = "log2(Mean Degree Ratio)",
    y = "-log10(p-value)"
  )

ggsave("results/plots/digs/volcano_sjs_vs_ctrl.png", p_volcano, width = 8, height = 7, bg = "white")

# --- 2.1 Cluster Volcano Plots ---
pseudocount <- 0.01
for (cl in 1:4) {
  cl_file <- paste0("results/digs/cluster_", cl, "_vs_others_all_genes.csv")
  if (file.exists(cl_file)) {
    cl_res <- read_csv(cl_file) %>% add_symbols()
    
    volcano_data <- cl_res %>%
      mutate(
        log2FC = log2((mean_target + pseudocount) / (mean_ref + pseudocount)),
        logP = -log10(p_val),
        Significant = ifelse(padj < 0.05 & abs(log2FC) > 0.5, "Significant", "Not Significant")
      )

    top_up <- volcano_data %>%
      filter(padj < 0.05 & log2FC > 0) %>%
      arrange(p_val) %>%
      head(15)

    top_down <- volcano_data %>%
      filter(padj < 0.05 & log2FC < 0) %>%
      arrange(p_val) %>%
      head(15)

    top_genes <- bind_rows(top_up, top_down)

    p_volcano <- ggplot(volcano_data, aes(x = log2FC, y = logP, color = Significant)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_text(data = top_genes, aes(label = symbol), vjust = 1.5, color = "black", check_overlap = TRUE) +
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
      theme_minimal() +
      labs(
        title = paste0("Volcano Plot: Cluster ", cl, " vs Others (Differential Interactivity)"),
        x = "log2(Mean Degree Ratio)",
        y = "-log10(p-value)"
      )

    ggsave(paste0("results/plots/digs/volcano_cluster_", cl, ".png"), p_volcano, width = 8, height = 7, bg = "white")
  }
}

# Heatmap of Top Significant DIGs
# Get top 50 significant genes by p-value
sig_genes <- sjs_vs_ctrl_res %>%
  filter(padj < 0.05) %>%
  arrange(p_val) %>%
  head(50)

# Prepare matrix for heatmap
# Replace Ensembl IDs with Symbols for the matrix
plot_mat_raw <- deg_mat[sig_genes$gene, ]
rownames(plot_mat_raw) <- sig_genes$symbol

# Annotation for heatmap
anno_col <- meta_sjs_ctrl %>%
  dplyr::select(ID, Condition) %>%
  as.data.frame()
rownames(anno_col) <- anno_col$ID
anno_col$ID <- NULL

# Ensure matrix columns match metadata and are in the same order
common_samples <- intersect(colnames(deg_mat), rownames(anno_col))
plot_mat <- plot_mat_raw[, common_samples]
anno_col <- anno_col[common_samples, , drop = FALSE]

# Apply log2 transformation (using same pseudocount as volcano)
plot_mat_log <- log2(plot_mat + pseudocount)

# Standardize rows (Z-score) for better visualization of relative connectivity
plot_mat_scaled <- t(apply(plot_mat_log, 1, scale))
colnames(plot_mat_scaled) <- colnames(plot_mat)

# Save Heatmap
png("results/plots/digs/heatmap_top50_sjs_vs_ctrl.png", width = 1000, height = 1200, res = 150)
pheatmap(
  plot_mat_scaled,
  annotation_col = anno_col,
  show_colnames = FALSE,
  main = "Top 50 DIGs (SJS vs Control)",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  fontsize_row = 8
)
dev.off()

# --- 3. Refined DIGs Signature Heatmap (Top 20 per Cluster) ---

clusters_list <- 1:4
all_res_list <- list()
for (cl in clusters_list) {
  cl_file <- paste0("results/digs/cluster_", cl, "_vs_others_all_genes.csv")
  if (file.exists(cl_file)) {
    all_res_list[[as.character(cl)]] <- read_csv(cl_file)
  }
}

top_genes_by_cluster <- list()
for (cl in names(all_res_list)) {
  top_genes_by_cluster[[cl]] <- all_res_list[[cl]] %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(20) %>%
    pull(gene)
}

union_top_genes <- unique(unlist(top_genes_by_cluster))

sig_matrix_df <- data.frame(gene = union_top_genes) %>%
  add_symbols() %>%
  mutate(rowname = make.unique(ifelse(is.na(symbol), gene, symbol))) %>%
  column_to_rownames("rowname")

for (cl in names(all_res_list)) {
  cl_data <- all_res_list[[cl]] %>%
    dplyr::select(gene, estimate)
  
  # Join simples: se o gene não for sig/presente, ficará NA e depois 0
  sig_matrix_df[[paste0("Cluster_", cl)]] <- cl_data$estimate[match(union_top_genes, cl_data$gene)]
}

sig_matrix_final <- sig_matrix_df %>%
  dplyr::select(starts_with("Cluster_")) %>%
  replace(is.na(.), 0) %>%
  as.matrix()

palette_rev <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
sig_matrix_ordered <- sig_matrix_final[, paste0("Cluster_", 1:4)]

sig_matrix_log <- sign(sig_matrix_ordered) * log1p(abs(sig_matrix_ordered))


limit <- quantile(abs(sig_matrix_log), 0.98) 
sig_matrix_log[sig_matrix_log > limit] <- limit
sig_matrix_log[sig_matrix_log < -limit] <- -limit

palette_rev <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

row_anno <- data.frame(gene = union_top_genes) %>%
  left_join(
    bind_rows(lapply(names(top_genes_by_cluster), function(cl) {
      data.frame(gene = top_genes_by_cluster[[cl]], Cluster_Origem = paste0("Cluster ", cl))
    })),
    by = "gene"
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  add_symbols() %>%
  mutate(rowname = make.unique(ifelse(is.na(symbol), gene, symbol))) %>%
  dplyr::select(rowname, Cluster_Origem) %>%
  column_to_rownames("rowname")

# Definir cores fixas para os clusters na anotação
ann_colors = list(
  Cluster_Origem = c("Cluster 1" = "#E41A1C", "Cluster 2" = "#377EB8", 
                     "Cluster 3" = "#4DAF4A", "Cluster 4" = "#984EA3")
)

row_gaps <- cumsum(table(row_anno$Cluster_Origem)[unique(row_anno$Cluster_Origem)])
row_gaps <- row_gaps[-length(row_gaps)]

png("results/plots/digs/heatmap_signature_annotated.png", width = 1200, height = 1600, res = 150)
pheatmap(
  sig_matrix_log,
  color = palette_rev,
  breaks = seq(-limit, limit, length.out = 101),
  cluster_cols = FALSE,
  cluster_rows = FALSE, 
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  gaps_row = row_gaps,
  main = "DIG Signature Fingerprint\nAnnotated by Cluster of Origin",
  fontsize_row = 7,
  cellwidth = 40,
  angle_col = 0
)
dev.off()
