library(dplyr)
library(ggplot2)
library(vroom)
library(readr)
library(tibble)
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
    x = "log2(Mean Degree Ratio)",
    y = "-log10(p-value)"
  )

ggsave("results/plots/digs/volcano_sjs_vs_ctrl.svg", p_volcano, width = 8, height = 7, bg = "white")

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
        x = "log2(Mean Degree Ratio)",
        y = "-log10(p-value)"
      )

    ggsave(paste0("results/plots/digs/volcano_cluster_", cl, ".svg"), p_volcano, width = 8, height = 6, bg = "white")
  }
}
