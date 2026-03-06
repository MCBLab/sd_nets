library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Create output directories
dir.create("results/enrichment", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots/enrichment", showWarnings = FALSE, recursive = TRUE)

# Function to run GO enrichment for a given set of genes
run_go_enrichment <- function(sig_file, output_prefix) {
  message("Running GO enrichment for: ", output_prefix)
  
  df <- tryCatch({
    read.csv(sig_file)
  }, error = function(e) {
    message("Could not read file: ", sig_file)
    return(NULL)
  })
  
  if (is.null(df) || nrow(df) < 5) {
    message("Too few genes (n < 5) for enrichment in: ", output_prefix)
    return(NULL)
  }
  
  my_digs <- df$gene
  my_digs_clean <- gsub("\\..*$", "", my_digs)
  
  ego <- tryCatch({
    enrichGO(gene          = my_digs_clean,
             OrgDb         = org.Hs.eg.db,
             keyType       = 'ENSEMBL',
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05,
             readable      = TRUE)
  }, error = function(e) {
    message("Error in enrichGO for ", output_prefix, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No significant GO enrichment found for: ", output_prefix)
    return(NULL)
  }
  
  # Save CSV results
  write.csv(as.data.frame(ego), 
            paste0("results/enrichment/GO_", output_prefix, ".csv"), 
            row.names = FALSE)
  
  # Generate and save barplot
  p <- barplot(ego, showCategory=10, font.size = 10) +
    theme_bw() +
    labs(title = paste("GO Biological Processes:", gsub("_", " ", output_prefix)),
         x = "Gene Count")
  
  ggsave(paste0("results/plots/enrichment/GO_barplot_", output_prefix, ".png"), 
         p, width = 8, height = 6)
  
  return(ego)
}

# Identify all significant DIG files in results/digs/
sig_files <- list.files("results/digs", pattern = "_significant_padj05.csv$", full.names = TRUE)

# Loop through all found significant result files
for (f in sig_files) {
  prefix <- gsub("_significant_padj05.csv", "", basename(f))
  run_go_enrichment(f, prefix)
}

# --- The following parts (KEGG/Metabolic) are kept for the main comparison only ---
# We'll use the sjs_vs_ctrl_significant_padj05.csv result as the primary one
primary_sig_file <- "results/digs/sjs_vs_ctrl_significant_padj05.csv"

if (file.exists(primary_sig_file)) {
  message("Performing KEGG and Metabolic analysis for primary result (sjs_vs_ctrl)...")
  
  df_primary <- read.csv(primary_sig_file)
  my_digs_clean <- gsub("\\..*$", "", df_primary$gene)
  
  # For KEGG, we need Entrez IDs
  entrez_ids <- tryCatch({
    bitr(my_digs_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  }, error = function(e) return(NULL))
  
  if (!is.null(entrez_ids)) {
    ekegg <- tryCatch({
      enrichKEGG(gene         = entrez_ids$ENTREZID,
                 organism     = 'hsa', 
                 pvalueCutoff = 0.05)
    }, error = function(e) return(NULL))
    
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      write.csv(as.data.frame(ekegg), "results/enrichment/KEGG_sjs_vs_ctrl.csv", row.names = FALSE)
      
      p1f <- cnetplot(setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID'), 
                      circular = TRUE, 
                      colorEdge = TRUE, 
                      showCategory = 5) +
        labs(title = "KEGG Pathway-Gene Network: sjs_vs_ctrl")
      
      ggsave("results/plots/enrichment/KEGG_cnet_sjs_vs_ctrl.png", p1f, width = 8, height = 6)
    }
  }
}
