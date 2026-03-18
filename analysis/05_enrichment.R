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