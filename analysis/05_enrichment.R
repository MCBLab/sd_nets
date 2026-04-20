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
  
  df_all <- tryCatch({
    read.csv(sig_file)
  }, error = function(e) {
    message("Could not read file: ", sig_file)
    return(NULL)
  })
  
  if (is.null(df_all) || nrow(df_all) < 1) {
    message("Empty file: ", sig_file)
    return(NULL)
  }
  
  # Split into up and down regulated
  df_up <- df_all %>% filter(estimate > 0)
  df_down <- df_all %>% filter(estimate < 0)
  
  # Run enrichment for up, down and all
  results_list <- list(
    "all" = list(df = df_all, suffix = ""),
    "up" = list(df = df_up, suffix = "_up"),
    "down" = list(df = df_down, suffix = "_down")
  )
  
  for (res_type in names(results_list)) {
    df <- results_list[[res_type]]$df
    suffix <- results_list[[res_type]]$suffix
    current_prefix <- paste0(output_prefix, suffix)
    
    if (nrow(df) < 5) {
      message("Too few genes (n < 5) for ", res_type, " enrichment in: ", output_prefix)
      next
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
      message("Error in enrichGO for ", current_prefix, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      message("No significant GO enrichment found for: ", current_prefix)
      next
    }
    
    # Save CSV results
    write.csv(as.data.frame(ego), 
              paste0("results/enrichment/GO_", current_prefix, ".csv"), 
              row.names = FALSE)
    
    # Generate and save barplot
    p <- barplot(ego, showCategory=10, font.size = 10) +
      theme_bw() +
      labs(title = paste("GO Biological Processes:", gsub("_", " ", current_prefix)),
           subtitle = paste(res_type, "regulated DIGs"),
           x = "Gene Count")
    
    ggsave(paste0("results/plots/enrichment/GO_barplot_", current_prefix, ".png"), 
           p, width = 8, height = 6)
  }
}

# Identify all significant DIG files in results/digs/
sig_files <- list.files("results/digs", pattern = "_significant_padj05.csv$", full.names = TRUE)

# Loop through all found significant result files
for (f in sig_files) {
  prefix <- gsub("_significant_padj05.csv", "", basename(f))
  run_go_enrichment(f, prefix)
}