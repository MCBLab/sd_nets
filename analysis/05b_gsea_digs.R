library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(vroom)

dir.create("results/enrichment/gsea", showWarnings = FALSE, recursive = TRUE)

# GMT files
gmt_files <- list(
  hallmark = "data/h.all.v2026.1.Hs.symbols.gmt",
  curated  = "data/c2.all.v2026.1.Hs.symbols.gmt",
  ontology = "data/c5.all.v2026.1.Hs.symbols.gmt"
)

gmts <- lapply(gmt_files, read.gmt)

# Function to run GSEA for a given DIG result file
run_gsea_for_file <- function(file_path, gmts, output_prefix) {
  message("Processing: ", output_prefix)
  
  # Read all genes
  df <- vroom(file_path, show_col_types = FALSE)
  
  # Clean IDs and Map to Symbols
  df$clean_id <- gsub("[.].*$", "", df$gene)
  
  # Mapping
  mapping <- data.frame(
    ensembl_gene_id = unique(df$clean_id),
    symbol = mapIds(org.Hs.eg.db, keys = unique(df$clean_id), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
    stringsAsFactors = FALSE
  )
  
  df_mapped <- df %>%
    left_join(mapping, by = c("clean_id" = "ensembl_gene_id")) %>%
    filter(!is.na(symbol)) %>%
    filter(!is.na(estimate))
  
  # Handle duplicates (keep one with highest absolute estimate)
  df_ranked <- df_mapped %>%
    arrange(symbol, desc(abs(estimate))) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    ungroup()
  
  # Create ranked list
  gene_list <- df_ranked$estimate
  names(gene_list) <- df_ranked$symbol
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA for each GMT
  for (name in names(gmts)) {
    message("  - Running GSEA for ", name)
    gsea_res <- tryCatch({
      GSEA(gene_list, 
           TERM2GENE = gmts[[name]], 
           pvalueCutoff = 0.05, 
           pAdjustMethod = "BH",
           eps = 0) # Use eps = 0 for better precision in p-values
    }, error = function(e) {
      message("    Error in GSEA (", name, "): ", e$message)
      return(NULL)
    })
    
    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      write.csv(as.data.frame(gsea_res), 
                paste0("results/enrichment/gsea/GSEA_", name, "_", output_prefix, ".csv"), 
                row.names = FALSE)
      message("    Saved ", nrow(as.data.frame(gsea_res)), " significant pathways.")
    } else {
      message("    No significant pathways found for ", name)
    }
  }
}

# Identify all DIG result files
all_files <- list.files("results/digs", pattern = "_all_genes.csv$", full.names = TRUE)

# Loop through files
for (f in all_files) {
  prefix <- gsub("_all_genes.csv", "", basename(f))
  run_gsea_for_file(f, gmts, prefix)
}

message("GSEA analysis completed.")
