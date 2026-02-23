library(lionessR)
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(vroom)
library(dplyr)
library(tibble)
library(httr)
library(jsonlite)
library(parallel) # Added for potential speedup

# 1. Load Metadata
metadata <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

# 2. Load Expression Data
exp_cs <- vroom("data/precisesads/Transcriptome/CS_Transcriptome_counts.csv", delim = ";")
exp_i  <- vroom("data/precisesads/Transcriptome/I_Transcriptome_counts.csv", delim = ";", locale = locale(decimal_mark = ","))

exp <- inner_join(exp_cs, exp_i, by = "ID") %>%
  column_to_rownames("ID")

exp <- exp[, colnames(exp) %in% metadata$ID]
metadata <- metadata[metadata$ID %in% colnames(exp), ]
exp <- exp[, metadata$ID]
exp_rounded <- round(exp)

# 3. Create SummarizedExperiment and Normalize
rowData <- DataFrame(row.names = rownames(exp_rounded), gene = rownames(exp_rounded))
colData <- DataFrame(row.names = metadata$ID, sample = metadata$ID, condition = metadata$Condition)

se <- SummarizedExperiment(assays = list(counts = as.matrix(exp_rounded)), 
                           colData = colData, rowData = rowData)

dds <- DESeqDataSet(se, design = ~ condition)
vst_data <- vst(dds, blind = TRUE)
normalized_counts <- assay(vst_data)

# 4. Filter for Protein Coding Genes
keys <- keys(org.Hs.eg.db, keytype="ENSEMBL")
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys=keys, columns=c("ENSEMBL", "GENETYPE"), keytype="ENSEMBL")
protein_coding_genes <- gene_info$ENSEMBL[gene_info$GENETYPE == "protein-coding"]
protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes)]

valid_genes <- intersect(rownames(normalized_counts), protein_coding_genes)
se_filtered <- vst_data[valid_genes, ]

# Custom LIONESS function from Diego
lioness_cor <- function(se) {
  X <- assay(se)
  n <- ncol(X)
  p <- nrow(X)
  
  S  <- rowSums(X)
  SS <- rowSums(X^2)
  SP <- X %*% t(X) 
  
  num_full <- SP - outer(S, S) / n
  den_full <- sqrt((SS - S^2 / n) %o% (SS - S^2 / n))
  C_full <- num_full / den_full
  
  ut_idx <- which(upper.tri(C_full), arr.ind = TRUE)
  edge_names <- paste(rownames(X)[ut_idx[, 1]], rownames(X)[ut_idx[, 2]], sep = "_")
  m <- length(edge_names)
  
  out <- matrix(NA_real_, nrow = m, ncol = n, dimnames = list(edge_names, colnames(X)))
  v_full <- C_full[ut_idx]
  
  for (q in seq_len(n)) {
    cat("Processing sample", q, "of", n, "...\n")
    xq <- X[, q]
    S_q  <- S  - xq
    SS_q <- SS - xq^2
    SP_q <- SP - xq %*% t(xq)
    
    num_q <- SP_q - outer(S_q, S_q) / (n - 1)
    den_q <- sqrt((SS_q - S_q^2 / (n - 1)) %o% (SS_q - S_q^2 / (n - 1)))
    C_q <- num_q / den_q
    
    v_loo <- C_q[ut_idx]
    out[, q] <- n * (v_full - v_loo) + v_loo
  }
  return(out)
}

cat("Running LIONESS on full protein-coding set (", length(valid_genes), "genes)...\n")
W <- lioness_cor(se_filtered)

# 5. Correlation Threshold (Mean + 2*SD)
# Optimized thresholding for very large matrices
cat("Applying Mean + 2*SD threshold...\n")
for (i in 1:ncol(W)) {
  w_col <- abs(W[, i])
  thr <- mean(w_col, na.rm=TRUE) + 2 * sd(w_col, na.rm=TRUE)
  W[w_col <= thr, i] <- NA
}

# 6. Differential Connectivity Analysis
sjs_idx <- metadata$Condition == "Sjogrens"
ctrl_idx <- metadata$Condition == "Control"

cat("Running per-edge t-tests...\n")
# Using mclapply if on Linux HPC for parallelization (optional)
# Change mc.cores to your allocated CPU count
pvals <- apply(W, 1, function(x) {
  x_sjs <- x[sjs_idx]
  x_ctrl <- x[ctrl_idx]
  x_sjs <- x_sjs[!is.na(x_sjs)]
  x_ctrl <- x_ctrl[!is.na(x_ctrl)]
  
  if (length(x_sjs) < 2 || length(x_ctrl) < 2) return(NA_real_)
  
  var_sjs <- var(x_sjs)
  var_ctrl <- var(x_ctrl)
  if (is.na(var_sjs) || is.na(var_ctrl) || var_sjs == 0 || var_ctrl == 0) return(NA_real_)
  
  tryCatch({ return(t.test(x_sjs, x_ctrl)$p.value) }, error = function(e) { return(NA_real_) })
})

padj <- p.adjust(pvals, method = "BH")
res_edges <- data.frame(edge = rownames(W), pval = pvals, padj = padj)
sig_edges <- subset(res_edges, !is.na(padj) & padj < 0.05)

cat("Significant edges found:", nrow(sig_edges), "\n")

# 7. Build and Save Results
cat("Organizing results table...\n")

# Combine results into a single dataframe
final_results <- data.frame(
  edge = rownames(W),
  pval = pvals,
  padj = padj,
  stringsAsFactors = FALSE
)

# Extract gene names from the edge strings for easier filtering later
edge_parts <- strsplit(final_results$edge, "_")
final_results$gene1 <- sapply(edge_parts, `[`, 1)
final_results$gene2 <- sapply(edge_parts, `[`, 2)

# Reorder columns for better readability
final_results <- final_results %>% 
  select(edge, gene1, gene2, pval, padj)

# Filter for significance
sig_edges <- final_results %>% 
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)

cat("Saving results to disk...\n")

# A. Save the significant edges only (CSV - typically small/medium)
write.csv(sig_edges, "results/lioness_significant_edges.csv", row.names = FALSE)

# B. Save all p-values and statistics (RDS - compressed and fast for R)
# This includes all 200M+ tested edges
saveRDS(final_results, "results/lioness_full_differential_connectivity_results.rds")

# C. Save the weight matrix W (Optional - Very large file)
# Uncomment the next line only if you need the individual sample weights for every edge
# saveRDS(W, "lioness_weight_matrix_full.rds")

cat("Analysis complete. Results saved.\n")
cat("Total significant edges found:", nrow(sig_edges), "\n")

rm(W)
gc()
