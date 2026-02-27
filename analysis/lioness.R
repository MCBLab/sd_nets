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

# --- Revised Custom LIONESS function (Diego + Optimization) ---
lioness_cor_optimized <- function(se) {
  X <- assay(se)
  n <- ncol(X); p <- nrow(X)
  S <- rowSums(X); SS <- rowSums(X^2); SP <- X %*% t(X) 
  
  num_full <- SP - outer(S, S) / n
  den_full <- sqrt((SS - S^2 / n) %o% (SS - S^2 / n))
  C_full <- num_full / den_full
  
  ut_idx <- which(upper.tri(C_full), arr.ind = TRUE)
  genes <- rownames(X)
  
  # Keep these as separate vectors to avoid strsplit later
  reg_names <- genes[ut_idx[, 1]]
  tar_names <- genes[ut_idx[, 2]]
  
  v_full <- C_full[ut_idx]
  out <- matrix(NA_real_, nrow = length(reg_names), ncol = n)
  colnames(out) <- colnames(X)
  
  for (q in seq_len(n)) {
    cat("Sample", q, "/", n, "\n")
    xq <- X[, q]
    S_q <- S - xq; SS_q <- SS - xq^2; SP_q <- SP - xq %*% t(xq)
    num_q <- SP_q - outer(S_q, S_q) / (n - 1)
    den_q <- sqrt((SS_q - S_q^2 / (n - 1)) %o% (SS_q - S_q^2 / (n - 1)))
    C_q <- num_q / den_q
    out[, q] <- n * (v_full - C_q[ut_idx]) + C_q[ut_idx]
  }
  return(list(W = out, reg = reg_names, tar = tar_names))
}

cat("Running LIONESS...\n")
lioness_out <- lioness_cor_optimized(se_filtered)
W <- lioness_out$W
reg_names <- lioness_out$reg
tar_names <- lioness_out$tar
rm(lioness_out)
gc()

# --- 4.5. Local STRING Filtering (HPC Optimized) ---
cat("Loading local STRING files and mapping IDs...\n")

# 1. Load Aliases with specific filtering
# We need to map 9606.ENSP... to ENSG...
aliases <- vroom("data/9606.protein.aliases.v12.0.txt.gz", 
                 comment = "#", col_names = c("string_id", "alias", "source"),
                 col_types = "ccc")

cat("Creating ENSP to ENSG dictionary...\n")
ensp_to_ensg_df <- aliases %>%
  filter(grepl("Ensembl_gene", source)) %>% # Strictly target Gene ID sources
  mutate(string_id = gsub("9606\\.", "", string_id)) %>% # Strip 9606. prefix
  select(string_id, alias) %>%
  distinct(string_id, .keep_all = TRUE)

ensp_to_ensg <- ensp_to_ensg_df$alias
names(ensp_to_ensg) <- ensp_to_ensg_df$string_id

rm(aliases, ensp_to_ensg_df); gc()

# 2. Load STRING links
cat("Loading experimental links...\n")
links <- vroom("data/9606.protein.links.detailed.v12.0.txt.gz", delim = " ") %>%
  filter(experimental > 0) %>%
  mutate(protein1 = gsub("9606\\.", "", protein1),
         protein2 = gsub("9606\\.", "", protein2)) %>%
  select(protein1, protein2)

cat("Mapping STRING protein IDs to gene IDs...\n")
links$gene1 <- ensp_to_ensg[links$protein1]
links$gene2 <- ensp_to_ensg[links$protein2]

# Remove interactions that didn't map to an ENSG ID
links <- links %>% filter(!is.na(gene1) & !is.na(gene2))

# Create keys for STRING edges
cat("Building lookup keys for experimental edges...\n")
s_keys <- unique(paste(pmin(links$gene1, links$gene2), 
                       pmax(links$gene1, links$gene2), sep = "_"))

rm(links, ensp_to_ensg); gc()

# 3. Filter LIONESS Matrix W
cat("Matching LIONESS edges against experimental data...\n")
# Ensure LIONESS keys are also sorted
l_keys <- paste(pmin(reg_names, tar_names), pmax(reg_names, tar_names), sep = "_")

valid_mask <- l_keys %in% s_keys
cat("Found", sum(valid_mask), "matching experimental edges.\n")

if (sum(valid_mask) == 0) {
  # Sanity check if it still fails
  cat("LIONESS ID Example:", reg_names[1], "\n")
  cat("STRING mapped ID Example:", s_keys[1], "\n")
  stop("Error: Zero edges matched. Please check ID formats printed above.")
}

# Apply Filter
W <- W[valid_mask, , drop = FALSE]
reg_all <- reg_names[valid_mask]
tar_all <- tar_names[valid_mask]
rownames(W) <- paste(reg_all, tar_all, sep = "_")

rm(l_keys, s_keys, valid_mask); gc()

# --- 5. Build Degree Matrix (HPC Optimized) ---
cat("Starting Degree Matrix construction...\n")
all_genes_in_network <- unique(c(reg_all, tar_all))

deg_list <- mclapply(seq_len(ncol(W)), function(i) {
  col_w <- abs(W[, i])
  thr <- mean(col_w, na.rm = TRUE) + (2 * sd(col_w, na.rm = TRUE))
  keep <- !is.na(col_w) & col_w > thr
  
  # Use the pre-existing vectors
  active_genes <- c(reg_all[keep], tar_all[keep])
  gene_counts <- table(active_genes)
  
  sample_deg <- setNames(integer(length(all_genes_in_network)), all_genes_in_network)
  sample_deg[names(gene_counts)] <- as.integer(gene_counts)
  return(sample_deg)
}, mc.cores = 32)

deg_mat <- as.data.frame(do.call(cbind, deg_list))
colnames(deg_mat) <- colnames(W)
rownames(deg_mat) <- all_genes_in_network
saveRDS(deg_mat, "results/lioness_gene_degree_matrix.rds")

# Cleanup temporary objects to free up RAM before moving to t-tests
rm(edge_parts, reg_all, tar_all, deg_list)
gc()

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

write.csv(sig_edges, "results/lioness_significant_edges.csv", row.names = FALSE)

saveRDS(final_results, "results/lioness_full_differential_connectivity_results.rds")

# Uncomment the next line only if you need the individual sample weights for every edge
# saveRDS(W, "lioness_weight_matrix_full.rds")

cat("Analysis complete. Results saved.\n")
cat("Total significant edges found:", nrow(sig_edges), "\n")

rm(W)
gc()

