library(vroom)
library(dplyr)
library(tibble)
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)

metadata <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

exp_cs <- vroom("data/precisesads/Transcriptome/CS_Transcriptome_counts.csv", delim = ";")
exp_i  <- vroom("data/precisesads/Transcriptome/I_Transcriptome_counts.csv", delim = ";", locale = locale(decimal_mark = ","))

exp <- inner_join(exp_cs, exp_i, by = "ID") %>%
  column_to_rownames("ID")

exp <- exp[, colnames(exp) %in% metadata$ID]
metadata <- metadata[metadata$ID %in% colnames(exp), ]
exp <- exp[, metadata$ID]
exp_rounded <- round(exp)

rowData <- DataFrame(row.names = rownames(exp_rounded), gene = rownames(exp_rounded))
colData <- DataFrame(row.names = metadata$ID, sample = metadata$ID, condition = metadata$Condition, sex = metadata$Sex)

se <- SummarizedExperiment(assays = list(counts = as.matrix(exp_rounded)), 
                           colData = colData, rowData = rowData)

dds <- DESeqDataSet(se, design = ~ condition + sex)
vst_data <- vst(dds, blind = TRUE)
norm_mat <- assay(vst_data)

ensembl_ids <- gsub("\\..*$", "", rownames(norm_mat))

mapping <- AnnotationDbi::select(org.Hs.eg.db, 
                                 keys = ensembl_ids, 
                                 columns = "SYMBOL", 
                                 keytype = "ENSEMBL")

map_df <- data.frame(ensembl_ver = rownames(norm_mat), 
                     ensembl_clean = ensembl_ids, 
                     stringsAsFactors = FALSE) %>%
  left_join(mapping, by = c("ensembl_clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL))

norm_mat_filtered <- norm_mat[map_df$ensembl_ver, ]
rownames(norm_mat_filtered) <- map_df$SYMBOL

# Handle Duplicate Symbols (Keep the one with highest mean expression)
if(any(duplicated(rownames(norm_mat_filtered)))) {
  gene_means <- rowMeans(norm_mat_filtered)
  norm_mat_filtered <- norm_mat_filtered[order(gene_means, decreasing = TRUE), ]
  norm_mat_filtered <- norm_mat_filtered[!duplicated(rownames(norm_mat_filtered)), ]
}

write.csv(norm_mat_filtered, "results/vst_normalized_counts_symbols.csv")
saveRDS(norm_mat_filtered, "results/vst_normalized_counts_symbols.rds")
