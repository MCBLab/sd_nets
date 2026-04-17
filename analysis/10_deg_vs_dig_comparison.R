library(DESeq2)
library(vroom)
library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(UpSetR)
library(ggplot2)

# --- 1. Differential Expression Analysis (SjS vs Control) ---

metadata <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

exp_cs <- vroom("data/precisesads/Transcriptome/CS_Transcriptome_counts.csv", delim = ";")
exp_i  <- vroom("data/precisesads/Transcriptome/I_Transcriptome_counts.csv", delim = ";", locale = locale(decimal_mark = ","))

exp <- inner_join(exp_cs, exp_i, by = "ID") %>%
  column_to_rownames("ID")

exp <- exp[, colnames(exp) %in% metadata$ID]
metadata <- metadata[metadata$ID %in% colnames(exp), ]
exp <- exp[, metadata$ID]
exp_rounded <- round(exp)

colData <- data.frame(row.names = metadata$ID, condition = factor(metadata$Condition, levels = c("Control", "Sjogrens")))
dds <- DESeqDataSetFromMatrix(countData = exp_rounded, colData = colData, design = ~ condition)

keys <- keys(org.Hs.eg.db, keytype="ENSEMBL")
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys=keys, columns=c("ENSEMBL", "GENETYPE"), keytype="ENSEMBL")
protein_coding_genes <- gene_info$ENSEMBL[gene_info$GENETYPE == "protein-coding"]
protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes)]

keep <- rownames(dds) %in% protein_coding_genes
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, name="condition_Sjogrens_vs_Control")

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(clean_id = gsub("\\..*$", "", gene))

res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$clean_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

write.csv(res_df, "results/deg_sjs_vs_ctrl_all_genes.csv", row.names = FALSE)

sig_degs <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5)

write.csv(sig_degs, "results/deg_sjs_vs_ctrl_significant_padj05.csv", row.names = FALSE)
cat("DE analysis complete. Significant DEGs found:", nrow(sig_degs), "\n")


# --- 2. Intersection Analysis (DEGs vs DIGs) ---

digs <- vroom("results/digs/sjs_vs_ctrl_significant_padj05.csv")
dig_genes <- digs$gene
deg_genes <- sig_degs$clean_id

list_input <- list(DIGs = dig_genes, DEGs = deg_genes)

intersection_genes <- intersect(dig_genes, deg_genes)
cat("Number of DIGs:", length(dig_genes), "\n")
cat("Number of DEGs:", length(deg_genes), "\n")
cat("Number of overlapping genes:", length(intersection_genes), "\n")

intersection_df <- sig_degs %>%
  filter(clean_id %in% intersection_genes) %>%
  dplyr::select(clean_id, symbol, log2FoldChange, padj) %>%
  rename(gene = clean_id, padj_deg = padj) %>%
  left_join(digs %>% dplyr::select(gene, padj) %>% rename(padj_dig = padj), by = "gene")

write.csv(intersection_df, "results/intersection_deg_dig_sjs_vs_ctrl.csv", row.names = FALSE)

# Plot UpSetR
pdf("results/plots/upset_deg_dig_intersection.pdf", width = 8, height = 6)
upset(fromList(list_input), order.by = "freq", main.bar.color = "steelblue", sets.bar.color = "darkred")
dev.off()
cat("Intersection analysis complete. Plot saved to results/plots/upset_deg_dig_intersection.pdf\n")


# --- 3. Filter DIGs that are not DEGs ---

digs_not_degs <- digs %>%
  filter(!(gene %in% deg_genes))

# Save the results
write.csv(digs_not_degs, "results/digs_not_degs_sjs_vs_ctrl.csv", row.names = FALSE)

cat("Number of DIGs that are not DEGs:", nrow(digs_not_degs), "\n")
cat("Table saved to: results/digs_not_degs_sjs_vs_ctrl.csv\n")
