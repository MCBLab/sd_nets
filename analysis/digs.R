library(dplyr)
library(tibble)
library(parallel)
library(vroom)
library(biomaRt)

deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
metadata <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

common_samples <- intersect(colnames(deg_mat), metadata$ID)
deg_mat <- deg_mat[, common_samples]
metadata <- metadata %>% filter(ID %in% common_samples) %>% arrange(match(ID, common_samples))

if(!all(colnames(deg_mat) == metadata$ID)) stop("Sample alignment failed!")

sjs_idx <- which(metadata$Condition == "Sjogrens")
ctrl_idx <- which(metadata$Condition == "Control")
cat("Comparing", length(sjs_idx), "Sjogren's vs", length(ctrl_idx), "Controls...\n")

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 32))

dig_results_list <- mclapply(seq_len(nrow(deg_mat)), function(i) {
  gene_name <- rownames(deg_mat)[i]
  deg_vec <- as.numeric(deg_mat[i, ])
  
  group_sjs <- deg_vec[sjs_idx]
  group_ctrl <- deg_vec[ctrl_idx]
  
  if (sd(group_sjs) == 0 && sd(group_ctrl) == 0) return(NULL)
  
  res <- tryCatch({
    t_test <- t.test(group_sjs, group_ctrl, alternative = "two.sided")
    data.frame(
      gene = gene_name,
      mean_deg_sjs = mean(group_sjs),
      mean_deg_ctrl = mean(group_ctrl),
      diff_deg = mean(group_sjs) - mean(group_ctrl),
      p_val = t_test$p.value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) return(NULL))
  return(res)
}, mc.cores = n_cores)

dig_results <- do.call(rbind, dig_results_list)

dig_results <- dig_results %>%
  mutate(padj = p.adjust(p_val, method = "BH")) %>%
  arrange(p_val)

DIGs <- dig_results %>% 
  filter(padj < 0.05)

translate_ids <- function(df) {
  if (nrow(df) == 0) return(df)
  
  # Strip version numbers from Ensembl IDs if they exist (e.g., .13)
  clean_ids <- gsub("\\..*$", "", df$gene)
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = clean_ids,
    mart = ensembl
  )
  
  df$clean_id <- clean_ids
  df <- df %>%
    left_join(mapping, by = c("clean_id" = "ensembl_gene_id")) %>%
    rename(symbol = external_gene_name) %>%
    dplyr::select(-clean_id)
  
  return(df)
}

DIGs <- translate_ids(DIGs)

write.csv(dig_results, "results/differential_interactivity_all_genes.csv", row.names = FALSE)
write.csv(DIGs, "results/DIGs_significant_padj05.csv", row.names = FALSE)

cat("Total DIGs found:", nrow(DIGs), "\n")