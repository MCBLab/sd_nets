library(dplyr)
library(tibble)
library(parallel)
library(vroom)
library(biomaRt)
library(car) # Para o teste de Levene
library(org.Hs.eg.db)

# Carregar dados
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
n_cores <- 6


translate_ids <- function(df) {
  if (nrow(df) == 0) return(df)
  
  # Limpar IDs
  clean_ids <- gsub("\\..*$", "", df$gene)
  
  # Mapeamento local (muito mais rápido e offline)
  mapping <- data.frame(
    ensembl_gene_id = clean_ids,
    symbol = mapIds(org.Hs.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
    description = mapIds(org.Hs.eg.db, keys = clean_ids, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first"),
    stringsAsFactors = FALSE
  )
  
  df$clean_id <- clean_ids
  df <- df %>%
    left_join(mapping, by = c("clean_id" = "ensembl_gene_id")) %>%
    dplyr::select(-clean_id)
  
  return(df)
}

# Função de análise baseada em Teste-t e Levene
run_dig_analysis_ttest <- function(deg_mat, metadata, group_var, target_level, output_prefix) {
  # Alinhamento de amostras
  common_samples <- intersect(colnames(deg_mat), metadata$ID)
  local_deg_mat <- deg_mat[, common_samples]
  local_metadata <- metadata %>% 
    filter(ID %in% common_samples) %>% 
    arrange(match(ID, common_samples))
  
  if(!all(colnames(local_deg_mat) == local_metadata$ID)) stop("Falha no alinhamento!")
  
  message("Executando Teste-t: ", output_prefix, " com ", length(common_samples), " amostras.")
  
  # Preparar fator de agrupamento
  group_factor <- factor(local_metadata[[group_var]])
  
  results_list <- mclapply(seq_len(nrow(local_deg_mat)), function(i) {
    gene_name <- rownames(local_deg_mat)[i]
    deg_vec <- as.numeric(local_deg_mat[i, ])
    
    # Remover NAs para o teste
    valid_idx <- !is.na(deg_vec) & !is.na(group_factor)
    v_deg <- deg_vec[valid_idx]
    v_group <- group_factor[valid_idx]
    
    if (length(unique(v_group)) < 2 || var(v_deg) == 0) return(NULL)
    
    res <- tryCatch({
      # 1. Teste de Levene para homogeneidade de variância
      # Se p < 0.05, variâncias são desiguais
      lev_test <- car::leveneTest(v_deg ~ v_group)
      var_equal <- lev_test[["Pr(>F)"]][1] > 0.05
      
      # 2. Teste-t de duas caudas
      # var.equal = TRUE -> Teste-t de Student clássico
      # var.equal = FALSE -> Teste-t de Welch (ajustado para variâncias desiguais)
      ttest_res <- t.test(v_deg ~ v_group, var.equal = var_equal)
      
      mean_target <- mean(v_deg[v_group == target_level])
      mean_ref <- mean(v_deg[v_group != target_level])
      
      data.frame(
        gene = gene_name,
        mean_target = mean_target,
        mean_ref = mean_ref,
        estimate = mean_target - mean_ref, # Diferença de médias
        levene_p = lev_test[["Pr(>F)"]][1],
        variances_equal = var_equal,
        p_val = ttest_res$p.value,
        stringsAsFactors = FALSE
      )
    }, error = function(e) return(NULL))
    return(res)
  }, mc.cores = n_cores)
  
  results <- do.call(rbind, results_list)
  if (is.null(results) || nrow(results) == 0) return(NULL)
  
  results <- results %>%
    mutate(padj = p.adjust(p_val, method = "BH")) %>%
    arrange(p_val)
  
  dir.create("results/digs", showWarnings = FALSE, recursive = TRUE)
  write.csv(results, paste0("results/digs/", output_prefix, "_all_genes.csv"), row.names = FALSE)
  
  sigs <- results %>% filter(padj < 0.05)
  if (nrow(sigs) > 0) {
    sigs <- translate_ids(sigs)
    write.csv(sigs, paste0("results/digs/", output_prefix, "_significant_padj05.csv"), row.names = FALSE)
  }
  
  message("Concluído ", output_prefix, ". Encontrados ", nrow(sigs), " genes significativos.")
  return(results)
}

# --- Execução das Comparações ---

# 1. Sjögren vs Controle
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv")
run_dig_analysis_ttest(deg_mat, meta_sjs_ctrl, "Condition", "Sjogrens", "sjs_vs_ctrl")

# 2. Aquisição Tardia vs Precoce (Mulheres)
meta_detailed <- vroom("data/precisesads/metadata_sjs_detailed.csv")
meta_women <- meta_detailed %>%
  filter(Sex == "Female") %>%
  mutate(AcquisitionAge = Age - DiseaseDuration_Years,
         AcquisitionGroup = ifelse(AcquisitionAge > 48, "Late", "Early"))

run_dig_analysis_ttest(deg_mat, meta_women, "AcquisitionGroup", "Late", "late_vs_early_acquisition")

# 3. Clusters Moleculares (Um vs Todos)
clusters <- vroom("data/precisesads/sjogren_clusters.csv")
meta_clusters <- meta_detailed %>% inner_join(clusters, by = "ID")

for (cl in unique(meta_clusters$PREDICTION)) {
  current_meta <- meta_clusters %>%
    mutate(TargetCluster = ifelse(PREDICTION == cl, "Yes", "No"))
  
  run_dig_analysis_ttest(
    deg_mat, 
    current_meta, 
    "TargetCluster", 
    "Yes", 
    paste0("cluster_", cl, "_vs_others")
  )
}