library(dplyr)
library(tibble)
library(parallel)
library(vroom)
library(biomaRt)

# Load data
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
n_cores <- 6

# Helper function to translate IDs
translate_ids <- function(df) {
  if (nrow(df) == 0) return(df)
  
  # Strip version numbers from Ensembl IDs if they exist (e.g., .13)
  clean_ids <- gsub("\\..*$", "", df$gene)
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- tryCatch({
    getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "description"),
      filters = "ensembl_gene_id",
      values = clean_ids,
      mart = ensembl
    )
  }, error = function(e) {
    message("BiomaRt error: ", e$message)
    return(data.frame(ensembl_gene_id = character(), external_gene_name = character(), description = character()))
  })
  
  df$clean_id <- clean_ids
  df <- df %>%
    left_join(mapping, by = c("clean_id" = "ensembl_gene_id")) %>%
    rename(symbol = external_gene_name) %>%
    dplyr::select(-clean_id)
  
  return(df)
}

# Core analysis function
run_dig_analysis <- function(deg_mat, metadata, formula_str, term_to_test, output_prefix) {
  # Sample alignment
  common_samples <- intersect(colnames(deg_mat), metadata$ID)
  local_deg_mat <- deg_mat[, common_samples]
  local_metadata <- metadata %>% filter(ID %in% common_samples) %>% arrange(match(ID, common_samples))
  
  if(!all(colnames(local_deg_mat) == local_metadata$ID)) stop("Sample alignment failed!")
  
  message("Running comparison: ", output_prefix, " with ", length(common_samples), " samples.")
  
  # Create a data frame for each gene model
  results_list <- mclapply(seq_len(nrow(local_deg_mat)), function(i) {
    gene_name <- rownames(local_deg_mat)[i]
    deg_vec <- as.numeric(local_deg_mat[i, ])
    
    # Combine degree with metadata
    dat <- cbind(data.frame(degree = deg_vec), local_metadata)
    
    if (var(dat$degree, na.rm = TRUE) == 0) return(NULL)
    
    res <- tryCatch({
      fit <- lm(as.formula(formula_str), data = dat)
      summ <- summary(fit)
      coefs <- summ$coefficients
      
      if (!(term_to_test %in% rownames(coefs))) return(NULL)
      
      var_name <- NA
      target_level <- NA
      # Basic heuristic for common cases:
      if (grepl("Condition", term_to_test)) { var_name <- "Condition"; target_level <- "Sjogrens" }
      if (grepl("AcquisitionGroup", term_to_test)) { var_name <- "AcquisitionGroup"; target_level <- "Late" }
      if (grepl("TargetCluster", term_to_test)) { var_name <- "TargetCluster"; target_level <- "Yes" }
      
      mean_target <- NA
      mean_ref <- NA
      if (!is.na(var_name)) {
        mean_target <- mean(dat$degree[dat[[var_name]] == target_level], na.rm = TRUE)
        mean_ref <- mean(dat$degree[dat[[var_name]] != target_level], na.rm = TRUE)
      }

      data.frame(
        gene = gene_name,
        mean_target = mean_target,
        mean_ref = mean_ref,
        estimate = coefs[term_to_test, "Estimate"],
        p_val = coefs[term_to_test, "Pr(>|t|)"],
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
  
  # Save all results
  dir.create("results/digs", showWarnings = FALSE, recursive = TRUE)
  write.csv(results, paste0("results/digs/", output_prefix, "_all_genes.csv"), row.names = FALSE)
  
  # Save significant results with translation
  sigs <- results %>% filter(padj < 0.05)
  if (nrow(sigs) > 0) {
    sigs <- translate_ids(sigs)
    write.csv(sigs, paste0("results/digs/", output_prefix, "_significant_padj05.csv"), row.names = FALSE)
  }
  
  message("Finished ", output_prefix, ". Found ", nrow(sigs), " significant genes.")
  
  return(results)
}

# --- COMPARISON 1: Sjogrens vs Control (Original) ---
metadata_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv")
metadata_sjs_ctrl$Condition <- factor(metadata_sjs_ctrl$Condition, levels = c("Control", "Sjogrens"))
metadata_sjs_ctrl$Sex <- as.factor(metadata_sjs_ctrl$Sex)

run_dig_analysis(
  deg_mat, 
  metadata_sjs_ctrl, 
  "degree ~ Condition + Sex", 
  "ConditionSjogrens", 
  "sjs_vs_ctrl"
)

# --- COMPARISON 2: Late vs Early Acquisition (Women only) ---
metadata_detailed <- vroom("data/precisesads/metadata_sjs_detailed.csv")
metadata_women <- metadata_detailed %>%
  filter(Sex == "Female") %>%
  mutate(AcquisitionAge = Age - DiseaseDuration_Years) %>%
  mutate(AcquisitionGroup = ifelse(AcquisitionAge > 48, "Late", "Early")) %>%
  mutate(AcquisitionGroup = factor(AcquisitionGroup, levels = c("Early", "Late")))

run_dig_analysis(
  deg_mat, 
  metadata_women, 
  "degree ~ AcquisitionGroup", 
  "AcquisitionGroupLate", 
  "late_vs_early_acquisition"
)

# --- COMPARISON 3: Sjogrens Biomolecular Clusters ---
clusters <- vroom("data/precisesads/sjogren_clusters.csv")
# Join clusters with metadata_detailed for Sex covariate if needed
metadata_clusters <- metadata_detailed %>%
  inner_join(clusters, by = "ID") %>%
  mutate(PREDICTION = as.factor(PREDICTION))

for (cl in levels(metadata_clusters$PREDICTION)) {
  current_meta <- metadata_clusters %>%
    mutate(TargetCluster = factor(ifelse(PREDICTION == cl, "Yes", "No"), levels = c("No", "Yes")))
  
  # Include Sex as covariate if it has more than 1 level in this subset
  formula_cl <- "degree ~ TargetCluster"
  if (length(unique(na.omit(current_meta$Sex))) > 1) {
    formula_cl <- "degree ~ TargetCluster + Sex"
  }
  
  run_dig_analysis(
    deg_mat, 
    current_meta, 
    formula_cl, 
    "TargetClusterYes", 
    paste0("cluster_", cl, "_vs_others")
  )
}
