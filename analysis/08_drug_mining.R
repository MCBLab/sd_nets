library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(RSQLite)
library(org.Hs.eg.db)
library(vroom)

query_drugs_by_uniprot <- function(uniprot_ids, con) {
  # Safety: wrap each ID in quotes
  idvec <- paste0("('", paste(uniprot_ids, collapse = "', '"), "')")

  query <- paste0("
    SELECT d.molregno,
           e.chembl_id AS drug_chembl_id,
           e.pref_name AS drug_name,
           d.mechanism_of_action AS moa,
           d.action_type AS action_type,
           e.first_approval AS first_approval,
           a.chembl_id AS target_chembl_id,
           d.tid AS tid,
           c.accession AS uniprot_id,
           c.description AS description,
           c.organism AS organism,
           GROUP_CONCAT(DISTINCT (f.mesh_id || ': ' || f.mesh_heading)) AS mesh_indication

    FROM target_dictionary AS a
    LEFT JOIN target_components AS b ON a.tid = b.tid
    LEFT JOIN component_sequences AS c ON b.component_id = c.component_id
    LEFT JOIN drug_mechanism AS d ON a.tid = d.tid
    LEFT JOIN molecule_dictionary AS e ON d.molregno = e.molregno
    LEFT JOIN drug_indication AS f ON e.molregno = f.molregno

    WHERE c.accession IN ", idvec, "
    GROUP BY d.molregno, c.accession
    ORDER BY c.accession, c.description
  ")

  # Run the query
  result <- dbGetQuery(con, query)
  return(result)
}

deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
ssgsea_scores <- readRDS("results/ssgsea/ssgsea_hallmark_scores.rds")

clusters_meta <- vroom("data/precisesads/sjogren_clusters.csv") %>%
  dplyr::select(ID, Cluster = PREDICTION)

common_samples <- intersect(colnames(deg_mat), colnames(ssgsea_scores)) %>%
  intersect(clusters_meta$ID)

deg_mat <- deg_mat[, common_samples]
ssgsea_scores <- ssgsea_scores[, common_samples]
clusters_meta <- clusters_meta %>% 
  filter(ID %in% common_samples) %>% 
  arrange(match(ID, common_samples))

# ChEMBL
con_chembl <- dbConnect(RSQLite::SQLite(), "data/chembl_36_sqlite/chembl_36.db")

# DrugBank
drugbank_targets <- read_csv("data/drugbank_targets.csv")

mine_cluster_targets <- function(target_cl) {
  message("\n==> Cluster: ", target_cl)
  
  digs_file <- paste0("results/digs/cluster_", target_cl, "_vs_others_significant_padj05.csv")
  path_file <- paste0("results/ssgsea/limma_pathway_cluster_", target_cl, "_vs_others.csv")
  
  if(!file.exists(digs_file) | !file.exists(path_file)) return(NULL)
  
  digs <- read_csv(digs_file)
  significant_pathways <- read_csv(path_file) %>% filter(padj < 0.05)
  
  if(nrow(significant_pathways) == 0) return(NULL)
  
  # Correlation (DIG degree vs Pathway Score)
  available_genes <- intersect(digs$gene, rownames(deg_mat))
  available_paths <- intersect(significant_pathways$pathway, rownames(ssgsea_scores))
  
  if(length(available_genes) == 0 | length(available_paths) == 0) return(NULL)

  message("Correlating ", length(available_genes), " genes with ", length(available_paths), " pathways...")

  # t() because cor() expects variables in columns
  cor_mat <- cor(t(deg_mat[available_genes, , drop=FALSE]), 
                 t(ssgsea_scores[available_paths, , drop=FALSE]), 
                 method = "spearman")
  
  # Convert matrix to long format
  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("gene", "pathway_associated", "rho")
  cor_df$gene <- as.character(cor_df$gene)
  cor_df$pathway_associated <- as.character(cor_df$pathway_associated)
  
  # Add symbols
  candidates <- cor_df %>%
    left_join(dplyr::select(digs, gene, symbol), by = "gene")
  
  if(nrow(candidates) == 0) return(NULL)

  # Map all significant DIGs to UniProt
  uniprot_map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(candidates$gene), 
                                      column = "UNIPROT", keytype = "ENSEMBL") %>%
    filter(!is.na(UNIPROT))
  
  if(nrow(uniprot_map) == 0) return(NULL)

  # Filter candidates to only those with UniProt
  candidates_druggable <- candidates %>%
    inner_join(uniprot_map, by = c("gene" = "ENSEMBL"))

  # Check which of these UniProts are in our drug databases
  druggable_uniprots <- intersect(candidates_druggable$UNIPROT, 
                                  unique(c(drugbank_targets$target_uniprot_id, 
                                           dbGetQuery(con_chembl, "SELECT accession FROM component_sequences")$accession)))
  
  if(length(druggable_uniprots) == 0) return(NULL)
  
  # Final candidates that actually have a drug match possibility
  candidates_final_pre <- candidates_druggable %>%
    filter(UNIPROT %in% druggable_uniprots)

  message("Calculating p-values for ", nrow(candidates_final_pre), " pathway-gene pairs...")

  # Add p-values only for these final candidates
  candidates_final <- candidates_final_pre %>%
    rowwise() %>%
    mutate(p_val = tryCatch({
      cor.test(as.numeric(deg_mat[gene, ]), 
               as.numeric(ssgsea_scores[pathway_associated, ]), 
               method = "spearman")$p.value
    }, error = function(e) 1)) %>%
    ungroup()

  # 1. DrugBank
  db_matches <- drugbank_targets %>%
    filter(target_uniprot_id %in% candidates_final$UNIPROT) %>%
    inner_join(candidates_final, by = c("target_uniprot_id" = "UNIPROT")) %>%
    dplyr::select(gene_name = symbol, drug_name, drug_id = drugbank_id, 
           mechanism = action, rho, p_val, pathway_associated) %>%
    mutate(source = "DrugBank")
  
  # 2. ChEMBL
  chembl_matches <- query_drugs_by_uniprot(na.omit(unique(candidates_final$UNIPROT)), con_chembl) %>%
    filter(!is.na(drug_name)) %>%
    inner_join(candidates_final, by = c("uniprot_id" = "UNIPROT")) %>%
    dplyr::select(gene_name = symbol, drug_name, drug_id = drug_chembl_id, 
           mechanism = moa, rho, p_val, pathway_associated) %>%
    mutate(source = "ChEMBL")
  
  final_targets <- bind_rows(db_matches, chembl_matches) %>% 
    distinct(gene_name, drug_name, drug_id, mechanism, pathway_associated, .keep_all = TRUE) %>%
    arrange(desc(abs(rho)))
  
  message(nrow(final_targets))
  
  return(final_targets)
}

dir.create("results/repurposing", showWarnings = FALSE)
all_clusters <- unique(clusters_meta$Cluster)

cluster_results <- map(all_clusters, ~mine_cluster_targets(.x))
names(cluster_results) <- paste0("Cluster_", all_clusters)

saveRDS(cluster_results, "results/repurposing/final_drug_prioritization.rds")

# Save individual CSV files for each cluster
iwalk(cluster_results, function(res, name) {
  if (!is.null(res)) {
    file_name <- paste0("results/repurposing/", tolower(name), "_drug_prioritization.csv")
    write_csv(res, file_name)
    message("Saved: ", file_name)
  }
})

dbDisconnect(con_chembl)

