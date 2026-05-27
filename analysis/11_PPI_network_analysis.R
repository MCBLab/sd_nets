library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(vroom)
library(ggraph)
library(RColorBrewer)
library(tidyr)

dir.create("results/ppi_network", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots/ppi_network", showWarnings = FALSE, recursive = TRUE)

digs_file <- "results/digs/sjs_vs_ctrl_significant_padj05.csv"
if (!file.exists(digs_file)) {
  stop("Significant DIGs file not found. Please run analysis/02_digs.R first.")
}
sig_digs <- vroom(digs_file)
dig_genes <- unique(sig_digs$gene)
dig_genes_clean <- gsub("\\..*$", "", dig_genes)

aliases <- vroom("data/9606.protein.aliases.v12.0.txt.gz",
  comment = "#", col_names = c("string_id", "alias", "source"),
  col_types = "ccc"
)

ensp_to_ensg_df <- aliases %>%
  filter(grepl("Ensembl_gene", source)) %>%
  mutate(string_id = gsub("9606\\.", "", string_id)) %>%
  dplyr::select(string_id, alias) %>%
  distinct(string_id, .keep_all = TRUE)

ensp_to_ensg <- ensp_to_ensg_df$alias
names(ensp_to_ensg) <- ensp_to_ensg_df$string_id
rm(aliases, ensp_to_ensg_df)

links <- vroom("data/9606.protein.links.detailed.v12.0.txt.gz", delim = " ") %>%
  filter(experimental > 0) %>%
  mutate(
    protein1 = gsub("9606\\.", "", protein1),
    protein2 = gsub("9606\\.", "", protein2)
  ) %>%
  dplyr::select(protein1, protein2, experimental)

links$gene1 <- ensp_to_ensg[links$protein1]
links$gene2 <- ensp_to_ensg[links$protein2]

ppi_validated <- links %>%
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  filter(gene1 %in% dig_genes_clean & gene2 %in% dig_genes_clean) %>%
  dplyr::select(gene1, gene2, experimental)

rm(links, ensp_to_ensg)
gc()

if (nrow(ppi_validated) == 0) {
  stop("No PPI-validated interactions found between the provided DIGs.")
}

ppi_graph <- graph_from_data_frame(d = ppi_validated, directed = FALSE)
E(ppi_graph)$weight <- ppi_validated$experimental
ppi_graph <- igraph::simplify(ppi_graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "max")

set.seed(42)
clusters_raw <- cluster_louvain(ppi_graph, resolution = 1)
V(ppi_graph)$cluster_raw <- membership(clusters_raw)

cluster_sizes <- table(V(ppi_graph)$cluster_raw)
large_clusters <- names(cluster_sizes[cluster_sizes >= 20])

V(ppi_graph)$symbol <- mapIds(org.Hs.eg.db,
  keys = V(ppi_graph)$name,
  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
)
V(ppi_graph)$estimate <- sig_digs$estimate[match(V(ppi_graph)$name, gsub("\\..*$", "", sig_digs$gene))]

cluster_info_list <- list()
enrichment_results_list <- list()

for (cl_id in large_clusters) {
  cluster_genes <- V(ppi_graph)$name[V(ppi_graph)$cluster_raw == cl_id]

  ego <- tryCatch(
    {
      enrichGO(
        gene = cluster_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        readable = TRUE
      )
    },
    error = function(e) NULL
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    res_df <- as.data.frame(ego) |>
      arrange(desc(Count))

    # Store enrichment results in list
    enrichment_results_list[[as.character(cl_id)]] <- res_df

    # Save raw cluster enrichment
    write.csv(res_df, paste0("results/ppi_network/GO_enrichment_raw_cluster_", cl_id, ".csv"), row.names = FALSE)

    cluster_info_list[[as.character(cl_id)]] <- data.frame(
      raw_cluster = cl_id,
      top_term = res_df$Description[1], # Use the top term as the "General Process"
      size = length(cluster_genes)
    )
  }
}

cluster_summary <- do.call(rbind, cluster_info_list)

# Collapse clusters sharing the same top_term
collapsed_mapping <- cluster_summary %>%
  group_by(top_term) %>%
  mutate(collapsed_name = paste0(top_term, " (", paste(raw_cluster, collapse = ","), ")")) %>%
  ungroup()

# Update node attributes with collapsed cluster names
V(ppi_graph)$general_process <- collapsed_mapping$top_term[match(V(ppi_graph)$cluster_raw, collapsed_mapping$raw_cluster)]

# Final summary of collapsed clusters
final_summary <- collapsed_mapping %>%
  group_by(top_term) %>%
  summarise(
    original_clusters = paste(raw_cluster, collapse = ", "),
    total_size = sum(size),
    .groups = "drop"
  )

write.csv(final_summary, "results/ppi_network/collapsed_cluster_summary.csv", row.names = FALSE)

library(dplyr)
library(graphlayouts)

set.seed(42)
# 1. Calculate the base backbone layout
# Keep remains low to focus on the strongest internal cluster edges
bb <- layout_as_backbone(ppi_graph, keep = 0.05)
xy <- bb$xy

E(ppi_graph)$col <- FALSE
E(ppi_graph)$col[bb$backbone] <- TRUE

# 2. THE FIX: Cluster-based coordinate expansion
# We extract the cluster IDs and calculate the mean position of each cluster

layout_df <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  cluster = V(ppi_graph)$cluster_raw
)

# Calculate centroids, excluding NAs to avoid pushing the 'background' noise too far
centroids <- layout_df %>%
  filter(as.character(cluster) %in% large_clusters) %>%
  group_by(cluster) %>%
  summarise(cx = mean(x), cy = mean(y))

# Push nodes away from the center (0,0) based on their cluster's position
# Increase the 'expansion_factor' until the islands separate to your liking
expansion_factor <- 8

for (cl in centroids$cluster) {
  idx <- which(V(ppi_graph)$cluster_raw == cl)
  # Vector from center to centroid
  vec_x <- centroids$cx[centroids$cluster == cl]
  vec_y <- centroids$cy[centroids$cluster == cl]

  # Apply displacement to all nodes in that cluster
  xy[idx, 1] <- xy[idx, 1] + (vec_x * expansion_factor)
  xy[idx, 2] <- xy[idx, 2] + (vec_y * expansion_factor)
}

# 3. Plot with the modified manual coordinates
V(ppi_graph)$degree <- degree(ppi_graph)
V(ppi_graph)$label <- NA

top_nodes <- data.frame(
  id = 1:vcount(ppi_graph),
  degree = V(ppi_graph)$degree,
  cluster = V(ppi_graph)$cluster_raw
) %>%
  filter(as.character(cluster) %in% large_clusters) %>%
  group_by(cluster) %>%
  slice_max(order_by = degree, n = 10, with_ties = FALSE) %>%
  pull(id)

V(ppi_graph)$label[top_nodes] <- V(ppi_graph)$symbol[top_nodes]

p1 <- ggraph(ppi_graph, layout = "manual", x = xy[, 1], y = xy[, 2]) +
  geom_edge_link0(aes(edge_colour = as.factor(col)),
    width = 0.05,
    alpha = 0.1, # Keep non-backbone edges very faint
    show.legend = FALSE
  ) +
  geom_node_point(aes(fill = general_process, size = abs(estimate)), alpha = 0.8, shape = 21, color = "#1f1f1f", stroke = 0.2) +
  shadowtext::geom_shadowtext(aes(x = x, y = y, label = label),
    size = 5.5,
    fontface = "bold",
    color = "black",
    bg.color = "white",
    bg.r = 0.1,
    check_overlap = TRUE
  ) +
  # geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf) +
  scale_size_continuous(range = c(5, 15), name = "|Delta Degree|") +
  scale_color_discrete(na.value = "gray80", name = "General Process") +
  scale_edge_color_manual(values = c("FALSE" = NA, "TRUE" = "black")) + # Optional: Hide non-backbone edges entirely
  theme_graph() +
  theme(legend.position = "right")

# p1

ggsave("results/plots/ppi_network/ppi_collapsed_dig_clusters.png", p1, width = 16, height = 12, dpi = 300)

# ==========================================
# 4. Hub Correlation Analysis for Highlighted Processes
# ==========================================
message("\nStarting Hub Correlation Analysis...")

# Load required data
deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
ssgsea_scores <- readRDS("results/ssgsea/ssgsea_ontology_scores.rds")
meta_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)
output_dir <- "results/plots/correlations/verify"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

common_samples <- intersect(colnames(deg_mat), colnames(ssgsea_scores)) %>%
  intersect(meta_sjs_ctrl$ID)

deg_mat_subset <- deg_mat[, common_samples]
ssgsea_scores_subset <- ssgsea_scores[, common_samples]

# Process each unique general_process (top_term)
unique_processes <- unique(final_summary$top_term)

all_hub_results <- list()

for (process in unique_processes) {
  # Format to GOBP
  target_pathway <- paste0("GOBP_", toupper(gsub(" ", "_", process)))

  if (!target_pathway %in% rownames(ssgsea_scores_subset)) {
    message("Pathway '", target_pathway, "' not found in ssGSEA scores. Skipping.")
    next
  }

  # Get genes in this specific cluster
  valid_nodes <- which(!is.na(V(ppi_graph)$general_process) & V(ppi_graph)$general_process == process)
  cluster_nodes <- V(ppi_graph)[valid_nodes]
  cluster_symbols <- cluster_nodes$symbol
  cluster_ensembl <- cluster_nodes$name

  # Keep only genes present in the degree matrix
  valid_idx <- cluster_ensembl %in% rownames(deg_mat_subset)
  valid_ensembl <- cluster_ensembl[valid_idx]
  valid_symbols <- cluster_symbols[valid_idx]

  if (length(valid_ensembl) == 0) {
    message("No genes from '", process, "' found in degree matrix.")
    next
  }

  # Split samples by condition
  sjs_samples <- meta_sjs_ctrl$ID[meta_sjs_ctrl$Condition == "Sjogrens"] %>% intersect(common_samples)
  ctrl_samples <- meta_sjs_ctrl$ID[meta_sjs_ctrl$Condition == "Control"] %>% intersect(common_samples)

  # Compute correlation for these genes
  path_vec <- as.numeric(ssgsea_scores_subset[target_pathway, common_samples])
  path_vec_sjs <- as.numeric(ssgsea_scores_subset[target_pathway, sjs_samples])
  path_vec_ctrl <- as.numeric(ssgsea_scores_subset[target_pathway, ctrl_samples])

  cor_results_list <- list()
  for (i in seq_along(valid_ensembl)) {
    g_id <- valid_ensembl[i]
    sym <- valid_symbols[i]

    deg_vec_sjs <- as.numeric(deg_mat_subset[g_id, sjs_samples])
    res_sjs <- tryCatch(cor.test(deg_vec_sjs, path_vec_sjs, method = "spearman", exact = FALSE), error = function(e) NULL)

    deg_vec_ctrl <- as.numeric(deg_mat_subset[g_id, ctrl_samples])
    res_ctrl <- tryCatch(cor.test(deg_vec_ctrl, path_vec_ctrl, method = "spearman", exact = FALSE), error = function(e) NULL)

    if (!is.null(res_sjs) && !is.null(res_ctrl)) {
      cor_results_list[[i]] <- data.frame(
        ensembl = g_id,
        symbol = sym,
        rho_sjs = res_sjs$estimate,
        p_val_sjs = res_sjs$p.value,
        rho_ctrl = res_ctrl$estimate,
        p_val_ctrl = res_ctrl$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(cor_results_list) == 0) next

  cor_df <- do.call(rbind, cor_results_list) %>%
    mutate(
      padj_sjs = p.adjust(p_val_sjs, method = "BH"),
      padj_ctrl = p.adjust(p_val_ctrl, method = "BH")
    )

  # Integrate topology and correlation
  hub_analysis <- cor_df %>%
    left_join(
      data.frame(
        ensembl = V(ppi_graph)$name,
        estimate = V(ppi_graph)$estimate
      ),
      by = "ensembl"
    ) %>%
    mutate(
      rho = rho_sjs,
      padj = padj_sjs,
      abs_rho = abs(rho_sjs),
      abs_estimate = abs(estimate),
      activity_score = abs_estimate * abs_rho
    )

  if (nrow(hub_analysis) == 0) next

  # Classify hubs within this cluster
  act_cutoff <- quantile(hub_analysis$activity_score, 0.90, na.rm = TRUE)
  est_cutoff <- quantile(hub_analysis$abs_estimate, 0.75, na.rm = TRUE)
  rho_cutoff <- quantile(hub_analysis$abs_rho, 0.25, na.rm = TRUE)

  hub_analysis <- hub_analysis %>%
    mutate(
      category = case_when(
        activity_score >= act_cutoff & padj < 0.05 ~ "Active Hub",
        abs_estimate >= est_cutoff & abs_rho < rho_cutoff ~ "Silent Hub",
        abs_estimate < est_cutoff ~ "Non-Hub",
        TRUE ~ "Intermediate"
      ),
      process = process,
      target_pathway = target_pathway
    )

  all_hub_results[[process]] <- hub_analysis

  selected_hubs <- hub_analysis %>%
    filter(category %in% c("Active Hub", "Silent Hub", "Intermediate"))

  if (nrow(selected_hubs) == 0) {
    message("No Selected Hubs found for ", process)
    next
  }

  message("Found ", nrow(selected_hubs), " hubs to plot for ", process)

  # Plot all selected hubs
  for (j in 1:nrow(selected_hubs)) {
    target_gene <- selected_hubs$symbol[j]
    g_id <- selected_hubs$ensembl[j]
    hub_cat <- selected_hubs$category[j]

    deg_vec <- as.numeric(deg_mat_subset[g_id, common_samples])

    cor_data <- data.frame(
      Degree = deg_vec,
      PathwayScore = path_vec,
      ID = common_samples
    ) %>%
      left_join(meta_sjs_ctrl, by = "ID")

    rho_sjs_val <- selected_hubs$rho_sjs[j]
    pval_sjs_val <- selected_hubs$p_val_sjs[j]
    rho_ctrl_val <- selected_hubs$rho_ctrl[j]
    pval_ctrl_val <- selected_hubs$p_val_ctrl[j]

    subtitle_text <- sprintf(
      "Category: %s\nSJS Rho: %.3f (p=%s) | Ctrl Rho: %.3f (p=%s)",
      hub_cat,
      rho_sjs_val, format.pval(pval_sjs_val, digits = 2),
      rho_ctrl_val, format.pval(pval_ctrl_val, digits = 2)
    )

    p <- ggplot(cor_data, aes(x = Degree, y = PathwayScore, color = Condition, fill = Condition)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, linetype = "dashed") +
      theme_minimal() +
      labs(
        title = paste(target_gene, "vs", target_pathway),
        subtitle = subtitle_text,
        x = paste("Gene Degree (", target_gene, ")"),
        y = "Pathway Activity (ssGSEA)",
        color = "Condition", fill = "Condition"
      )

    # Format the filename to include the category
    clean_cat <- gsub(" ", "_", hub_cat)
    plot_filename <- file.path(output_dir, paste0(clean_cat, "_", target_gene, "_", target_pathway, ".svg"))
    ggsave(plot_filename, p, width = 10, height = 7, bg = "white")
  }
}

if (length(all_hub_results) > 0) {
  final_hub_df <- do.call(rbind, all_hub_results)
  sig_correlations <- final_hub_df %>% filter(padj_sjs < 0.05 | padj_ctrl < 0.05)
  write.csv(sig_correlations, file.path("results/significant_hub_correlations.csv"), row.names = FALSE)
}

message("Hub correlation plotting complete.")
