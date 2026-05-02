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
                  size = 3.5,
                  fontface = "bold",
                  color = "black",
                  bg.color = "white",
                  bg.r = 0.1,
                  check_overlap = TRUE) +
  # geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf) +
  scale_size_continuous(range = c(3, 10), name = "|Delta Degree|") +
  scale_color_discrete(na.value = "gray80", name = "General Process") +
  scale_edge_color_manual(values = c("FALSE" = NA, "TRUE" = "black")) + # Optional: Hide non-backbone edges entirely
  theme_graph() +
  theme(legend.position = "right")

p1

ggsave("results/plots/ppi_network/ppi_collapsed_dig_clusters.svg", p1, width = 16, height = 12, dpi = 300)
