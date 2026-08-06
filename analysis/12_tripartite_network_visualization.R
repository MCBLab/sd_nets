#!/usr/bin/env Rscript
# Script to generate a tripartite network visualization: Hub -> First-Order Neighbor -> Drug

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(igraph)
})

# 1. Load Data
df_hubs <- read.csv("results/significant_hub_correlations.csv", stringsAsFactors = FALSE)
df_drugs <- read.csv("results/repurposing/sjogrens_vs_control_drug_prioritization.csv", stringsAsFactors = FALSE)
df_db_approval <- read.csv("data/drugbank_approval_status.csv", stringsAsFactors = FALSE)
df_drugs <- df_drugs %>% left_join(df_db_approval, by = c("drug_id" = "drugbank_id"))
ppi_g <- readRDS("results/ppi_network/ppi_collapsed_dig_graph.rds")

# Common function to generate network
create_hub_network <- function(hubs_to_plot, output_suffix, only_approved = FALSE) {
  
  valid_hubs_in_graph <- intersect(hubs_to_plot$ensembl, V(ppi_g)$name)
  
  if(length(valid_hubs_in_graph) == 0) {
    message(sprintf("No valid hub genes found in the PPI network for %s.", output_suffix))
    return()
  }
  
  message(sprintf("Found %d valid Hub Genes in the network for %s.", length(valid_hubs_in_graph), output_suffix))
  
  # 3. Extract First-Order Neighbors
  neighbors_list <- list()
  for(hub_id in valid_hubs_in_graph) {
    hub_symbol <- V(ppi_g)$symbol[V(ppi_g)$name == hub_id]
    hub_process <- hubs_to_plot$process[hubs_to_plot$ensembl == hub_id][1]
    
    nei <- neighbors(ppi_g, hub_id, mode="all")
    nei_symbols <- V(ppi_g)$symbol[nei]
    nei_symbols <- nei_symbols[!is.na(nei_symbols)]
    
    if(length(nei_symbols) > 0) {
      neighbors_list[[length(neighbors_list) + 1]] <- data.frame(
        hub = hub_symbol,
        neighbor = nei_symbols,
        process = hub_process,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df_neighbors <- do.call(rbind, neighbors_list) %>%
    distinct()
  
  # 4. Map Drug Associations
  # Identify valid drugs
  valid_drug_data <- df_drugs %>%
    filter(!is.na(mechanism), mechanism != "", mechanism != "NA", source == "DrugBank")
  
  if (only_approved) {
    valid_drug_data <- valid_drug_data %>% filter(is_approved == "True")
  }
  
  # Filter neighbors to only those targeted by a valid drug
  neighbors_with_drugs <- intersect(df_neighbors$neighbor, valid_drug_data$gene_name)
  
  if(length(neighbors_with_drugs) == 0) {
    message(sprintf("None of the first-order neighbors have mapped drugs for %s.", output_suffix))
    return()
  }
  
  message(sprintf("Filtered down to %d Neighbors with clinical drug mappings for %s.", length(neighbors_with_drugs), output_suffix))
  
  # Filter the edge lists
  edges_hub_neighbor <- df_neighbors %>%
    filter(neighbor %in% neighbors_with_drugs) %>%
    select(from = hub, to = neighbor) %>%
    mutate(edge_type = "Hub-Neighbor")
  
  edges_neighbor_drug <- valid_drug_data %>%
    filter(gene_name %in% neighbors_with_drugs) %>%
    select(from = gene_name, to = drug_name) %>%
    distinct(from, to) %>%
    mutate(edge_type = "Neighbor-Drug")
  
  edges <- bind_rows(edges_hub_neighbor, edges_neighbor_drug)
  
  # 5. Create Nodes and map attributes
  # Determine module mapping for hubs
  hub_process_map <- df_neighbors %>%
    select(hub, process) %>%
    distinct(hub, .keep_all = TRUE)
  
  nodes_hub <- data.frame(name = unique(edges_hub_neighbor$from), type = "Hub Gene") %>%
    left_join(hub_process_map, by = c("name" = "hub")) %>%
    mutate(module = process) %>%
    select(-process)
  
  nodes_neighbor <- data.frame(name = unique(c(edges_hub_neighbor$to, edges_neighbor_drug$from)), type = "Neighbor Gene") %>%
    mutate(module = "Neutral") # Requested to be neutral
  
  nodes_drug <- data.frame(name = unique(edges_neighbor_drug$to), type = "Drug") %>%
    mutate(module = "Neutral")
  
  nodes <- bind_rows(nodes_hub, nodes_neighbor, nodes_drug)
  
  # 6. Build Graph and Layout
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  nodes_df <- as_tibble(g) %>%
    mutate(id = row_number())
  
  space_y <- function(n) {
    if (n == 1) return(5)
    seq(0, 10, length.out = n)
  }
  
  nodes_df <- nodes_df %>%
    group_by(type) %>%
    mutate(
      x = case_when(
        type == "Hub Gene" ~ 1,
        type == "Neighbor Gene" ~ 2,
        type == "Drug" ~ 3
      ),
      y = space_y(n())
    ) %>%
    ungroup() %>%
    arrange(id)
  
  layout <- create_layout(g, layout = "manual", x = nodes_df$x, y = nodes_df$y)
  
  # 7. Aesthetics and Formatting
  # Custom colors for the specific hub processes
  pathway_colors <- c(
    "defense response to virus" = "#8B9A46", # Olive
    "B cell activation" = "#E27D60",         # Coral
    "Neutral" = "#F5F5F5"                    # Light Gray
  )
  # Dynamic fallback if other pathways exist
  for(pw in unique(nodes_hub$module)) {
    if(!(pw %in% names(pathway_colors))) pathway_colors[pw] <- scales::hue_pal()(1)
  }
  
  # Plot
  p <- ggraph(layout) +
    geom_edge_link(aes(color = edge_type), edge_width = 0.8, edge_alpha = 0.6, show.legend = FALSE) +
    scale_edge_color_manual(values = c("Hub-Neighbor" = "#4A8B9E", "Neighbor-Drug" = "gray75")) +
    
    geom_node_point(aes(fill = module, shape = type), size = 5.5, color = "gray20", stroke = 0.5) +
    scale_shape_manual(values = c("Hub Gene" = 22, "Neighbor Gene" = 21, "Drug" = 24)) + 
    scale_fill_manual(values = pathway_colors) +
    
    geom_node_text(aes(label = name, filter = type == "Hub Gene"), hjust = 1, nudge_x = -0.05, size = 3.5, fontface = "bold") +
    geom_node_text(aes(label = name, filter = type == "Neighbor Gene"), vjust = -1.5, size = 3) +
    geom_node_text(aes(label = name, filter = type == "Drug"), hjust = 0, nudge_x = 0.05, size = 3) +
    
    theme_graph(base_family = "sans") +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(20, 80, 20, 80)
    ) +
    guides(
      shape = guide_legend(override.aes = list(fill = "gray50")),
      fill = guide_legend(override.aes = list(shape = 22))
    )
  
  # Calculate dynamic height based on the maximum number of nodes in a column
  max_nodes <- max(table(nodes_df$type))
  plot_height <- max(10, min(max_nodes * 0.3, 30)) # At least 10, max 30 to avoid gigantic files
  
  pdf_out <- sprintf("results/repurposing/tripartite_network_hub_neighbor_drug_%s.pdf", output_suffix)
  png_out <- sprintf("results/repurposing/tripartite_network_hub_neighbor_drug_%s.png", output_suffix)
  
  ggsave(pdf_out, p, width = 12, height = plot_height, device = cairo_pdf, limitsize = FALSE)
  ggsave(png_out, p, width = 12, height = plot_height, dpi = 300, limitsize = FALSE)
  
  message(sprintf("Successfully saved plot to %s with dynamic height %.1f", pdf_out, plot_height))
}

# Run 1: Original target processes
target_processes <- c("defense response to virus", "B cell activation")
hubs_all <- df_hubs %>%
  filter(process %in% target_processes, category %in% c("Active Hub", "Silent Hub")) %>%
  distinct(ensembl, symbol, process)

create_hub_network(hubs_all, "all")

# Run 2: Specific genes
specific_genes <- c("OASL", "LGALS3BP", "ISG15", "UBE2L6")

# Even if they are not explicitly "Active Hub" in the requested pathways, we pull them from df_hubs as long as they are present
hubs_specific <- df_hubs %>%
  filter(symbol %in% specific_genes) %>%
  distinct(ensembl, symbol, process)

create_hub_network(hubs_specific, "specific")
create_hub_network(hubs_specific, "specific_approved", only_approved = TRUE)
