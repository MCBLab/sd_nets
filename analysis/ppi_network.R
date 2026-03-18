library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(vroom)

digs_df <- read.csv("results/digs/sjs_vs_ctrl_significant_padj05.csv")

digs_df <- digs_df %>% rename(diff_deg = estimate)

digs_df <- digs_df %>% filter(!is.na(symbol))

aliases <- vroom("data/9606.protein.aliases.v12.0.txt.gz", 
                 comment = "#", col_names = c("string_id", "alias", "source"),
                 show_col_types = FALSE)
links <- vroom("data/9606.protein.links.detailed.v12.0.txt.gz", 
               delim = " ", show_col_types = FALSE)

# Mapping gene symbols to STRING IDs
symbol_to_string_map <- aliases %>%
  filter(alias %in% digs_df$symbol) %>%
  mutate(string_id = gsub("9606\\.", "", string_id)) %>%
  select(symbol = alias, STRING_id = string_id) %>%
  distinct()

mapped_genes <- digs_df %>%
  inner_join(symbol_to_string_map, by = "symbol")

score_threshold <- 400
dig_string_ids <- unique(mapped_genes$STRING_id)

interactions <- links %>%
  filter(combined_score >= score_threshold) %>%
  mutate(
    protein1 = gsub("9606\\.", "", protein1),
    protein2 = gsub("9606\\.", "", protein2)
  ) %>%
  filter(protein1 %in% dig_string_ids & protein2 %in% dig_string_ids) %>%
  select(from = protein1, to = protein2, combined_score)

edges <- interactions %>% distinct()
node_ids_in_network <- unique(c(edges$from, edges$to))

nodes <- mapped_genes %>%
  filter(STRING_id %in% node_ids_in_network) %>%
  group_by(STRING_id) %>%
  summarise(
    symbol = first(symbol),
    diff_deg = first(diff_deg),
    p_val = first(p_val),
    padj = first(padj),
    .groups = "drop"
  ) %>%
  rename(name = STRING_id)

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree(),
         betweenness = centrality_betweenness())

write.csv(as.data.frame(edges), "results/PPI_Edges_DIGs.csv", row.names = FALSE)
write.csv(as.data.frame(nodes), "results/PPI_Nodes_DIGs.csv", row.names = FALSE)

final_layout <- easylayout::easylayout(graph)

p_ppi <- ggraph(graph, layout = final_layout) +
  geom_edge_link(aes(alpha = combined_score), color = "grey", show.legend = FALSE) +
  geom_node_point(aes(color = diff_deg), size = 5, alpha = 0.8) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Δ Degree") +
  theme_graph() +
  labs(title = "Rede de Interação PPI (DIGs)",
       subtitle = "Layout: easylayout | Cores: Mudança na conectividade")

ggsave("results/plots/ppi_network_digs.png", p_ppi, width = 10, height = 8)
