#### StringDB export

library(easylayout)
library(igraph)
library(dplyr)
library(ggrepel)
library(tidyr)

setwd("/Documents and Settings/diego.coelho/Documents/RVU501/")

string <- read.table("sigs.string_network_coordinates.tsv", fill = T)
# sigs.string_network_coordinates.tsv is a file exported from StringDB with the following columns:
#node1	node2	node1_string_id	node2_string_id	neighborhood_on_chromosome	gene_fusion	phylogenetic_cooccurrence	homology	coexpression	experimentally_determined_interaction	database_annotated	automated_textmining	combined_score
#Abhd5	Pparg	10090.ENSMUSP00000122274	10090.ENSMUSP00000000450	0	0	0	0	0.083	0	0	0.497	0.519
#Abhd5	Aloxe3	10090.ENSMUSP00000122274	10090.ENSMUSP00000021268	0	0	0	0	0	0	0	0.589	0.589
#Abhd5	Fabp4	10090.ENSMUSP00000122274	10090.ENSMUSP00000029041	0	0	0	0	0.066	0	0	0.637	0.646


FCs <- readxl::read_excel("XXX_XXX_Lung_RNAseq_meta.xlsx")
sigs <- FCs %>% filter(padj_Vehicle_vs_XXX_12d < 0.05 | padj_Vehicle_vs_XXX_15d < 0.05 |
                         padj_Vehicle_vs_XXX_17_18d < 0.05) %>% pull(log2FoldChange_Vehicle_vs_XXX_17_18d)

ppi_graph <- graph_from_data_frame(d = string[,1:2], directed = F)
connected_nodes <- V(ppi_graph)[degree(ppi_graph) > 0]

### Day 12d

V(ppi_graph)$expression <- FCs$log2FoldChange_Vehicle_vs_XXX_12d[match(V(ppi_graph)$name, FCs$gene_name)]

# Convert igraph to data frame
layout <- layout_with_fr(ppi_graph) %>% as.data.frame()
colnames(layout) <- c("x", "y")
layout$protein <- V(ppi_graph)$name
layout$degree <- degree(ppi_graph)
layout$expression <- V(ppi_graph)$expression
layout$highlight <- FCs$padj_Vehicle_vs_XXX_12d[match(V(ppi_graph)$name, FCs$gene_name)] < 0.05

edges <- as.data.frame(get.edgelist(ppi_graph))
colnames(edges) <- c("source", "target")

# Merge coordinates with edges
edges <- edges %>%
  left_join(layout, by = c("source" = "protein")) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("target" = "protein")) %>%
  dplyr::rename(x2 = x, y2 = y)

# Plot using ggplot2
png(filename = "figures/PPI.DGE.VehicleVsXXX.12d.FC.png", res = 300, height = 1500, width = 2000)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Vehicle vs. XXX - Day 12",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

### Day 15d

V(ppi_graph)$expression <- FCs$log2FoldChange_Vehicle_vs_XXX_15d[match(V(ppi_graph)$name, FCs$gene_name)]

# Convert igraph to data frame
layout$expression <- V(ppi_graph)$expression
layout$highlight <- FCs$padj_Vehicle_vs_XXX_15d[match(V(ppi_graph)$name, FCs$gene_name)] < 0.05

edges <- as.data.frame(get.edgelist(ppi_graph))
colnames(edges) <- c("source", "target")

# Merge coordinates with edges
edges <- edges %>%
  left_join(layout, by = c("source" = "protein")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("target" = "protein")) %>%
  rename(x2 = x, y2 = y)

# Plot using ggplot2
png(filename = "figures/PPI.DGE.VehicleVsXXX.15d.FC.png", res = 300, height = 1500, width = 2000)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Vehicle vs. XXX - Day 15",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

### Day 17d

V(ppi_graph)$expression <- FCs$log2FoldChange_Vehicle_vs_XXX_17_18d[match(V(ppi_graph)$name, FCs$gene_name)]

# Convert igraph to data frame
layout$expression <- V(ppi_graph)$expression
layout$highlight <- FCs$padj_Vehicle_vs_XXX_17_18d[match(V(ppi_graph)$name, FCs$gene_name)] < 0.05

edges <- as.data.frame(get.edgelist(ppi_graph))
colnames(edges) <- c("source", "target")

# Merge coordinates with edges
edges <- edges %>%
  left_join(layout, by = c("source" = "protein")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("target" = "protein")) %>%
  rename(x2 = x, y2 = y)

# Plot using ggplot2
png(filename = "figures/PPI.DGE.VehicleVsXXX.17.18d.FC.png", res = 300, height = 1500, width = 2000)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Vehicle vs. XXX - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

#############

string <- read.table("all.responders.string_network_coordinates.tsv", fill = T)

FCs <- readxl::read_excel("58321_OIV2024_Lung_RNAseq_meta.xlsx")
sigs <- FCs %>% filter(padj_non_vs_responder_17_18d < 0.05 | padj_vehicle_vs_non_responder_17_18d < 0.05 |
                         padj_vehicle_vs_responder_17_18d < 0.05) %>% pull(log2FoldChange_non_vs_responder_17_18d)

ppi_graph <- graph_from_data_frame(d = string[,1:2], directed = F)
# Identify nodes with at least one connection
connected_nodes <- V(ppi_graph)[degree(ppi_graph) > 0]
# Subset the graph to include only connected nodes
ppi_graph <- induced_subgraph(ppi_graph, connected_nodes)

# Apply Louvain clustering
clusters <- cluster_louvain(ppi_graph, resolution = 0.3)

# Add cluster membership to the node attributes
V(ppi_graph)$cluster <- membership(clusters)


### Non vs. Responder

V(ppi_graph)$expression <- FCs$log2FoldChange_non_vs_responder_17_18d[match(V(ppi_graph)$name, FCs$gene_name)]

# Convert igraph to data frame
layout <- layout_with_fr(ppi_graph) %>% as.data.frame()
colnames(layout) <- c("x", "y")
layout$protein <- V(ppi_graph)$name
layout$expression <- V(ppi_graph)$expression
layout$highlight <- FCs$padj_non_vs_responder_17_18d[match(V(ppi_graph)$name, FCs$gene_name)] < 0.05
layout$cluster <- V(ppi_graph)$cluster
keep <- names(table(layout$cluster)[!is.na(ifelse(table(layout$cluster) < 10,
                                                  NA,
                                                  table(layout$cluster)))])
layout$cluster <- ifelse(layout$cluster %in% keep, layout$cluster, NA)

clusters[[1]] %>% write.table(file = "test_cluster1.txt", quote = F, row.names = F, col.names = F)
clusters[[2]] %>% write.table(file = "test_cluster2.txt", quote = F, row.names = F, col.names = F)
clusters[[7]] %>% write.table(file = "test_cluster7.txt", quote = F, row.names = F, col.names = F)

edges <- as.data.frame(get.edgelist(ppi_graph))
colnames(edges) <- c("source", "target")

# Merge coordinates with edges
edges <- edges %>%
  left_join(layout, by = c("source" = "protein")) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("target" = "protein")) %>%
  dplyr::rename(x2 = x, y2 = y)

# Plot using ggplot2
png(filename = "figures/PPI.DGE.ResponderVsNonResponder.17.18d.FC.png", res = 300, height = 2400, width = 3200)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Responder vs. Non-responder - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

# Plot using ggplot2
png(filename = "figures/PPI.DGE.Responders_clusters.17.18d.FC.png", res = 300, height = 2400, width = 3200)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = factor(cluster), size = 20)) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_manual(values = rainbow(length(unique(layout$cluster)))) +
  # scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with clusters - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

### Sub cluster1
ppi_graph_c1 <- induced_subgraph(ppi_graph, V(ppi_graph)[ifelse(ifelse(is.na(layout$cluster), 0, layout$cluster) == 1, T,F)])
layout_sub <- easylayout(ppi_graph_c1)
colnames(layout_sub) <- c("x", "y")
layout_sub <- cbind(layout_sub, layout[ifelse(ifelse(is.na(layout$cluster), 0, layout$cluster) == 1, T,F),][,-c(1,2)])

edges_sub <- as.data.frame(get.edgelist(ppi_graph_c1))
colnames(edges_sub) <- c("source", "target")

# Merge coordinates with edges
edges_sub <- edges_sub %>%
  left_join(layout_sub, by = c("source" = "protein")) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  left_join(layout_sub, by = c("target" = "protein")) %>%
  dplyr::rename(x2 = x, y2 = y)


png(filename = "figures/PPI.DGE.ResponderVsNonResponder_clus1.17.18d.FC.png", res = 300, height = 2000, width = 2500)
ggplot() +
  geom_segment(data = edges_sub, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout_sub, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout_sub %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout_sub, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Responder vs. Non-responder (Cluster1) - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",dim(layout_sub)[1])) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

### Vehicle vs. Non-responder

V(ppi_graph)$expression <- FCs$log2FoldChange_vehicle_vs_non_responder_17_18d[match(V(ppi_graph)$name, FCs$gene_name)]

# Convert igraph to data frame
layout$expression <- V(ppi_graph)$expression
layout$highlight <- FCs$padj_vehicle_vs_non_responder_17_18d[match(V(ppi_graph)$name, FCs$gene_name)] < 0.05

edges <- as.data.frame(get.edgelist(ppi_graph))
colnames(edges) <- c("source", "target")

# Merge coordinates with edges
edges <- edges %>%
  left_join(layout, by = c("source" = "protein")) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("target" = "protein")) %>%
  dplyr::rename(x2 = x, y2 = y)

# Plot using ggplot2
png(filename = "figures/PPI.DGE.VehicleVsNonResponder.17.18d.FC.png", res = 300, height = 2400, width = 3200)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Non-Responder vs. Vehicle - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

### Sub cluster1
layout_sub$expression <- layout[ifelse(ifelse(is.na(layout$cluster), 0, layout$cluster) == 1, T,F),]$expression
layout_sub$highlight <- layout[ifelse(ifelse(is.na(layout$cluster), 0, layout$cluster) == 1, T,F),]$highlight

png(filename = "figures/PPI.DGE.VehicleVsNonResponder_clus1.17.18d.FC.png", res = 300, height = 2000, width = 2500)
ggplot() +
  geom_segment(data = edges_sub, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout_sub, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout_sub %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout_sub, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Non-Responder vs. Vehicle (Cluster 1) - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",dim(layout_sub)[1])) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()


### Vehicle vs. Responder

V(ppi_graph)$expression <- FCs$log2FoldChange_vehicle_vs_responder_17_18d[match(V(ppi_graph)$name, FCs$gene_name)]

# Convert igraph to data frame
layout$expression <- V(ppi_graph)$expression
layout$highlight <- FCs$padj_vehicle_vs_responder_17_18d[match(V(ppi_graph)$name, FCs$gene_name)] < 0.05

edges <- as.data.frame(get.edgelist(ppi_graph))
colnames(edges) <- c("source", "target")

# Merge coordinates with edges
edges <- edges %>%
  left_join(layout, by = c("source" = "protein")) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("target" = "protein")) %>%
  dplyr::rename(x2 = x, y2 = y)

# Plot using ggplot2
png(filename = "figures/PPI.DGE.VehicleVsResponder.17.18d.FC.png", res = 300, height = 2400, width = 3200)
ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Responder vs. Vehicle - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",length(connected_nodes))) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()

### Sub cluster1
layout_sub$expression <- layout[ifelse(ifelse(is.na(layout$cluster), 0, layout$cluster) == 1, T,F),]$expression
layout_sub$highlight <- layout[ifelse(ifelse(is.na(layout$cluster), 0, layout$cluster) == 1, T,F),]$highlight

png(filename = "figures/PPI.DGE.VehicleVsResponder_clus1.17.18d.FC.png", res = 300, height = 2000, width = 2500)
ggplot() +
  geom_segment(data = edges_sub, aes(x = x1, y = y1, xend = x2, yend = y2), color = "grey") +
  geom_point(data = layout_sub, aes(x = x, y = y, color = expression, size = 20)) +
  geom_point(data = layout_sub %>% filter(highlight), aes(x = x, y = y), 
             shape = 21, size = 4, color = "#351c75", alpha = 0.7, fill = NA, stroke = 1.5) +
  geom_text_repel(data = layout_sub, aes(x = x, y = y, label = protein), size = 3, max.overlaps = 10) +  # Avoid overlapping
  scale_color_gradientn(limits = c(-5,5), colors = c("dodgerblue2", "#ededed", "firebrick2"), na.value = "#ededed") +
  guides(size = "none") +  # Remove the size legend
  labs(title = "    PPI Network with DEGs from Responder vs. Vehicle (Cluster 1) - Day 17/18",
       color = "Fold-Change    ",
       caption = paste0("* Circles in purple are DEGs. Total genes = ",dim(layout_sub)[1])) +
  theme_void() + theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
dev.off()
