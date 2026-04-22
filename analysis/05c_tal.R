library(dplyr)
library(rrvgo)
library(ggraph)
library(GOSemSim)
library(tidyverse)
library(tidygraph)
library(TreeAndLeaf)
library(RColorBrewer)
library(org.Hs.eg.db)

# Output directory for plots
dir.create("results/plots/enrichment", showWarnings = FALSE, recursive = TRUE)


input_file <- "results/enrichment/GO_sjs_vs_ctrl.csv"
output_csv <- "results/enrichment/GO_semantic_reduced_sjs_vs_ctrl.csv"
output_plot <- "results/plots/enrichment/GO_treemap_sjs_vs_ctrl.png"

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

enrich_df <- read_csv(input_file)

if (nrow(enrich_df) == 0) {
  stop("Input enrichment dataframe is empty.")
}

unique_go_ids <- unique(enrich_df$ID)


scores_df <- enrich_df %>%
  group_by(ID) %>%
  summarise(min_padj = min(p.adjust, na.rm = TRUE)) %>%
  mutate(score = -log10(min_padj))

scores <- setNames(scores_df$score, scores_df$ID)

simMatrix <- calculateSimMatrix(
  unique_go_ids,
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel" 
)

reducedTerms <- reduceSimMatrix(
  simMatrix,
  scores,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)

if ("Cluster" %in% colnames(enrich_df)) {
  enrich_reduced_df <- enrich_df %>%
    left_join(reducedTerms, by = c("ID" = "go")) %>%
    dplyr::select(
      Cluster, 
      ID, 
      Description, 
      p.adjust, 
      parent,       
      parentTerm,   
      everything()  
    )
} else {
  enrich_reduced_df <- enrich_df %>%
    left_join(reducedTerms, by = c("ID" = "go")) %>%
    dplyr::select(
      ID, 
      Description, 
      p.adjust, 
      parent,       
      parentTerm,   
      everything()  
    )
}

write_csv(enrich_reduced_df, output_csv)

png(output_plot, width = 1000, height = 800)
treemapPlot(reducedTerms)
dev.off()

semData <- godata('org.Hs.eg.db', ont="BP")
mSim <- mgoSim(
  unique_go_ids,
  unique_go_ids,
  semData = semData,
  measure = "Wang",
  combine = NULL
)

enrich_prepared <- enrich_df %>%
  mutate(
    path_length = as.integer(sapply(strsplit(GeneRatio, "/"), "[", 2)),
    ratio = Count / path_length
  )

hc <- hclust(dist(mSim), "average")
tal <- treeAndLeaf(hc)

max_gg_size <- 8
min_gg_size <- 2

tal_tidy <- as_tbl_graph(tal) %>%
  activate(nodes) %>%
  left_join(enrich_prepared, by = c("name" = "ID")) %>%
  mutate(
    color_bin = ntile(ratio, 5),
    nodeColor = case_when(
      is.na(ratio) ~ "#191919", 
      TRUE ~ colorRampPalette(brewer.pal(9, "OrRd"))(5)[color_bin]
    ),
    size_bin = ntile(Count, 7), 
    nodeSize = case_when(
      isLeaf & !is.na(Count) ~ min_gg_size + (size_bin * 0.8),
      isLeaf & is.na(Count)  ~ 1,
      TRUE ~ 0
    ),
    
    nodeLabel = ifelse(isLeaf, Description, NA)
  )

layout <- easylayout(tal)

p <- ggraph(tal_tidy, layout = 'manual', x = layout[, 1], y = layout[, 2]) + 
  geom_edge_diagonal(color = "#191919", alpha = 0.5) + 
  geom_node_point(aes(color = nodeColor, size = nodeSize)) +
  geom_node_text(aes(label = nodeLabel), 
                 repel = TRUE, size = 3, max.overlaps = 50) +
  scale_color_identity() + 
  scale_size_identity() +
  theme_void() +
  labs(title = "SJS vs Control GO Enrichment",
       subtitle = "Node Size: Ranked Count | Node Color: Ratio Quantiles")

p

ggsave(
  filename = "results/plots/enrichment/GO_TreeAndLeaf_sjs_vs_ctrl.pdf",
  plot = p,
  device = "pdf",
  width = 15,
  height = 15,
  units = "in",
  limitsize = FALSE
)
