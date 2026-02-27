library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

my_digs <- read.csv("results/DIGs_significant_padj05.csv")$gene
my_digs_clean <- gsub("\\..*$", "", my_digs)

ego <- enrichGO(gene          = my_digs_clean,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.csv(as.data.frame(ego), "results/GO_Enrichment_Results.csv", row.names = FALSE)

p1d <- barplot(ego, showCategory=10, font.size = 10) +
  theme_bw() +
  labs(title = "GO Biological Processes of DIGs",
       x = "Gene Count")

p1d

metabolic_terms <- ego@result %>%
  filter(grepl("metabolic", Description, ignore.case = TRUE)) %>%
  arrange(pvalue)

write.csv(metabolic_terms, "results/GO_Metabolic_Processes_Only.csv", row.names = FALSE)

p1e <- ggplot(head(metabolic_terms, 10), aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(title = "Top Metabolic Processes of DIGs", x = "Process", y = "-log10(p-value)")

p1e

entrez_ids <- bitr(my_digs_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa', 
                    pvalueCutoff = 0.05)

write.csv(as.data.frame(ekegg), "results/KEGG_Pathway_Results.csv", row.names = FALSE)

p1f <- cnetplot(setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID'), 
                circular = TRUE, 
                colorEdge = TRUE, 
                showCategory = 5) +
  labs(title = "KEGG Pathway-Gene Network")

p1f

ggsave("results/plots/enrich_plot.png", p1d, width = 8, height = 6)
ggsave("results/plots/enrich_metabolic_plot.png", p1e, width = 8, height = 6)
