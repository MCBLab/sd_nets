library(GSVA)
library(GSEABase)
library(dplyr)
library(tibble)
library(ggplot2)
library(vroom)


normalized_counts <- readRDS("results/vst_normalized_counts_symbols.rds")

gene_sets <- getGmt("data/h.all.v2026.1.Hs.symbols.gmt")

params <- ssgseaParam(
  exprData = as.matrix(normalized_counts),
  geneSets = gene_sets,
  minSize = 10,
  maxSize = 500
)

ssgsea_res <- gsva(params)

dir.create("results/ssgsea", showWarnings = FALSE)
saveRDS(ssgsea_res, "results/ssgsea/ssgsea_hallmark_scores.rds")
write.csv(as.data.frame(ssgsea_res), "results/ssgsea/ssgsea_scores_matrix.csv")

# Exemplo: SIGLEC1 (Hub do Cluster 1/IFN) vs Pathway Interferon Alpha Response
hub_gene <- "ENSG00000088827" # SIGLEC1
target_pathway <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"

deg_mat <- readRDS("results/lioness_gene_degree_matrix.rds")
metadata_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv")

# Extrair grau do hub da sua deg_mat
hub_degrees <- as.numeric(deg_mat[hub_gene, ])
pathway_scores <- as.numeric(ssgsea_res[target_pathway, ])

cor_data <- data.frame(
  ID = colnames(deg_mat),
  Degree = hub_degrees,
  PathwayScore = pathway_scores
) %>% left_join(metadata_sjs_ctrl, by = "ID")

# Cálculo da correlação de Spearman (não paramétrica)
cor_test <- cor.test(cor_data$Degree, cor_data$PathwayScore, method = "spearman")

p_cor <- ggplot(cor_data, aes(x = Degree, y = PathwayScore, color = Condition)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste("Conectividade de", hub_gene, "vs", target_pathway),
    subtitle = paste("Spearman Rho:", round(cor_test$estimate, 3), "| p-val:", format.pval(cor_test$p.value)),
    x = "Grau (Conectividade na Rede)",
    y = "Score ssGSEA (Atividade da Via)"
  )

p_cor

ggsave("results/plots/hub_pathway_correlation.png", p_cor, width = 7, height = 6)