library(dplyr)
library(tibble)
library(vroom)
library(limma)
library(readr)

# --- 1. Carregamento de Dados ---
# Substitua pelo caminho real do seu CSV do UCell
ucell_raw <- vroom("UCell/ucell_scores.csv") 

# Converter para matriz: Amostras nas colunas, Caminhos nas linhas
# O snippet mostra que a primeira coluna é o ID da amostra
ucell_mat <- ucell_raw %>%
  column_to_rownames(var = colnames(ucell_raw)[1]) %>%
  t()

# Limpar nomes das linhas (remover o sufixo _UCell para facilitar a leitura)
rownames(ucell_mat) <- gsub("_UCell$", "", rownames(ucell_mat))

# Metadados
metadata_sjs_ctrl <- vroom("data/precisesads/metadata_sjs_ctrl.csv", show_col_types = FALSE)
clusters <- vroom("data/precisesads/sjogren_clusters.csv", show_col_types = FALSE)

dir.create("results/ucell_limma", showWarnings = FALSE)

# --- 2. Preparação das Comparações ---

# Sincronizar amostras
common_samples <- intersect(colnames(ucell_mat), metadata_sjs_ctrl$ID)
ucell_mat <- ucell_mat[, common_samples]

# --- 3. Análise por Cluster (One-vs-All) ---
metadata_clusters <- metadata_sjs_ctrl %>% 
  inner_join(clusters, by = "ID") %>%
  filter(ID %in% colnames(ucell_mat))

if (nrow(metadata_clusters) > 0) {
  ucell_clusters <- ucell_mat[, metadata_clusters$ID]
  group <- factor(paste0("Cluster_", metadata_clusters$PREDICTION))
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(ucell_clusters, design)
  
  for (cl in levels(group)) {
    others <- levels(group)[levels(group) != cl]
    contrast_formula <- paste0(cl, " - (", paste(others, collapse = " + "), ") / ", length(others))
    contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit_contrast <- contrasts.fit(fit, contrast_matrix)
    fit_contrast <- eBayes(fit_contrast)
    
    res <- topTable(fit_contrast, coef = 1, number = Inf) %>%
      rownames_to_column("pathway") %>%
      rename(diff_score = logFC, p_val = P.Value, padj = adj.P.Val, t_stat = t) %>%
      arrange(p_val)
    
    clean_name <- gsub("Cluster_", "cluster_", cl)
    write_csv(res, paste0("results/ucell_limma/limma_ucell_", clean_name, "_vs_others.csv"))
  }
}

# --- 4. Sjögren vs Control ---
metadata_general <- metadata_sjs_ctrl %>%
  filter(ID %in% colnames(ucell_mat))

if (nrow(metadata_general) > 0) {
  ucell_general <- ucell_mat[, metadata_general$ID]
  # Ajustar nome da condição para bater com o CSV (ex: Sjögren/Sjogrens)
  condition <- factor(metadata_general$Condition)
  
  # Criar design e garantir nomes válidos
  design_sjs <- model.matrix(~ 0 + condition)
  colnames(design_sjs) <- make.names(levels(condition))
  
  fit_sjs <- lmFit(ucell_general, design_sjs)
  
  # Identificar o nome das colunas dinamicamente (ajuste se necessário)
  target_col <- grep("Sjogren|Sjs", colnames(design_sjs), value = TRUE, ignore.case = TRUE)
  ref_col <- grep("Control", colnames(design_sjs), value = TRUE, ignore.case = TRUE)
  
  if (length(target_col) > 0 && length(ref_col) > 0) {
    contrast_sjs <- makeContrasts(paste0(target_col, " - ", ref_col), levels = design_sjs)
    fit_sjs <- contrasts.fit(fit_sjs, contrast_sjs)
    fit_sjs <- eBayes(fit_sjs)
    
    res_sjs <- topTable(fit_sjs, coef = 1, number = Inf) %>%
      rownames_to_column("pathway") %>%
      rename(diff_score = logFC, p_val = P.Value, padj = adj.P.Val, t_stat = t) %>%
      arrange(p_val)
    
    write_csv(res_sjs, "results/ucell_limma/limma_ucell_sjs_vs_ctrl.csv")
  }
}

message("Análise UCell concluída. Resultados em results/ucell_limma/")

library(vroom)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork) # Para combinar os gráficos

# --- 1. Carregamento e Preparação ---
# Substitua pelo nome correto do seu arquivo CSV do UCell
ucell_raw <- vroom("UCell/ucell_scores.csv") 
colnames(ucell_raw)[1] <- "ID"

# Limpar nomes das colunas (remover o sufixo _UCell)
ucell_clean <- ucell_raw
colnames(ucell_clean) <- gsub("_UCell$", "", colnames(ucell_clean))

# Carregar Metadados dos Clusters
clusters_meta <- vroom("data/precisesads/sjogren_clusters.csv") %>%
  rename(Cluster = PREDICTION) %>%
  mutate(Cluster = factor(paste0("Cluster ", Cluster)))

# Sincronizar dados
common_ids <- intersect(ucell_clean$ID, clusters_meta$ID)
ucell_data <- ucell_clean %>% filter(ID %in% common_ids) %>% arrange(match(ID, common_ids))
clusters_meta <- clusters_meta %>% filter(ID %in% common_ids) %>% arrange(match(ID, common_ids))

# Matriz apenas com os scores
mat <- ucell_data %>% column_to_rownames("ID") %>% as.matrix()

# --- 2. Painel A: PCA (Agrupamento Global) ---
pca_res <- prcomp(mat, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("ID") %>%
  left_join(clusters_meta, by = "ID")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  theme_minimal() +
  labs(title = "A. Global Functional Space (PCA)", 
       x = "PC1", y = "PC2") +
  theme(legend.position = "none")

# --- 3. Painel B: Boxplots de Vias de Interesse ---
# Selecionamos uma via de Célula B (MS4A1) e uma de Apresentação (Inata)
target_pathways <- c("B cell receptor signaling pathway", 
                    "Antigen processing and presentation")

p_boxes <- ucell_data %>%
  select(ID, all_of(target_pathways)) %>%
  pivot_longer(-ID, names_to = "Pathway", values_to = "Score") %>%
  left_join(clusters_meta, by = "ID") %>%
  ggplot(aes(x = Cluster, y = Score, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  facet_wrap(~Pathway, scales = "free_y") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  theme_minimal() +
  labs(title = "B. Key Pathway Activities") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# --- 4. Painel C: Heatmap de Assinatura (Top Variação) ---
# Selecionar as 30 vias com maior variância entre clusters
top_vars <- apply(mat, 2, var)
top_pathways <- names(sort(top_vars, decreasing = TRUE))[1:30]
h_mat <- t(mat[, top_pathways])

# Escalonar para Z-score
h_mat_scaled <- t(apply(h_mat, 1, scale))
colnames(h_mat_scaled) <- rownames(mat)

# Cores e Anotação
anno_col <- clusters_meta %>% select(ID, Cluster) %>% column_to_rownames("ID")
ann_colors <- list(Cluster = c("Cluster 1" = "#E41A1C", "Cluster 2" = "#377EB8", 
                              "Cluster 3" = "#4DAF4A", "Cluster 4" = "#984EA3"))

# Gerar o Heatmap como um objeto grob para combinar
library(grid)
p_heatmap <- pheatmap(
  h_mat_scaled,
  annotation_col = anno_col,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  cluster_cols = TRUE,
  main = "C. Top 30 Differential KEGG Pathways",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  silent = TRUE
)$gtable

# --- 5. Combinar e Salvar ---
# Organiza os plots: PCA e Boxplots no topo, Heatmap embaixo
final_plot <- (p_pca | p_boxes) / wrap_elements(p_heatmap)

ggsave("results/ucell_exploratory_figure.png", final_plot, width = 12, height = 10, bg = "white")