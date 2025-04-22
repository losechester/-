library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(SCP)
library(slingshot)

# 分出oligo子集，并进行PCA聚类降维
oligo_obj <- subset(sc_obj, subset = cell_type == "oligo")
all_genes <- rownames(oligo_obj)
oligo_obj <- ScaleData(oligo_obj, features = all_genes)
oligo_obj <- FindVariableFeatures(oligo_obj, selection.method = "vst", nfeatures = 2000)
oligo_obj <- RunPCA(oligo_obj, features = VariableFeatures(object = oligo_obj))
oligo_obj <- FindNeighbors(oligo_obj, dims = 1:10)
oligo_obj <- FindClusters(oligo_obj, resolution = 0.5)
oligo_obj <- RunUMAP(oligo_obj, dims = 1:10)

# 分出astro子集，并进行PCA聚类降维
astro_obj <- subset(sc_obj, subset = cell_type == "astro")
all_genes <- rownames(astro_obj)
astro_obj <- ScaleData(astro_obj, features = all_genes)
astro_obj <- FindVariableFeatures(astro_obj, selection.method = "vst", nfeatures = 1000)
astro_obj <- RunPCA(astro_obj, features = VariableFeatures(object = astro_obj))
astro_obj <- FindNeighbors(astro_obj, dims = 1:10)
astro_obj <- FindClusters(astro_obj, resolution = 0.5)
astro_obj <- RunUMAP(astro_obj, dims = 1:10)

# 按细胞亚群、是否患病分别聚类上色
p2 <- CellDimPlot(
  #srt = astro_obj,
  #srt = oligo_obj,
  srt = sc_obj,
  group.by = c(
    "seurat_clusters"
    #"condition"
    #"subclust"
    #"cell_type"
  ),
  #label = TRUE,
  reduction = "UMAP",
 # reduction = "harmony",
  legend.position = "right",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  #title = "GSE138852_oligo"
)
p2

#拟时序分析，计算分化轨迹
reduced_data <- Embeddings(astro_obj, reduction = "umap")
clusters <- astro_obj$subclust
sling <- slingshot(reduced_data, clusters)
pseudotime <- slingPseudotime(sling)
for (i in seq_len(ncol(pseudotime))) {  # 用列数来循环
  col_name <- paste0("Lineage", i)
  astro_obj[[col_name]] <- pseudotime[, i]  # 取每一条 lineage 的 pseudotime
}


# 获取轨迹信息
curves <- slingCurves(sling)
p2 <- FeatureDimPlot(
  astro_obj,
  features = paste0("Lineage", 1),
  # lineages = paste0("Lineage", 1),
  reduction = "harmony",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  subtitle = ""
)

# 按需要组合生成图片
library(patchwork)
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# 显示组合后的图
print(combined_plot)

###oligo
reduced_data <- Embeddings(oligo_obj, reduction = "umap")
clusters <- oligo_obj$subclust
sling <- slingshot(reduced_data, clusters)
pseudotime <- slingPseudotime(sling)
for (i in seq_len(ncol(pseudotime))) {  # 用列数来循环
  col_name <- paste0("Lineage", i)
  oligo_obj[[col_name]] <- pseudotime[, i]  # 取每一条 lineage 的 pseudotime
}

# 获取轨迹信息
curves <- slingCurves(sling)
p4 <- FeatureDimPlot(
  oligo_obj,
  features = paste0("Lineage", 2),
  # lineages = paste0("Lineage", 1),
  reduction = "harmony",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  #subtitle = "GSE138852_oligo"
  subtitle = " "
)
p10

# 特征基因计算，目前没实现
astro_obj <- RunDynamicFeatures(
  srt = astro_obj,
  lineages = c("Lineage1"),
  n_candidates = 100
)

features_sel <- names(astro_obj@tools$DynamicFeatures$family)
astro_obj<- RunDynamicFeatures(
  srt = astro_obj@tools,
  features = c(
    features_sel[1:49],
    "ADAM32"
  ),
  lineages = c("Lineage1"),
  n_candidates = 100
)

ht2 <- DynamicHeatmap(
  srt = astro_obj,
  # features = c("ADAM32"),
  lineages = c("Lineage1"),
  show_row_names = TRUE,
  min_expcells = 0,
  r.sq = 0,
  dev.expl = 0,
  padjust = 1,
  # use_fitted = TRUE,
  n_split = 3,
  # reverse_ht = "Lineage1",
  species = "Homo_sapiens",
  # db = "GO_BP",
  db = "GO_BP",
  anno_terms = TRUE,
  # anno_keys = TRUE,
  # anno_features = TRUE,
  heatmap_palette = "viridis",
  cell_annotation = "subclust",
  separate_annotation = list("subclust"), # , c("ADAM32")
  separate_annotation_palette = c("Paired"), # , "Set1"
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  height = 7,
  width = 3
)
ht2$plot

