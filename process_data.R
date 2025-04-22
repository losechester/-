if (F) {
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("zhanghao-njmu/SCP")
  devtools::install_local("D:/R-4.4.3/library/SCP",force = TRUE)
  
  remove.packages(c("Seurat", "SeuratObject"))
  install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'))
  install.packages('SeuratObject')
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(SCP)

# 读取数据
counts <- read.csv("D:/test4scp/GSE138852_counts.csv", row.names = 1)
metadata_path <- "D:/test4scp/GSE138852.csv.gz"
metadata <- read.csv(metadata_path, row.names = 1)
head(metadata)

dim(counts)    # 行：基因数；列：细胞数
dim(metadata) 

# 2. 构建 Seurat 对象
# 将表达矩阵和元数据结合，创建 Seurat 对象
sc_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)

# 3. 质量控制 (QC)
# 示例：计算线粒体基因百分比
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")

# 可视化 QC 指标
VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 根据实际数据选择合适的过滤阈值
# 保留检测到的基因数在 200 ~ 2500 之间，且线粒体百分比低于 5%
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

sc_obj$cell_type <- metadata[colnames(sc_obj), "oupSample.cellType"]
sc_obj$condition <- metadata[colnames(sc_obj), "oupSample.batchCond"]
sc_obj$subclust <- metadata[colnames(sc_obj), "oupSample.subclustID"]

VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 4. 数据归一化与高变基因识别
# 对数据进行归一化（LogNormalize 方法，每个细胞总表达量归一化到 10,000）
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 寻找高变基因（选择 2000 个）
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)

# 查看前10个高变基因
top10 <- head(VariableFeatures(sc_obj), 10)
print(top10)
my_colors <- c( "pink3", "brown")

# 绘制可变特征图
plot1 <- VariableFeaturePlot(sc_obj, cols = my_colors)
# 可视化高变基因
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# 5. 数据缩放和主成分分析 (PCA)
# 对所有基因进行数据缩放
all_genes <- rownames(sc_obj)
sc_obj <- ScaleData(sc_obj, features = all_genes)

# 进行 PCA 分析，使用高变基因
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
# 查看 PCA 结果（前几个主成分的主要基因）
print(sc_obj[["pca"]], dims = 1:5, nfeatures = 5)

# 可视化 PCA 负载
VizDimLoadings(sc_obj, dims = 1:2, reduction = "pca")
DimPlot(sc_obj, reduction = "pca")
  
# 绘制肘部图判断需要保留的主成分数
ElbowPlot(sc_obj)

# 6. 聚类分析
# 基于 PCA 结果构建邻近图（使用前10个主成分，根据肘部图调整）
sc_obj <- FindNeighbors(sc_obj, dims = 1:10)
# 进行聚类（分辨率参数可调，通常在 0.4~1.2 之间）
sc_obj <- FindClusters(sc_obj, resolution = 0.5)

# 查看聚类结果（保存在 meta.data 中，默认为 "seurat_clusters" 列）
table(sc_obj$seurat_clusters)
DimPlot(sc_obj, reduction = "pca")
# 绘制热图
DimHeatmap(sc_obj, dims = 1:15, cells = 500, balanced = TRUE) 

# 7. 降维可视化：UMAP/tSNE
# 运行 UMAP 降维（同样使用前10个主成分）
sc_obj <- RunUMAP(sc_obj, dims = 1:10)
# 绘制 UMAP 图，并显示细胞聚类
DimPlot(sc_obj, reduction = "umap", label = TRUE, label.size = 5) + ggtitle("UMAP 分布")

DimPlot(sc_obj, reduction = "umap",group.by="cell_type")+ ggtitle("UMAP 分布")
# 或者也可以运行 tSNE：
sc_obj <- RunTSNE(sc_obj, dims = 1:10)
DimPlot(sc_obj, reduction = "tsne",group.by="cell_type") + ggtitle("tSNE 分布")

# 8. 保存结果
# 保存 Seurat 对象以便后续分析
saveRDS(sc_obj, file = "D:/test4scp/sc_obj.rds")
# 活动身份改为cell_type
Idents(sc_obj) <- "cell_type"
markers <- FindAllMarkers(sc_obj, group.by = "cell_type", only.pos = TRUE, 
                          min.pct = 0.25, logfc.threshold = 0.25)

# 为所有聚类寻找marker基因（只保留上调基因）(未定义的)
#markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 查看每个聚类top marker基因

top_markers <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
head(top_markers)

# 用热图展示各细胞群的top marker表达
DoHeatmap(sc_obj,features = top_markers$gene,size = 3,angle = 0,label = F)
