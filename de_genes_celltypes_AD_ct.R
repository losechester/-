library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggsci)

# 只有astro和oligo细胞有单独的GO 分析代码
# 1. 处理 astro 细胞类型
astro_obj <- subset(sc_obj, subset = cell_type == "astro")
Idents(astro_obj) <- "oupSample.batchCond"
de_genes_astro <- FindMarkers(astro_obj, ident.1 = "AD", ident.2 = "ct", 
                              min.pct = 0.1, logfc.threshold = 0.25)
volcano_data_astro <- data.frame(
  gene = rownames(de_genes_astro),
  log2FC = de_genes_astro$avg_log2FC,
  pvalue = de_genes_astro$p_val_adj
)
volcano_data_astro$significance <- ifelse(volcano_data_astro$log2FC > 0 & volcano_data_astro$pvalue < 0.05, "Up",
                                          ifelse(volcano_data_astro$log2FC < 0 & volcano_data_astro$pvalue < 0.05, "Down", "Not significant"))
volcano_data_astro$to_label <- (volcano_data_astro$pvalue < 0.05) & (abs(volcano_data_astro$log2FC) > 1)


# 选择最显著的前15个上调和下调基因
top_up_astro <- volcano_data_astro %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_astro <- volcano_data_astro %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_astro <- rbind(top_up_astro, top_down_astro)

# 修改绘图代码
library(viridis)  # 使用 viridis 调色板

# 美化火山图，颜色渐变
p_astro <- ggplot(volcano_data_astro, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = log2FC), size = 2) +  # 按 log2FC 渐变颜色，增大点大小
  scale_color_viridis(option = "magma", name = "Log2 Fold Change", 
                      direction = -1) +  # 使用 magma 调色板，负值到正值渐变
  geom_text_repel(data = to_label_genes_astro, aes(label = gene), 
                  size = 3, max.overlaps = 30, 
                  box.padding = 0.5, point.padding = 0.5, 
                  segment.color = "grey50") +  # 优化标签排斥
  labs(x = "Log2 Fold Change", 
       y = "-Log10(Adjusted p-value)", 
       title = "Volcano Plot of Differentially Expressed Genes in Astro Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.1, face = "bold", size = 14),  # 美化标题
        axis.title = element_text(size = 12),  # 调整轴标题
        axis.text = element_text(size = 10),  # 调整轴文本
        legend.position = "right")  # 图例位置

# 打印火山图
print(p_astro)

# 提取显著差异表达基因
sig_genes <- volcano_data_astro %>% 
  filter(significance != "Not significant") %>% 
  pull(gene)

# 将基因符号转换为 Entrez ID
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = sig_genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

############new
ego<-enrichGO(gene=entrez_ids,
              OrgDb=org.Hs.eg.db, 
              keyType='ENTREZID', 
              ont="ALL", 
              pAdjustMethod="BH",  #对P值进行矫正
              pvalueCutoff=0.05,  #pvalue<0.05
              qvalueCutoff=0.2,   #qvalue<0.2
              readable=TRUE)

ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID') #转换ENTREZID到Symbol
nrow(ego) #富集到的条目数量
sum(ego$ONTOLOGY=="BP") #Biological process
sum(ego$ONTOLOGY=="CC") #Cellular component
sum(ego$ONTOLOGY=="MF") #Molecular function
write.csv(as.data.frame(ego),"GO_enrich_astro.csv",row.names =F)  #输出为表格保存

p<-barplot(ego,
        #title='GO Enrichment Analysis of Astro Cells', #可设置标题
        color='p.adjust', 
        showCategory=5,  #显示前几条
        font.size = 10, #条目字体大小
        split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")+  #按照BP、CC、MF分隔图片
  scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 60)) +#设置条形图框比例
ggtitle("GO Enrichment Analysis for Astro Cells") +  # 添加标题
theme(plot.title = element_text(hjust = 0.9, face = "bold", size = 14))
#p+ scale_fill_distiller(palette = "Blues", direction = 1)
p + viridis::scale_fill_viridis(option = "plasma", begin = 0.2, end = 1)
cnetplot(ego,
         showCategory = 5, #显示前几条
         circular = T,  #圆形设置为TURE
         colorEdge = TRUE,
         node_label ="category", ## category | gene | all | none
         #cex_category = 5,
         cex_gene = 1,
         cex_label_category = 1)


# 进行 GO 富集分析
#ego <- enrichGO(gene = entrez_ids,
#OrgDb = org.Hs.eg.db,
 #               ont = "BP",  # 可以选择 BP（生物学过程）、MF（分子功能）、CC（细胞组分）
  #              pAdjustMethod = "BH",
   #             qvalueCutoff = 0.05)

# 输出 GO 富集分析结果
#print(ego)

# 可视化 GO 富集分析结果（点图）
dotplot(ego, showCategory = 10)
ggsave("GO_enrichment_dotplot.png", width = 8, height = 8, dpi = 300)

# 2. 处理 doublet 细胞类型
doublet_obj <- subset(sc_obj, subset = cell_type == "doublet")
Idents(doublet_obj) <- "oupSample.batchCond"
de_genes_doublet <- FindMarkers(doublet_obj, ident.1 = "AD", ident.2 = "ct", 
                                min.pct = 0.1, logfc.threshold = 0.25)
volcano_data_doublet <- data.frame(
  gene = rownames(de_genes_doublet),
  log2FC = de_genes_doublet$avg_log2FC,
  pvalue = de_genes_doublet$p_val_adj
)
volcano_data_doublet$significance <- ifelse(volcano_data_doublet$log2FC > 0 & volcano_data_doublet$pvalue < 0.05, "Up",
                                            ifelse(volcano_data_doublet$log2FC < 0 & volcano_data_doublet$pvalue < 0.05, "Down", "Not significant"))
volcano_data_doublet$to_label <- (volcano_data_doublet$pvalue < 0.05) & (abs(volcano_data_doublet$log2FC) > 1)
# 选择最显著的前15个上调和下调基因
top_up_doublet <- volcano_data_doublet %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_doublet <- volcano_data_doublet %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_doublet <- rbind(top_up_doublet, top_down_doublet)

# 修改绘图代码
p_doublet <- ggplot(volcano_data_doublet, aes(x = log2FC, y = -log10(pvalue), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "gray")) +
  geom_text_repel(data = to_label_genes_doublet, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)", title = "Volcano Plot of Differentially Expressed Genes in doublet cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste0("Volcano_doublet.png"),
  plot = p_doublet,
  width = 8,
  height = 10,
  dpi = 300
)

# 3. 处理 endo 细胞类型
endo_obj <- subset(sc_obj, subset = cell_type == "endo")
Idents(endo_obj) <- "oupSample.batchCond"
de_genes_endo <- FindMarkers(endo_obj, ident.1 = "AD", ident.2 = "ct", 
                             min.pct = 0.1, logfc.threshold = 0.25)
volcano_data_endo <- data.frame(
  gene = rownames(de_genes_endo),
  log2FC = de_genes_endo$avg_log2FC,
  pvalue = de_genes_endo$p_val_adj
)
volcano_data_endo$significance <- ifelse(volcano_data_endo$log2FC > 0 & volcano_data_endo$pvalue < 0.05, "Up",
                                         ifelse(volcano_data_endo$log2FC < 0 & volcano_data_endo$pvalue < 0.05, "Down", "Not significant"))
volcano_data_endo$to_label <- (volcano_data_endo$pvalue < 0.05) & (abs(volcano_data_endo$log2FC) > 1)
# 选择最显著的前15个上调和下调基因
top_up_endo <- volcano_data_endo %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_endo <- volcano_data_endo %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_endo <- rbind(top_up_endo, top_down_endo)

# 修改绘图代码
p_endo <- ggplot(volcano_data_endo, aes(x = log2FC, y = -log10(pvalue), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "gray")) +
  geom_text_repel(data = to_label_genes_endo, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)", title = "Volcano Plot of Differentially Expressed Genes in endo cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste0("Volcano_endo.png"),
  plot = p_endo,
  width = 8,
  height = 10,
  dpi = 300
)

# 4. 处理 mg 细胞类型
mg_obj <- subset(sc_obj, subset = cell_type == "mg")
Idents(mg_obj) <- "oupSample.batchCond"
de_genes_mg <- FindMarkers(mg_obj, ident.1 = "AD", ident.2 = "ct", 
                           min.pct = 0.1, logfc.threshold = 0.25)
volcano_data_mg <- data.frame(
  gene = rownames(de_genes_mg),
  log2FC = de_genes_mg$avg_log2FC,
  pvalue = de_genes_mg$p_val_adj
)
volcano_data_mg$significance <- ifelse(volcano_data_mg$log2FC > 0 & volcano_data_mg$pvalue < 0.05, "Up",
                                       ifelse(volcano_data_mg$log2FC < 0 & volcano_data_mg$pvalue < 0.05, "Down", "Not significant"))
volcano_data_mg$to_label <- (volcano_data_mg$pvalue < 0.05) & (abs(volcano_data_mg$log2FC) > 1)
# 选择最显著的前15个上调和下调基因
top_up_mg <- volcano_data_mg %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_mg <- volcano_data_mg %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_mg <- rbind(top_up_mg, top_down_mg)

# 修改绘图代码
p_mg <- ggplot(volcano_data_mg, aes(x = log2FC, y = -log10(pvalue), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "gray")) +
  geom_text_repel(data = to_label_genes_mg, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)", title = "Volcano Plot of Differentially Expressed Genes in mg cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste0("Volcano_mg.png"),
  plot = p_mg,
  width = 8,
  height = 10,
  dpi = 300
)

# 5. 处理 neuron 细胞类型
neuron_obj <- subset(sc_obj, subset = cell_type == "neuron")
Idents(neuron_obj) <- "oupSample.batchCond"
de_genes_neuron <- FindMarkers(neuron_obj, ident.1 = "AD", ident.2 = "ct", 
                               min.pct = 0.1, logfc.threshold = 0.25)
volcano_data_neuron <- data.frame(
  gene = rownames(de_genes_neuron),
  log2FC = de_genes_neuron$avg_log2FC,
  pvalue = de_genes_neuron$p_val_adj
)
volcano_data_neuron$significance <- ifelse(volcano_data_neuron$log2FC > 0 & volcano_data_neuron$pvalue < 0.05, "Up",
                                           ifelse(volcano_data_neuron$log2FC < 0 & volcano_data_neuron$pvalue < 0.05, "Down", "Not significant"))
volcano_data_neuron$to_label <- (volcano_data_neuron$pvalue < 0.05) & (abs(volcano_data_neuron$log2FC) > 1)
# 选择最显著的前15个上调和下调基因
top_up_neuron <- volcano_data_neuron %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_neuron <- volcano_data_neuron %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_neuron <- rbind(top_up_neuron, top_down_neuron)

# 修改绘图代码
p_neuron <- ggplot(volcano_data_neuron, aes(x = log2FC, y = -log10(pvalue), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "gray")) +
  geom_text_repel(data = to_label_genes_neuron, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)", title = "Volcano Plot of Differentially Expressed Genes in neuron cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste0("Volcano_neuron.png"),
  plot = p_neuron,
  width = 8,
  height = 10,
  dpi = 300
)
########################################################################
# 6. 处理 oligo 细胞类型
oligo_obj <- subset(sc_obj, subset = cell_type == "oligo")
Idents(oligo_obj) <- "oupSample.batchCond"
de_genes_oligo <- FindMarkers(oligo_obj, ident.1 = "AD", ident.2 = "ct", 
                              min.pct = 0.1, logfc.threshold = 0.25)
volcano_data_oligo <- data.frame(
  gene = rownames(de_genes_oligo),
  log2FC = de_genes_oligo$avg_log2FC,
  pvalue = de_genes_oligo$p_val_adj
)

volcano_data_oligo$significance <- ifelse(volcano_data_oligo$log2FC > 0 & volcano_data_oligo$pvalue < 0.05, "Up",
                                          ifelse(volcano_data_oligo$log2FC < 0 & volcano_data_oligo$pvalue < 0.05, "Down", "Not significant"))
volcano_data_oligo$to_label <- (volcano_data_oligo$pvalue < 0.05) & (abs(volcano_data_oligo$log2FC) > 1)


# 选择最显著的前15个上调和下调基因
top_up_oligo <- volcano_data_oligo %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_oligo <- volcano_data_oligo %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_oligo <- rbind(top_up_oligo, top_down_oligo)

p_oligo <- ggplot(volcano_data_oligo, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = log2FC), size = 2) +  # 按 log2FC 渐变颜色，增大点大小
  scale_color_viridis(option = "viridis", name = "Log2 Fold Change", 
                      direction = -1) +  
  geom_text_repel(data = to_label_genes_oligo, aes(label = gene), 
                  size = 3, max.overlaps = 30, 
                  box.padding = 0.5, point.padding = 0.5, 
                  segment.color = "grey50") +  # 优化标签排斥
  labs(x = "Log2 Fold Change", 
       y = "-Log10(Adjusted p-value)",
       title = "Volcano Plot of Differentially Expressed Genes in Oligo Cells") +
   theme_minimal() +
  theme(plot.title = element_text(hjust = 0.1, face = "bold", size = 14),  # 美化标题
        axis.title = element_text(size = 12),  # 调整轴标题
        axis.text = element_text(size = 10),  # 调整轴文本
        legend.position = "right")  # 图例位置 + 
coord_cartesian(ylim = c(0, 300))

# 打印火山图
print(p_oligo)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggsci)
# 提取显著差异表达基因
sig_genes <- volcano_data_oligo %>% 
  filter(significance != "Not significant") %>% 
  pull(gene)

# 将基因符号转换为 Entrez ID
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = sig_genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

############new
ego<-enrichGO(gene=entrez_ids,
              OrgDb=org.Hs.eg.db, 
              keyType='ENTREZID', 
              ont="ALL", 
              pAdjustMethod="BH",  #对P值进行矫正
              pvalueCutoff=0.05,  #pvalue<0.05
              qvalueCutoff=0.2,   #qvalue<0.2
              readable=TRUE)

ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID') #转换ENTREZID到Symbol
nrow(ego) #富集到的条目数量
sum(ego$ONTOLOGY=="BP") #Biological process
sum(ego$ONTOLOGY=="CC") #Cellular component
sum(ego$ONTOLOGY=="MF") #Molecular function
write.csv(as.data.frame(ego),"GO_enrich_oligo.csv",row.names =F)  #输出为表格保存

p<-barplot(ego,
           #title='GO Enrichment Analysis of Astro Cells', #可设置标题
           color='p.adjust', 
           showCategory=5,  #显示前几条
           font.size = 10, #条目字体大小
           split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")+  #按照BP、CC、MF分隔图片
  scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 60)) +#设置条形图框比例
  ggtitle("GO Enrichment Analysis for Oligo Cells") +  # 添加标题
  theme(plot.title = element_text(hjust = 0.9, face = "bold", size = 14))
#p+ scale_fill_distiller(palette = "Blues", direction = 1)
p + viridis::scale_fill_viridis(option = "viridis", begin = 0.6, end = 1)
# 修改绘图代码
p_oligo <- ggplot(volcano_data_oligo, aes(x = log2FC, y = -log10(pvalue), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "gray")) +
  geom_text_repel(data = to_label_genes_oligo, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)", title = "Volcano Plot of Differentially Expressed Genes in oligo cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste0("Volcano_oligo.png"),
  plot = p_oligo,
  width = 8,
  height = 10,
  dpi = 300
)


# 7.OPC
OPC_obj <- subset(sc_obj, subset = cell_type == "OPC")

# 设置细胞身份标识
Idents(OPC_obj) <- "oupSample.batchCond"

# 进行差异表达分析
de_genes_OPC <- FindMarkers(OPC_obj, ident.1 = "AD", ident.2 = "ct", 
                            min.pct = 0.1, logfc.threshold = 0.25)

# 整理火山图数据
volcano_data_OPC <- data.frame(
  gene = rownames(de_genes_OPC),
  log2FC = de_genes_OPC$avg_log2FC,
  pvalue = de_genes_OPC$p_val_adj
)

# 添加显著标记
volcano_data_OPC$significance <- ifelse(volcano_data_OPC$log2FC > 0 & volcano_data_OPC$pvalue < 0.05, "Up",
                                        ifelse(volcano_data_OPC$log2FC < 0 & volcano_data_OPC$pvalue < 0.05, "Down", "Not significant"))

# 标记需要标注的基因
volcano_data_OPC$to_label <- (volcano_data_OPC$pvalue < 0.05) & (abs(volcano_data_OPC$log2FC) > 1)

# 选择最显著的前15个上调和下调基因
top_up_OPC <- volcano_data_OPC %>% 
  filter(significance == "Up") %>% 
  arrange(pvalue) %>% 
  head(15)

top_down_OPC <- volcano_data_OPC %>% 
  filter(significance == "Down") %>% 
  arrange(pvalue) %>% 
  head(15)

# 合并需要标记的基因
to_label_genes_OPC <- rbind(top_up_OPC, top_down_OPC)

# 绘制火山图
p_OPC <- ggplot(volcano_data_OPC, aes(x = log2FC, y = -log10(pvalue), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "gray")) +
  geom_text_repel(data = to_label_genes_OPC, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)", title = "Volcano Plot of Differentially Expressed Genes in OPC cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# 保存火山图
ggsave(
  filename = paste0("Volcano_OPC.png"),
  plot = p_OPC,
  width = 8,
  height = 12,
  dpi = 300
)



#####################同时GO分析，并在一张图上

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 构建每个细胞类型的显著基因列表（以符号形式）
geneList <- list(
  astro   = volcano_data_astro %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene),
  doublet = volcano_data_doublet %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene),
  endo    = volcano_data_endo %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene),
  mg      = volcano_data_mg %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene),
  neuron  = volcano_data_neuron %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene),
  oligo   = volcano_data_oligo %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene),
  OPC     = volcano_data_OPC %>% filter(significance != "Not significant", abs(log2FC) > 1) %>% pull(gene)
)

# 将基因符号转换为 Entrez ID（去除转换失败的 NA）
geneList_entrez <- lapply(geneList, function(genes) {
  ids <- mapIds(org.Hs.eg.db,
                keys = genes,
                column = "ENTREZID",
                keytype = "SYMBOL",
                multiVals = "first")
  ids[!is.na(ids)]
})

# 使用 compareCluster 进行 GO 生物学过程（BP）富集分析
cc <- compareCluster(geneCluster = geneList_entrez, 
                     fun = "enrichGO",
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",             # 也可以选择 MF 或 CC
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

library(stringr)  # 加载字符串处理包

p_compare <- dotplot(cc, showCategory = 10) + 
  theme_bw() +
  ggtitle("GO Enrichment Analysis Across Cell Types") +
  # 添加标签换行设置，width 根据实际标签长度调整（如 20 表示每行约 20 字符）
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) 
# 保存合并后的富集分析图
ggsave("GO_compareCluster_dotplot.png", p_compare, width = 10, height = 12, dpi = 300)
