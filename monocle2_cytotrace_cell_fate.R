library(Seurat)
library(dplyr)
library(monocle)
library(RColorBrewer)
library(cowplot)
library(harmony)
library(sva)
library(CytoTRACE)
library(ggplot2)
library(Matrix)
library(scater)
library(scibetR)
library(scran)

install.packages("devtools")
devtools::install_local("/Users/tongcen/Desktop/blast/CytoTRACE_0.3.3.tar")

 
getwd()
setwd('/Users/tongcen/R.studio/cell annotation/15_cluster_update/radicle_11subclusterfrom15clusters/monocle2/selected_cluster/24h')
install.packages(c("sparsesvd", "docopt"))
install.packages("/Users/tongcen/R.studio/RNA velocity/qlcMatrix", repos = NULL, type = "source")



#setwd('/hwfssz1/ST_EARTH/P20Z10200N0035/USER/lihuanjin/05.Arabidopsis/monocle2_20230623')
#thedata <- readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/lihuanjin/05.Arabidopsis/monocle3_20230622/rosette_combine-EPI_MES_VAS.rds")

THEDATA1 <- readRDS("/Users/tongcen/R.studio/cell annotation/15_cluster_update/radicle_11subclusterfrom15clusters/15clusters.harmony.radicle.adata_R.rds")
dim(THEDATA1)

#thedata <- subset(thedata,idents=c("9","10"))
thedata <- subset(THEDATA1, HAP=='24H')
dim(thedata)
selected_clusters <- c("4", "5", "6", "8", "1")
thedata <- subset(thedata, seurat_clusters %in% selected_clusters) #注释后要这样提取
dim(thedata)

# 新建路径
dir <- "monocle2"
dir.create(dir)
setwd(dir)
# setwd("..")

# pdf('Dimplot.pdf',width = 10, height = 8 )
# DimPlot(object = thedata, reduction = "tsne",label = TRUE, pt.size = 0.5) #caution
# DimPlot(object = thedata, reduction = "umap",label = TRUE, pt.size = 0.2)
# dev.off()

###提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(thedata@assays$RNA@counts), 'sparseMatrix')
dim(thedata@assays$RNA@counts)
dim(expr_matrix)
# expr_matrix <- as(as.matrix(thedata@assays$SCT@counts), 'sparseMatrix')
##提取表型信息到p data(phenotype data)里面
p_data <- thedata@meta.data
dim(p_data)
p_data$celltype <- thedata@active.ident
p_data$celltype
row.names(thedata)
DefaultAssay(thedata) <- "RNA"
DefaultAssay(thedata)
f_data <- data.frame(gene_short_name = row.names(thedata),row.names = row.names(thedata))
dim(f_data)
#创建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
head(expr_matrix)
head(fd)

#size factor标准化细胞之间的mRNA差异
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

saveRDS(cds,"trajectory.cds")

#num_gene, filter
cds <- detectGenes(cds,min_expr = 0.1) #在pd和fd那里加上两列num_cells_expressed
print(head(fData(cds)))
length(fData(cds)$gene_short_name)
expressed_genes <- row.names(subset(fData(cds), 
                                    num_cells_expressed >= 10)) #在小于10个细胞中表达的基因会被去掉,值越大保留的基因越少？
length(fData(cds)$gene_short_name)

#使用seurat选择的高变基因
# expressed_genes<- VariableFeatures(pbmc)
# cds<- setOrderingFilter(cds,expressed_genes)
# plot_ordering_genes(cds)

#使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(pbmc)
# expressed_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds,expressed_genes)
# plot_ordering_genes(cds)

#使用monocle选择的高变基因
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >=1 * dispersion_fit)$gene_id
# cds <- setOrderingFilter(cds,disp.genes)
# plot_ordering_genes(cds)

#使用自定义发育marker基因

thedata <- FindVariableFeatures(thedata)
expressed_genes <- VariableFeatures(thedata)
length(expressed_genes)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~seurat_clusters",cores=1)
#差异表达基因作为轨迹构建基因，差异基因的选择标准是qval<0.01, decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) #默认是0.01
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
length(deg)
write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F, sep="\t", quote=F)
#轨迹构建基因可视化
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
#得到想要的基因列表后，需要使用setOrderingFilter嵌入cds对象，后续的操作都是依赖这个list
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]
table(cds@featureData@data[["use_for_ordering"]])

pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是第二步计算的基因表达量大小和离散度分布的趋势（找到的基因属于离散度比较高的基因）

#select top400 genes
ordergene <- row.names(deg)[order(deg$qval)][1:400]
#reduction
cds <- reduceDimension(cds,max_components = 2, method = "DDRtree")
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 2) #先在第一轮的结果看一下是否符合预期结果,然后再手动设置
pdf("trajectory_result.pdf", width = 8, height = 8)
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "seurat_clusters")
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()
saveRDS(cds,"trajectory.cds")

cds <- readRDS("trajectory.cds")

#其他的类似
pdf("trajectory_result_cell_type_HAP.pdf", width = 12, height = 15)
  plot_cell_trajectory(cds, color_by = "HAP") +
      facet_wrap(~seurat_clusters, nrow = 4) #设置几行几列展示
dev.off()

pdf("trajectory_result_HAP.pdf", width = 12, height = 15)
plot_cell_trajectory(cds, color_by = "HAP") +
  facet_wrap(~HAP, nrow = 4) #设置几行几列展示
dev.off()

#寻找拟时相关的基因
#heatpmap of DEG
#这里是把排序基因 (ordergene) 提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] 
write.csv(Time_diff,"Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
head(Time_genes)
p=plot_pseudotime_heatmap(cds[head(Time_genes,200),], 
                          num_clusters=4, 
                          show_rownames=T, 
                          return_heatmap=T)
q=plot_pseudotime_heatmap(cds[Time_genes], 
                          num_clusters=4, 
                          show_rownames=T, 
                          return_heatmap=T)

ggsave("Time_all_heatmapAll.pdf", q, width = 5, height = 45)
ggsave("Time_all_heatmapAll_short.pdf", q, width = 5, height = 5)

ggsave("Time_heatmapAll.pdf", p, width = 5, height = 45)
ggsave("Time_heatmapAll_short.pdf", p, width = 5, height = 5)

p$tree_row
clusters <- cutree(p$tree_row, k=4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv (clustering, "Time_clustering_all.csv")

q$tree_row
clusters <- cutree(q$tree_row, k=4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv (clustering, "Time_all_clustering_all.csv")

#显著差异基因按热图结果排序保存
hp.gene <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.gene, c("gene_short_name", "pval","qval")]
write.csv(Time_diff_sig, "Time_diff_sig.csv", row.names = F)

#显著差异基因按热图结果排序保存
hp.gene <- q$tree_row$labels[q$tree_row$order]
Time_diff_sig <- Time_diff[hp.gene, c("gene_short_name", "pval","qval")]
write.csv(Time_diff_sig, "Time_800_diff_sig.csv", row.names = F)

#基因双向图


pdf("heatmap_2direction_from4.pdf", width = 10, height = 50)
plot_genes_branched_heatmap(cds[Time_genes],
                            branch_states = c(1,4),#选择要比较的cluster
                            num_clusters = 6,
                            cores = 1, 
                            show_rownames = T,
                            return_heatmap = T)
dev.off()

p <- plot_genes_branched_heatmap(cds[Time_genes],
                                   branch_states = c(1,4),#选择要比较的cluster
                                   num_clusters = 6,
                                   cores = 1, 
                                   show_rownames = T,
                                   return_heatmap = T)

p$ph_res$tree_row$labels

p$annotation_col
p$annotation_row

gene_group=p$annotation_row
gene_group$gene=rownames(gene_group)
gene_group$gene

str(p,2)
# pdf("heatmap_2direction_from2.pdf", width = 10, height = 50)
# plot_genes_branched_heatmap(cds[Time_genes],
#                             branch_point = 2, #选择分叉点
#                             num_clusters = 5,
#                             cores = 4,
#                             show_rownames = T)
# dev.off()

a = c("HORVU.MOREX.r3.6HG0564480", "HORVU.MOREX.r3.3HG0234040")
plot_genes_in_pseudotime(a, color_by = "seurat_clusters")


pdf("trajectory_gene.pdf", width = 8, height = 8)
plot_cell_trajectory(cds,color_by = "Pseudotime") #+scale_color_viridis_c()) ????ɫ????
plot_cell_trajectory(cds,color_by = "seurat_clusters")
plot_cell_trajectory(cds,markers = a)
dev.off()

#指定基因绘制热图
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("HORVU.MOREX.r3.4HG0417260","HORVU.MOREX.r3.7HG0663130")))
diff_test_res <- differentialGeneTest(cds[c("HORVU.MOREX.r3.4HG0410930", "HORVU.MOREX.r3.7HG0741770"),],
                                      fullModelFormulaStr ="~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[a,],
                        num_clusters = 2,
                        cores =1,
                        show_rownames = T)


head(Time_genes)                          
# pdf("trajectory_gene.pdf", width = 8, height = 8)
# plot_cell_trajectory(cds,color_by = "Pseudotime") #+scale_color_viridis_c()) ????ɫ????
# plot_cell_trajectory(cds,color_by = "seurat_clusters")
# plot_cell_trajectory(cds,markers = "Zm00001eb055920")
# dev.off()

library(clusterProfiler)
library("org.BMorex.V3.eg.db")
##GO enrichment of cluster gene

gene_group=p$annotation_row
gene_group$gene=rownames(gene_group)
gene_group$gene
gene_group$Cluster

# all_clusters_800genes <- read.csv("Time_800genes_clustering_all.csv", header = TRUE, sep = ",")
# head(all_clusters_800genes)

cluster_gene_list <- list()
df <- list()
gene_group$Cluster <- as.character(gene_group$Cluster)
for (k in unique(gene_group$Cluster)[1:7]){
  print(k)
  # GO
  cluster_gene_list[[k]] <- gene_group[gene_group$Cluster == k,]
  # print(cluster_gene_list[[i]]$gene)
  df[[k]]=enrichGO(gene=as.vector(unique(cluster_gene_list[[k]]$gene)), 
                   OrgDb=org.BMorex.V3.eg.db,
                   keyType = "GID",
                   ont = "ALL",
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)
  
  plot <- dotplot(df[[k]], title = paste("cluster_", k, sep=""), font.size=8, split="ONTOLOGY", label_format=50) + facet_grid(ONTOLOGY ~ ., scales = "free")
  ggsave(paste0("monocle_cluster_", k, ".png"), plot, width = 10, height = 8, units = "in")
  write.csv(df[[k]], file = paste0("monocle_cluster_", k, ".csv"))
}


##使用cytoTrace进行细胞亚群分化潜能打分
expr <- as.matrix(thedata@assays$RNA@counts)
results <- CytoTRACE(mat = expr, enableFast = F)
pheno <- thedata@meta.data$SCT_snn_res.0.5
pheno <- as.character(pheno)
names(pheno) <- thedata@assays$RNA@counts@Dimnames[[2]]
options(repr.plot.height = 7.5, repr.plot.width = 16)
pdf("cytotrace.pdf", height = 10, width = 10)
plotCytoGenes(results, numOfGenes = 10)
plotCytoTRACE(results, phenotype = pheno)
dev.off()


