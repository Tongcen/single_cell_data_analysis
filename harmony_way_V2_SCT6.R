args<-commandArgs(T)
library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(cowplot)
library(harmony)
library(getopt)
arg <- matrix(c(
                "input_file","i",1,"character","input txt file,include count_mtx.tsv.gz",
                "input_sample","s",1,"character","input txt file,include corresponding sample_name",
                "Feature_min","n",1,"integer","nFeature_RNA > ",
                "Feature_max","a",1,"integer","nFeature_RNA < ",
                "Count_min","c",1,"integer","nCount_RNA > ",
                "PercentageFeatureSet1","p",1,"character","Calculate proportion for quality control1",
                "PercentageFeatureSet2","f",1,"character","Calculate proportion for quality control2",
                "percent1","q",1,"integer","selected corresponding portion1",
                "percent2","g",1,"integer","selected corresponding portion2",
                "RunUMAP","U",1,"integer","each sample RunUMAP",
                "FindNeighbors","b",1,"integer","Computing nearest neighbor graph",
                "lambda","l","1","numeric"," RUNharmony lambda, default=1",
                "resolution","r",1,"numeric","Find Clusters resolution",
                "tissue","t","1","character","tissue,default=Embro"
                ),byrow=T,ncol=5)


opt = getopt(arg)
if (is.null(opt$input_file)){
    opt$input_file <- "count_mtx.filter.lst"
}
if (is.null(opt$input_sample)){
    opt$input_sample <- "sample.filter.list"
}
if (is.null(opt$Feature_min)){
    opt$Feature_min <-200
}
if (is.null(opt$Feature_max)){
    opt$Feature_max <- 3000
}
if (is.null(opt$Count_min)){
    opt$Count_min <- 1000
}
if (is.null(opt$PercentageFeatureSet1)){
    opt$PercentageFeatureSet1 <- "none"
}
if (is.null(opt$PercentageFeatureSet2)){
    opt$PercentageFeatureSet2 <- "none"
}
if (is.null(opt$percent1)){
    opt$percent1 <- 5
}
if (is.null(opt$percent2)){
    opt$percent2 <- 5
}
if (is.null(opt$RunUMAP)){
    opt$RunUMAP <- 30
}
if (is.null(opt$FindNeighbors)){
    opt$FindNeighbors <- 30
}
if (is.null(opt$resolution)){
    opt$FindNeighbors <- 1
}
if (is.null(opt$lambda)){
    opt$lambda <- 1
}
if (is.null(opt$tissue)){
        opt$tissue <- "Embryo"
}
print(opt)


Find_doublet <- function(data){
    sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    nExp_poi <- round(0.05*ncol(data))
    p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
    data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    data<-subset(data,subset=doublet_info=="Singlet")
    data
}


o.list<-list()
matrix_3_file <- read.table(opt$input_file)
sample_file <- read.table(opt$input_sample)
matrix_num <- length(matrix_3_file[[1]])
pdf('harmony.Vinplot_ElbowPlot.each.pdf',width = 16 , height = 24 )
for (i in 1:matrix_num) {
    mat <- Read10X(data.dir = matrix_3_file[[1]][i],gene.column = 1, cell.column = 1)
    rds<-CreateSeuratObject(mat,min.cells=3,min.features =200)
    #        rds[["percent.1"]] <- PercentageFeatureSet(rds, pattern = opt$PercentageFeatureSet1,assay = 'RNA')
    #        rds[["percent.2"]] <- PercentageFeatureSet(rds, pattern = opt$PercentageFeatureSet2,assay = 'RNA')
    #	p1 <- VlnPlot(rds, features = c("nFeature_RNA", "nCount_RNA", "percent.1","percent.2"), ncol = 4,pt.size = 0.1)
    p1 <- VlnPlot(rds, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0.1)
    #        rds <- subset(rds, subset = nFeature_RNA > opt$Feature_min & nFeature_RNA < opt$Feature_max & nCount_RNA > opt$Count_min & percent.1 < opt$percent1 & percent.2 < opt$percent2)
    rds <- subset(rds, subset = nFeature_RNA > opt$Feature_min & nFeature_RNA < opt$Feature_max & nCount_RNA > opt$Count_min ) 
    rds$batch <- sample_file[[1]][i]
    #        rds$organs <- gsub(pattern = "_iDrop_.*",replacement ="",x=sample_file[[1]][i])
    #        rds$sample<-gsub(pattern = ".*_iDrop",replacement ="iDrop",x=sample_file[[1]][i])


    rds <- NormalizeData(rds)
    rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 2000)
    rds <- ScaleData(rds)
    rds <- RunPCA(rds,npcs = 30)
    #p2 <- ElbowPlot(rds)
    #rds <- RunUMAP(rds, dims = 1:opt$RunUMAP)

    rds<-Find_doublet(rds)
    rds<-subset(rds,subset=doublet_info=="Singlet")
    o.list[[i]]<-rds
    print(p1)
    #print(p2)
}
dev.off()

adata<-merge(x=o.list[[1]],y=o.list[c(2:length(o.list))])

#adata<-NormalizeData(adata)
#adata<- FindVariableFeatures(adata, selection.method = "vst")
#adata<- ScaleData(adata)
adata<- SCTransform(adata, verbose = FALSE)
adata <- RunPCA(object = adata, npcs = 30, verbose = FALSE)
adata <-RunHarmony(adata,'batch',plot_convergenc=TRUE,lambda = opt$lambda, assay.use='SCT')

#DefaultAssay(adata)="RNA"
#all.genes <- rownames(adata)
#adata <- ScaleData(object = adata, features = all.genes, verbose = FALSE)

#adata <- RunPCA(object = adata, npcs = 50, verbose = FALSE)

pdf('harmony.ElbowPlot.pdf',width = 16 , height = 24)
ElbowPlot(adata)
dev.off()


adata <-FindNeighbors(object = adata,dims = 1:opt$FindNeighbors,reduction = "harmony")
adata <- FindClusters(object = adata,resolution = opt$resolution,reduction = "harmony")

pdf('harmony.Vinplot.all.pdf',width = 16 , height = 24 )
VlnPlot(adata, features = c("nFeature_RNA"), ncol = 1,pt.size = 0.1)
VlnPlot(adata, features = c("nCount_RNA"), ncol = 1,pt.size = 0.1)
dev.off()


pdf('harmony.Dimplot.pdf',width = 10, height = 8 )
adata <- RunUMAP(adata, reduction = "harmony", dims = 1:opt$FindNeighbors)
DimPlot(object = adata, reduction = "umap",label = TRUE, pt.size = 0.2)
DimPlot(object = adata, reduction = "umap",group.by = "batch", pt.size = 0.2)
#DimPlot(object = adata, reduction = "umap",group.by = "organs", pt.size = 0.2)

adata <- RunTSNE(adata, reduction = "harmony", dims = 1:opt$FindNeighbors)
DimPlot(object = adata, reduction = "tsne",label = TRUE, pt.size = 0.2)
DimPlot(object = adata, reduction = "tsne",group.by = "batch", pt.size = 0.2)
#DimPlot(object = adata, reduction = "tsne",group.by = "organs", pt.size = 0.2)
dev.off()

#pdf('harmony.Dimplot.batch.pdf',width = 32, height = 8 )
#DimPlot(object = adata, reduction = "umap",group.by = "batch", pt.size = 0.2)
#DimPlot(object = adata, reduction = "tsne",group.by = "batch", pt.size = 0.2)
#dev.off()
#
#pdf('harmony.Dimplot.split.pdf',width = 32, height = 8 )
#DimPlot(adata, reduction = "umap", split.by = "organs",label = TRUE)
#DimPlot(adata, reduction = "tsne", split.by = "organs",label = TRUE)
#dev.off()

#DefaultAssay(object = adata) <- "RNA"
adata.markers <- FindAllMarkers(adata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allmakers <- adata.markers %>% group_by(cluster)
write.csv(allmakers,file = "markergenes_list.harmony.csv")
pdf('harmony.top5_markergenes.pdf',width = 12, height = 12)
top5 <- adata.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- adata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(adata, features = top5$gene) + NoLegend()
dev.off()

pdf ("13_markergeneheatmap.pdf")
top50 <- adata.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(adata, features = top50$gene) + NoLegend() 
dev.off()
write.table(top50,file="13_markergeneheatmap.txt",sep="\t",quote=F,col.names = T,row.names = FALSE)

saveRDS(adata, file = "harmony.all_tissue.rds")


pdf('harmony.top5_maker_dotplot.unique.pdf',width = 48, height = 12)
DotPlot(object = adata, features = unique(top5$gene),cols = c("#C0C0C0", "#E41C12")) + RotatedAxis()
dev.off()


write.csv(top10,file = "top10_markergene_list.csv")


sink("stat.txt")
for(i in 1:matrix_num){
    print(paste0('data',i,'.cells = ',length(Cells(o.list[[i]]))))
}
print(paste0('harmony.data','.cells = ',length(Cells(adata))))

for(i in 1:matrix_num){
    print(paste0('data',i,'.genes = ',length(rownames(o.list[[i]]))))
}
print(paste0('harmony.data','.genes = ',length(rownames(adata))))
sink()
#pdf('harmony.top5.heatmap.pdf',width = 36,height = 8)
#DoHeatmap(adata, features = top5$gene) + scale_fill_gradient2(low = c("gray"),high = c('red'))
#dev.off()



#arth_genenames <- read.csv('/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/scRNA_seq/reference_data/Arabidopsis_th/TAIR10.gene2name_noredundant.csv', head=F, stringsAsFactors=F)
#markergenes_id <- read.csv("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/scRNA_seq/Find.cluster/markergenes/markergenelist1.csv", head=T, stringsAsFactor=FALSE)
#setlist1 <- c()
#setlist1 <- markergenes_id[,4]
#for( i in 1:length(setlist1) ){
# if( setlist1[i] %in% arth_genenames[,1] ){
#  setlist1[i] <- arth_genenames[arth_genenames[,1] == setlist1[i], 2 ]
# }
#}

pdf("10_cell_number.pdf", width = 12, height = 8)
dat<-as.data.frame(table(adata@active.ident))
p<-ggplot(dat,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position = "dodge") + geom_text(aes(label=Freq,vjust=-0.5) )
p
dev.off()

write.table(dat,file = "10_cell_number.txt",sep="\t",quote=F,row.names = F)


umap <- data.frame(cellid=adata@assays[["SCT"]]@data@Dimnames[[2]],adata@reductions$umap@cell.embeddings,clusters=adata$seurat_clusters,nFeature_RNA=adata@meta.data[["nFeature_RNA"]],nCount_RNA=adata@meta.data[["nCount_RNA"]])
pdf("dotplot_nfeature.pdf")
p <- ggplot(umap, mapping=aes(x=UMAP_1, y=UMAP_2, colour = nFeature_RNA))+geom_point(size = 0.05)
p + theme_cowplot() +scale_color_gradient(low = "cyan",high = "red")
dev.off()

pdf("dotplot_ncount.pdf")
p <- ggplot(umap, mapping=aes(x=UMAP_1, y=UMAP_2, colour = nCount_RNA))+geom_point(size = 0.05)
p + theme_cowplot() +scale_color_gradient(low = "cyan",high = "red")
dev.off()



