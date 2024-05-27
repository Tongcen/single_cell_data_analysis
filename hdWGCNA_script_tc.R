# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

getwd()
setwd("/Users/tongcen/R.studio/hdWGCNA_48HAI/48HAI_new/")

seurat_obj  = readRDS("48HAI_new_harmony.all_tissue.rds")
#sce@active.assay <- "RNA"
seurat_obj = FindClusters(seurat_obj, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1))
library('clustree')
P <- clustree(seurat_obj@meta.data, prefix = "SCT_snn_res.")
print(P)
dev.off()
table(seurat_obj$SCT_snn_res.0.05)

#创建misc slot，挑选基因
# gene_select parameter：
##variable: use the genes stored in the Seurat object’s VariableFeatures
##fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by
##custom: use genes that are specified in a custom list

#seurat_obj <- SetupForWGCNA(
  #seurat_obj,
  #gene_select = "fraction", # the gene selection approach
  #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  #wgcna_name = "tutorial" # the name of the hdWGCNA experiment
#)

seurat_obj <- SetupForWGCNA(seurat_obj, 
                            gene_select = "variable", 
                            wgcna_name = "PBMC")

str(seurat_obj@misc)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("SCT_snn_res.0.05"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'SCT_snn_res.0.05',
  min_cells = 50# set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

#seurat_obj@misc$PBMC$wgcna_metacell_obj

#seurat_obj@misc$PBMC$wgcna_metacell_obj@meta.data %>% head

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c(0,1,2,3,4,5),# can be all groups
  group.by='SCT_snn_res.0.05' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)


seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

plot_list <- PlotSoftPowers(seurat_obj)

wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)


seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=,
  setDatExpr=FALSE,
  overwrite_tom = TRUE
)

PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

seurat_obj@misc$PBMC$wgcna_modules %>% head

table(seurat_obj@misc$PBMC$wgcna_modules$module)

seurat_obj <- ModuleEigengenes(
  seurat_obj
)
seurat_obj@misc$PBMC$MEs %>% head

seurat_obj <- ModuleConnectivity(
  seurat_obj
)
seurat_obj@misc$PBMC$wgcna_modules[,c(-1,-2,-3)] %>% head

hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
hub_df

## 计算每个细胞对于每个模块hub基因的表达活性(module score),可使用seurat包或者Ucell包
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

seurat_obj@misc$PBMC$module_scores %>% head

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)
seurat_obj@misc$PBMC$module_scores %>% head

#每个细胞对于每个模块的特征值
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='MEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

wrap_plots(plot_list, ncol=2)

#每个细胞对于每个模块hub基因的表达活性
seurat_obj@misc$PBMC$module_scores %>% head
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
wrap_plots(plot_list, ncol=2)

#每个细胞对于每个模块hub基因的表达活性
MEs <- GetMEs(seurat_obj)
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
modules <- GetModules(seurat_obj)
saveRDS(modules,"modules.RDS")

mods <- colnames(MEs); mods <- mods[mods != 'grey']
modules_genes <- modules[modules$color != 'grey',]
modules_genes <- modules_genes[order(modules_genes$module),]
write.csv(modules_genes,"modules_genes.csv",row.names = F)

p <- DotPlot(seurat_obj, features = mods, group.by = 'SCT_snn_res.0.05')
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p


############################# network,can not save? ########################################
library(igraph)
ModuleNetworkPlot(seurat_obj)
# if no error, file ModuleNetworks has been saved

# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
g <- HubGeneNetworkPlot(seurat_obj)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 , # min distance between points in UMAP space
  wgcna_name = "PBMC"
)
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
mod_umap <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +umap_theme()

png("ModuleUMAPPlot.png", height=3, width=3, units='in', res=200)
mod_umap
dev.off()
# ModuleUMAPPlot(
#   seurat_obj,
#   edge.alpha=0.25,
#   sample_edges=TRUE,
#   edge_prop=0.1, # proportion of edges to sample (20% here)
#   label_hubs=2 ,# how many hub genes to plot per module?
#   keep_grey_edges=FALSE
#   
# )

# g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
# png("ModuleUMAPPlot.png", height=16, width=20, units='in', res=200)
# g
# dev.off()

saveRDS(seurat_obj@meta.data,"meta_data.rds")
saveRDS(seurat_obj, file='hdWGCNA_object.rds')

