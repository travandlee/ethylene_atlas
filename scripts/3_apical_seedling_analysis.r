#### Apical hook subset analysis

setwd("")
system("mkdir results/harmony/apical_subcluster")
system("mkdir results/harmony/apical_subcluster/L2_subcluster")

library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer) 
library(gplots)
library(RColorBrewer)


seedling = readRDS("results/harmony/harmony.rds")

subset.clusters = c(0,1,2,3,5,7,8,11,20)

seedling$seedling_clusters = seedling$seurat_clusters

apical_seedling <- subset(seedling, idents = subset.clusters)

apical_seedling = FindVariableFeatures(apical_seedling)
apical_seedling = RunPCA(apical_seedling, ndims=50)
apical_seedling <- RunUMAP(apical_seedling, dims = 1:50, verbose = TRUE)
apical_seedling <- FindNeighbors(apical_seedling, dims = 1:50, verbose = FALSE)
apical_seedling <- FindClusters(apical_seedling, verbose = FALSE, resolution = 1)

apical_seedling$apical_cluster = apical_seedling$seurat_clusters

huntrix = c("#423b9c", "#6ea0b4", "#aa82c8", "#8bd282", "#b4fac8", "#91beaa", "#8ca0a0", "#c8415a", "#d29f41", "#FFD700")
huntrix_pal <- colorRampPalette(huntrix)(n = length(levels(apical_seedling$seurat_clusters)))

pdf("apical_cluster_UMAP_recolor.pdf")
DimPlot(object = apical_seedling, reduction = "umap", label = TRUE, cols = huntrix_pal)
dev.off()


saveRDS(apical_seedling, "results/harmony/apical_subcluster/apical_clusters.rds")


###### Identify all cluster markers
cluster.markers <- FindAllMarkers(object = apical_seedling, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
cluster.markers %>% filter(avg_log2FC > 1)

markers.df = cluster.markers

desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
#desc[,1] <- gsub("\\..*", "", desc[,1]); desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
markers.df <- cbind(Desc=descv[gsub("_.*", "", markers.df$gene)], markers.df)

write.table(markers.df, "results/harmony/apical_subcluster/apical_subcluster_markers.xls", sep="\t", row.names=T, col.names=NA)


### top marker
top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


pdf("results/harmony/apical_subcluster/apical_cluster_markers.pdf")
DotPlot(apical_seedling, features = unique(top1$gene), dot.min=.03) + RotatedAxis()
dev.off()

#### expression of spatial markers
GH9B8 = "AT2G32990"
LACS2 = "AT1G49430"
EIR1 = "AT5G57090"
XTH7 = "AT4G37800"
AT4G13210 = "AT4G13210"

genes = c(GH9B8, LACS2, EIR1, XTH7, AT4G13210)

pdf("fig1F_apical_seedling_spatial_epidermis_markers.pdf")
DotPlot(hook_10x, features=genes)
dev.off()



PEL3 = c("AT5G23940")
ANT = c("AT4G37750")
AIL6 = "AT5G10510"
UGT85A1 = "AT1G22400"

pdf("figS4D_apical_seedling_spatial_markers.pdf")
FeaturePlot(hook_10x, features=c(PEL3, ANT), order=T)
DotPlot(hook_10x, features=c(PEL3, ANT))
DimPlot(hook_10x, label=T)
FeaturePlot(hook_10x, features=c(PEL3, ANT, AIL6, UGT85A1), order=T)
DotPlot(hook_10x, features=c(PEL3, ANT, AIL6, UGT85A1))

dev.off()


##### Iterative subclustering

for(i in seq(1:length(levels(apical_seedling@meta.data$seurat_clusters))))
{
	file_path = "results/harmony/apical_subcluster/L2_subcluster/"
	title = paste0("subcluster_cluster_", i-1, ".pdf")
	file_name = paste0(file_path, title)
	table_title = paste0("markers_apical_subcluster_cluster_", i-1, ".xls")
	table.name = paste0(file_path, table_title)
	
	L2.hook = NULL
	
	subset.clusters = c(i-1)
	
	L2.hook <- subset(apical_seedling, idents = subset.clusters)
	
	L2.hook = try(RunPCA(L2.hook, ndims=50) %>%
		RunUMAP(dims = 1:30, verbose = FALSE) %>%
		FindNeighbors(dims = 1:30, verbose = FALSE) %>%
		FindClusters(verbose = FALSE, resolution = .4))
	
	print(file_name)
	
	pdf(file_name)
	
	try(print(ElbowPlot(object = L2.hook, ndims=50)))
	p0 = try(DimPlot(object = L2.hook, reduction = "umap", label = TRUE))
	p1 <- try(DimPlot(object = L2.hook, reduction = "umap", group.by = "stim"))
	p2 <- try(DimPlot(object = L2.hook, reduction = "umap", label = TRUE))
	p3 <- try(DimPlot(object = L2.hook, reduction = "umap", group.by = "cond"))
	p4 <- try(DimPlot(object = L2.hook, reduction = "umap", group.by = "geno"))
	
	try(print(plot_grid(p0)))
	try(print(plot_grid(p1)))
	try(print(plot_grid(p2)))
	try(print(plot_grid(p3)))
	try(print(plot_grid(p4)))
	try(print(FeaturePlot(object = L2.hook, features = "nCount_RNA", min.cutoff = "q10", max.cutoff = "q90", order=T)))
	try(print(FeaturePlot(object = L2.hook, features = "nFeature_RNA", min.cutoff = "q10", max.cutoff = "q90", order=T)))
	
	cluster.markers <- try(FindAllMarkers(object = L2.hook, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25))
	
	markers.df = cluster.markers
	
	desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
	descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
	try(markers.df <- tryCatch(cbind(Desc=descv[gsub("_.*", "", markers.df$gene)], markers.df)))
	try(write.table(markers.df, paste0(dir, table_title), sep="\t", row.names=T, col.names=NA))
	
	top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
	
	try(print(DotPlot(L2.hook, features = rev(top1$gene), dot.min=.03) + RotatedAxis()))	
	try(print(file_name))
	try(dev.off())

	saveRDS(L2.hook, paste0(paste0(file_path, "subcluster_cluster_", i-1, ".rds"))

}


### expression of spatial markers in apical seedling subcluster

c11 = readRDS(paste0(file_path, "subcluster_cluster_11.rds"))

zoey = c("#7dc7c3", "#d7a445", "#79aefc", "#be8cd2", "#191964", "#e69cb5")

pdf("figS4E_subcluster_spatial_expression.pdf")
DimPlot(c11, label=T, cols=zoey)
DotPlot(c11, features=AIL6)
DotPlot(c11, features=UGT85A1)
dev.off()
