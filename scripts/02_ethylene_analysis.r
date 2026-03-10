#### Cluster level ethylene analysis

setwd("")

library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(reshape2)
library(dplyr)
library(clusterProfiler)
library(igraph)
library(AnnotationHub)
library(org.At.tair.db)
library(org.Athaliana.eg.db)
library(RColorBrewer) 
library(openxlsx)

###### Identify all ethylene regulated markers within clusters


seedling = readRDS("/ceph/ethylene/analysis_231120/results/harmony/harmony.rds")

seedling$celltype.stim = NULL

seedling$celltype.stim <- paste(seedling@meta.data$seurat_clusters, seedling$stim, sep = "_")
seedling$celltype <- Idents(object = seedling)
Idents(object = seedling) <- "celltype.stim"


#### Identify ethylene up and down regulated genes within clusters and filtered for pval < 0.05
for(i in seq(1:length(levels(seedling@meta.data$seurat_clusters))))
{
	et.cluster = paste0(i-1, "_ET")
	ctrl.cluster = paste0(i-1, "_CTRL")
	
	et.response <- FindMarkers(object = seedling, assay="RNA", ident.1 = et.cluster, ident.2 = ctrl.cluster, verbose = T)
	
	assign(paste("et.response.", i-1, sep = ""), et.response)
	assign(paste("et.up.", i-1, sep = ""), row.names(et.response[et.response$avg_log2FC > 0 & et.response$p_val < .05,]))
	assign(paste("et.down.", i-1, sep = ""), row.names(et.response[et.response$avg_log2FC < 0 & et.response$p_val < .05,]))
}

df.up = data.frame()
df.down = data.frame()

for(i in seq(1:length(levels(seedling@meta.data$seurat_clusters))))
{
	df.up = rbind(df.up, length(get(paste("et.up.", i-1, sep = ""))))
	df.down = rbind(df.down, length(get(paste("et.down.", i-1, sep = ""))))
}

et.cluster = cbind(df.up, -df.down)
colnames(et.cluster) = c("ET.up", "ET.down")
et.cluster$cluster <- seq.int(nrow(et.cluster))
et.cluster$cluster = et.cluster$cluster - 1

df.m = melt(et.cluster, id.vars = "cluster")

pdf("results/harmony/ET_pval_clusters.pdf")
p = ggplot(df.m, aes(x = cluster, y = value, fill = variable, width=.8)) 
p + geom_bar(stat="identity", position="identity") + 
	scale_x_continuous(breaks=0:(nrow(et.cluster)-1)) +
	theme_classic()
dev.off()


######## Overlap of ethylene upregulated genes

up.markers.list = NULL
up.markers.list = as.list(up.markers.list)

down.markers.list = NULL
down.markers.list = as.list(down.markers.list)

for(i in seq(1:length(levels(seedling@meta.data$seurat_clusters))))
{
	up.genes = paste0("et.up.", i-1)
	down.genes = paste0("et.down.", i-1)
	temp.list = NULL
	
	up.list = list(get(up.genes))
	names(up.list) = up.genes
	
	down.list = list(get(down.genes))
	names(down.list) = down.genes
	
	up.markers.list = c(up.markers.list, up.list)
	down.markers.list = c(down.markers.list, down.list)
}


et.up.markers.list = up.markers.list
et.down.markers.list = down.markers.list

##### raw number of overlapping marker genes
up.mat.raw <- sapply(up.markers.list, function(x) sapply(up.markers.list, function(y) length(intersect(x,y))))
down.mat.raw <- sapply(down.markers.list, function(x) sapply(down.markers.list, function(y) length(intersect(x,y))))

up.df.raw <- melt(up.mat.raw)
down.df.raw <- melt(down.mat.raw)

##### percentage of overlapping marker genes
up.mat <- sapply(up.markers.list, function(x) sapply(up.markers.list, function(y) length(intersect(x, y)) / length(unique(c(x, y))) * 100))
down.mat <- sapply(down.markers.list, function(x) sapply(down.markers.list, function(y) length(intersect(x, y)) / length(unique(c(x, y))) * 100))

up.df <- melt(up.mat)
down.df <- melt(down.mat)



palette.breaks <- seq(0, 100)

color.palette  <- colorRampPalette(c("white", "red"))(length(palette.breaks) - 1)
down.color.palette  <- colorRampPalette(c("white", "blue"))(length(palette.breaks) - 1)


pdf("results/harmony/et_cluster_similarity_percent.pdf")
ggplot(data=up.df, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "red") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=down.df, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "blue") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



pdf("results/harmony/et_cluster_similarity_raw.pdf")
ggplot(data=up.df.raw, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "red") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=down.df.raw, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "blue") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#### write tables of cluster level 
write.xlsx(setNames(as.list(lapply(up.markers.list, data.frame)), names(up.markers.list)), file="results/harmony/et_up_clusters.xlsx")
write.xlsx(setNames(as.list(lapply(down.markers.list, data.frame)), names(down.markers.list)), file="results/harmony/et_down_clusters.xlsx")


#########
### GO terms of cluster specific ethylene upregulated genes
#########

ck <- compareCluster(geneCluster = up.markers.list, fun = enrichGO, OrgDb = org.Athaliana.eg.db, keyType = 'GID', ont = "BP")
ck <- clusterProfiler::simplify(ck, cutoff=0.7, by="p.adjust", select_fun=min)


pdf("results/harmony/GO/ET/ET_UP_GO_clusters_compare.pdf", height=12)
dotplot(ck, font.size = 6)
dev.off()




##### Ethylene response in ein2-5 mutant cells 


ein2 = subset(seedling, geno == "ein2")

ein2$celltype.stim = NULL
ein2$celltype.stim <- paste(ein2@meta.data$seurat_clusters, seedling$stim, sep = "_")
ein2$celltype <- Idents(object = seedling)
Idents(object = ein2) <- "celltype.stim"




for(i in seq(1:length(levels(ein2@meta.data$seurat_clusters))))
{
	et.cluster = paste0(i-1, "_ET")
	#ctrl.cluster = paste0(i-1, "_AIR")
	ctrl.cluster = paste0(i-1, "_CTRL")
	
	et.response <- FindMarkers(object = ein2, assay="RNA", ident.1 = et.cluster, ident.2 = ctrl.cluster, verbose = T)
	
	assign(paste("et.response.", i-1, sep = ""), et.response)
	
	assign(paste("et.up.", i-1, sep = ""), row.names(et.response[et.response$avg_log2FC > 0 & et.response$p_val_adj < .01,]))
	assign(paste("et.down.", i-1, sep = ""), row.names(et.response[et.response$avg_log2FC < 0 & et.response$p_val_adj < .01,]))
}


df.up = data.frame()
df.down = data.frame()

for(i in seq(1:length(levels(ein2@meta.data$seurat_clusters))))
{
	df.up = rbind(df.up, length(get(paste("et.up.", i-1, sep = ""))))
	df.down = rbind(df.down, length(get(paste("et.down.", i-1, sep = ""))))
}

et.cluster = cbind(df.up, -df.down)
colnames(et.cluster) = c("ein2.ET.up", "ein2.ET.down")
et.cluster$cluster <- seq.int(nrow(et.cluster))
et.cluster$cluster = et.cluster$cluster - 1


df.m = melt(et.cluster, id.vars = "cluster")

pdf("results/harmony/ein2_ET_fdr_clusters.pdf")
p = ggplot(df.m, aes(x = cluster, y = value, fill = variable, width=.8)) 
p + geom_bar(stat="identity", position="identity") + 
	scale_x_continuous(breaks=0:(nrow(et.cluster)-1)) +
	ylim(-375, 375) +
	theme_classic()
dev.off()


######## Overlap of ethylene upregulated genes

up.markers.list = NULL
up.markers.list = as.list(up.markers.list)


down.markers.list = NULL
down.markers.list = as.list(down.markers.list)

for(i in seq(1:length(levels(seedling@meta.data$seurat_clusters))))
{
	up.genes = paste0("et.up.", i-1)
	down.genes = paste0("et.down.", i-1)
	
	temp.list = NULL
	
	up.list = list(get(up.genes))
	names(up.list) = up.genes
	
	down.list = list(get(down.genes))
	names(down.list) = down.genes
	
	up.markers.list = c(up.markers.list, up.list)
	down.markers.list = c(down.markers.list, down.list)
}

et.up.markers.list = up.markers.list
et.down.markers.list = down.markers.list


##### raw number of overlapping marker genes
up.mat.raw <- sapply(up.markers.list, function(x) sapply(up.markers.list, function(y) length(intersect(x,y))))
down.mat.raw <- sapply(down.markers.list, function(x) sapply(down.markers.list, function(y) length(intersect(x,y))))


##### percentage of overlapping marker genes
up.mat <- sapply(up.markers.list, function(x) sapply(up.markers.list, function(y) length(intersect(x, y)) / length(unique(c(x, y))) * 100))
down.mat <- sapply(down.markers.list, function(x) sapply(down.markers.list, function(y) length(intersect(x, y)) / length(unique(c(x, y))) * 100))

up.df.raw <- melt(up.mat.raw)
down.df.raw <- melt(down.mat.raw)

up.df <- melt(up.mat)
down.df <- melt(down.mat)


palette.breaks <- seq(0, 100)


color.palette  <- colorRampPalette(c("white", "red"))(length(palette.breaks) - 1)
down.color.palette  <- colorRampPalette(c("white", "blue"))(length(palette.breaks) - 1)


pdf("results/harmony/ein2_et_cluster_similarity_percent.pdf")
ggplot(data=up.df, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "red") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=down.df, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "blue") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("results/harmony/ein2_et_cluster_similarity_raw.pdf")
ggplot(data=up.df.raw, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "red") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=down.df.raw, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value), color="grey") +
scale_fill_gradient(low = "white", high = "blue") +
labs(x = "cluster number", y = "cluster number") +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




library(openxlsx)
write.xlsx(setNames(as.list(lapply(up.markers.list, data.frame)), names(up.markers.list)), file="results/harmony/ein2_et_up_clusters.xlsx")
write.xlsx(setNames(as.list(lapply(down.markers.list, data.frame)), names(down.markers.list)), file="results/harmony/ein2_et_down_clusters.xlsx")

