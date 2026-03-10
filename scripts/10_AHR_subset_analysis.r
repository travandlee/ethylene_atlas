##########
## GA analysis of c26
##########

system("mkdir GA")

subclust = CreateSeuratObject(merged[["SCT"]], assay="SCT", meta.data = merged@meta.data)
subclust[["RNA"]] = merged[["RNA"]]
subclust[["ATAC"]] = merged[["ATAC"]]
subclust[["peaks"]] = merged[["peaks"]]
subclust[["peaks_cluster"]] = merged[["peaks_cluster"]]
subclust[["chromvar"]] = merged[["chromvar_cluster"]]
subclust[["alra"]] = merged[["alra"]]




subclust = subset_opt(subclust, seurat_clusters=="26")

DefaultAssay(subclust) = "SCT"

spatial_genes = rownames(merged[["Vizgen"]])

subclust <- RunPCA(subclust, assay="SCT", ndims=50, features = spatial_genes)

subclust <- IntegrateLayers(
  object = subclust, method = HarmonyIntegration,
  orig.reduction = "pca", normalization.method = "SCT", new.reduction = "harmony", features=spatial_genes,
  verbose = TRUE
)

subclust <- subclust %>% 
    RunUMAP(reduction = "harmony", dims = 1:10) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
    FindClusters(resolution = 0.2) %>% 
    identity()

subclust = FindClusters(subclust, resolution = 0.3) %>% 
    identity()


pdf("GA/GA_c26_UMAP.pdf")
DimPlot(subclust, reduction = "umap", label=T, raster=FALSE)
DimPlot(subclust, split.by = "orig.ident", raster=FALSE)
DimPlot(subclust, group.by = "orig.ident", raster=FALSE)
DimPlot(subclust, group.by = "assay", raster=FALSE)
DimPlot(subclust, split.by = "assay", raster=FALSE)
DimPlot(subclust, split.by = "geno", raster=FALSE)

dev.off()


subclust$seurat_clusters = gsub("3", "1", subclust$seurat_clusters)
# subclust$seurat_clusters = gsub("5", "1", subclust$seurat_clusters)
Idents(subclust) = "seurat_clusters"


# rumi = c("#ca6dca", "#5826a3", "#ffa22b", "#000000", "#b6a364", "#f50501")
rumi = c("#ca6dca", "#5826a3", "#ffa22b", "#b6a364", "#f50501")


pdf("GA/GA_c26_UMAP_updated.pdf")
DimPlot(subclust, reduction = "umap", label=T, raster=FALSE)
DimPlot(subclust, reduction = "umap", label=T, raster=FALSE, cols = rumi)
DimPlot(subclust, split.by = "orig.ident", raster=FALSE)
DimPlot(subclust, group.by = "orig.ident", raster=FALSE)
DimPlot(subclust, group.by = "assay", raster=FALSE)
DimPlot(subclust, split.by = "assay", raster=FALSE)
DimPlot(subclust, split.by = "geno", raster=FALSE)

dev.off()


subclust = PrepSCTFindMarkers(subclust)
cluster.markers <- FindAllMarkers(object = subclust, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
#cluster.markers <- FindAllMarkers(object = subclust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)



markers.df = cluster.markers

desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
#desc[,1] <- gsub("\\..*", "", desc[,1]); desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
markers.df <- cbind(Desc=descv[gsub("_.*", "", markers.df$gene)], markers.df)

write.table(markers.df, "GA/GA_c26_subcluster_markers.xls", sep="\t", row.names=T, col.names=NA)





### top marker
top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)

// cluster.markers$ratio = cluster.markers$pct.1/cluster.markers$pct.2
// top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = ratio)

// top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = -p_val_adj) %>% top_n(n = 1, wt = avg_logFC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = p_val_adj)





pdf("GA/GA_c26_spatial_multiome_integration_dotplot_top_marker.pdf")
DotPlot(subclust, features = unique(top1$gene), dot.min=.03) + RotatedAxis()
DotPlot(subclust, features = unique(top10$gene), dot.min=.03) + RotatedAxis()
DotPlot(subclust, features = unique(top1$gene), dot.min=.03, assay="alra") + RotatedAxis()
DotPlot(subclust, features = unique(top10$gene), dot.min=.03, assay="alra") + RotatedAxis()

# dev.off()


pdf("GA/GA_c26_spatial_multiome_integration_features_top_marker.pdf")
FeaturePlot(object = subclust, features = top1$gene[1:4], order = T)
FeaturePlot(object = subclust, features = top1$gene[5:8], order = T)
FeaturePlot(object = subclust, features = top1$gene[9:12], order = T)

dev.off()



saveRDS(subclust, "GA/GA_subclust26_multiome_imputed_250825.rds")


subclust = readRDS("GA/GA_subclust26_multiome_imputed_250825.rds")


# subclust = RunALRA(subclust)
DefaultAssay(subclust) = "alra"


subclust$stim_subclust = paste0(subclust$seurat_clusters, "_", subclust$stim)
subclust$cond_subclust = paste0(subclust$seurat_clusters, "_", subclust$cond)



pdf("GA/GA_subclust26_chromatin.pdf")

CoveragePlot(subclust, region = GA_genes[1], features = GA_genes[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_genes[1], features = GA_genes[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_genes[1], features = GA_genes[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_genes[2], features = GA_genes[2], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_genes[2], features = GA_genes[2], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_genes[2], features = GA_genes[2], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_genes[3], features = GA_genes[3], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_genes[3], features = GA_genes[3], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_genes[3], features = GA_genes[3], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_genes[4], features = GA_genes[4], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_genes[4], features = GA_genes[4], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_genes[4], features = GA_genes[4], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_genes[5], features = GA_genes[5], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_genes[5], features = GA_genes[5], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_genes[5], features = GA_genes[5], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_repressors[1], features = GA_repressors[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_repressors[1], features = GA_repressors[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_repressors[1], features = GA_repressors[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_repressors[2], features = GA_repressors[2], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_repressors[2], features = GA_repressors[2], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_repressors[2], features = GA_repressors[2], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_repressors[3], features = GA_repressors[3], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_repressors[3], features = GA_repressors[3], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_repressors[3], features = GA_repressors[3], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_repressors[4], features = GA_repressors[4], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_repressors[4], features = GA_repressors[4], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_repressors[4], features = GA_repressors[4], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

CoveragePlot(subclust, region = GA_repressors[5], features = GA_repressors[5], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(subclust, region = GA_repressors[5], features = GA_repressors[5], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_subclust")
CoveragePlot(subclust, region = GA_repressors[5], features = GA_repressors[5], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_subclust")

dev.off()


pdf("GA/GA_subclust26_quantificaiton.pdf")
DimPlot(subclust, label=T)
DotPlot(subclust, features=c(GA_genes, GA_signaling, GA_repressors), assay="alra") + RotatedAxis()
DotPlot(subclust, features=c(GA_genes, GA_signaling, GA_repressors), assay="SCT") + RotatedAxis()
VlnPlot(subclust, features=c(GA_genes, GA_signaling, GA_repressors), assay="alra", pt.size=0) + RotatedAxis()
VlnPlot(subclust, features=c(GA_genes, GA_signaling, GA_repressors), assay="SCT", pt.size=0) + RotatedAxis()
FeatureScatter(subclust, feature1="alra_AT5G07200", feature2="alra_AT1G30040")

dev.off()


# Idents(subclust) <- factor(subclust$seurat_clusters,levels=c("6","0","3","2","1","7"))
Idents(subclust) <- factor(subclust$seurat_clusters,levels=c("4","0","2","1","5"))

pdf("GA/GA_subclust26_quantificaiton.pdf")

DimPlot(subclust, label=T, cols=rumi)
DotPlot(subclust, features=c(GA_genes, GA_repressors), assay="alra") + RotatedAxis()
DotPlot(subclust, features=c(GA_genes, GA_repressors), assay="SCT") + RotatedAxis()
VlnPlot(subclust, features=c(GA_genes, GA_repressors), assay="alra", pt.size=0, cols=rumi) + RotatedAxis()
VlnPlot(subclust, features=c(GA_genes, GA_repressors), assay="SCT", pt.size=0, cols=rumi) + RotatedAxis()
FeatureScatter(subclust, feature1="alra_AT5G07200", feature2="alra_AT1G30040")

dev.off()



gata10.motif <- ConvertMotifID(merged, name = 'GATA10', assay="peaks")
erf39.motif <- ConvertMotifID(merged, name = 'ERF039', assay="peaks")

pdf("GA/GA_subclust26_chromvar_gata10_erf39.pdf")
FeaturePlot(subclust, features=gata10.motif, min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(subclust, features=erf39.motif, min.cutoff = "q10", max.cutoff = "q90")

VlnPlot(subclust, features=c(gata10.motif, erf39.motif), pt.size=0, assay="chromvar", cols=rumi) + RotatedAxis()
DotPlot(subclust, features=c(gata10.motif, erf39.motif), assay="chromvar") + RotatedAxis()

VlnPlot(merged, features=c(gata10.motif, erf39.motif), pt.size=0, assay="chromvar_cluster", cols=huntrix_pal$cols) + RotatedAxis()
DotPlot(merged, features=c(gata10.motif, erf39.motif), assay="chromvar_cluster") + RotatedAxis()
dev.off()




genes = c("AT1G72290")

CoveragePlot(merged, region = genes[1], features = genes[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")
CoveragePlot(merged, region = genes[1], features = genes[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "stim_clust")
CoveragePlot(merged, region = genes[1], features = genes[1], extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "cond_clust")
dev.off()





pdf("GA/GA_subclust26_DELLA_RGA1.pdf")
DotPlot(subclust, features=c("sct_AT2G01570", "alra_AT2G01570")) + RotatedAxis()
VlnPlot(subclust, features=c("sct_AT2G01570", "alra_AT2G01570")) + RotatedAxis()

FeaturePlot(subclust, features="sct_AT2G01570", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(subclust, features="alra_AT2G01570", min.cutoff = "q10", max.cutoff = "q90")
CoveragePlot(subclust, region = "AT2G01570", features = "AT2G01570", extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "SCT", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")

dev.off()



pdf("GA/GA_subclust26_DELLA_GAI.pdf")
DotPlot(subclust, features=c("sct_AT1G14920", "alra_AT1G14920")) + RotatedAxis()
VlnPlot(subclust, features=c("sct_AT1G14920", "alra_AT1G14920")) + RotatedAxis()

FeaturePlot(subclust, features="sct_AT1G14920", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(subclust, features="alra_AT1G14920", min.cutoff = "q10", max.cutoff = "q90")
CoveragePlot(subclust, region = "AT1G14920", features = "AT1G14920", extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "SCT", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters")

dev.off()





DefaultAssay(subclust) = "ATAC"

c26_sub_da_peaks = FindAllMarkers(
  object = subclust,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks',
  assay = "peaks"
)





for(i in levels(factor(subclust$seurat_clusters)))
{

cluster_peaks = c26_sub_da_peaks[c26_sub_da_peaks$p_val < 0.05,]
cluster_peaks = cluster_peaks[cluster_peaks$cluster == i,]$gene

try(assign(paste0("cluster.motifs.", i), FindMotifs(
  object = subclust,
  features = cluster_peaks,
  background = 20000,
  assay = "peaks"
)))
}




enriched.motifs.list = list()

for(i in levels(factor(subclust$seurat_clusters)))
{

enriched.motifs = paste0("cluster.motifs.", i)

temp.list = NULL


motifs.list = list(get(enriched.motifs))
names(motifs.list) = enriched.motifs


enriched.motifs.list = c(enriched.motifs.list, motifs.list)

}


library(openxlsx)
write.xlsx(setNames(as.list(lapply(enriched.motifs.list, data.frame)), names(enriched.motifs.list)), file="GA/GA_subclust26_subclust_enriched_motifs.xlsx")




pdf("GA/GA_subclust26_da_peaks_motif.pdf", width = 10, height = 6)
for(i in levels(factor(subclust$seurat_clusters)))
{


try(print(MotifPlot(
  object = subclust,
  assay = "peaks",
  motifs = head(rownames(get(paste0("cluster.motifs.", i)))))
+ ggtitle(paste0("cluster_", i, "_differentially_accesible_peak_motifs")))
)
}
dev.off()



gene = "AT2G45050"
DefaultAssay(subclust) = "alra"
p1=FeaturePlot(subclust, gene, order=T, raster = FALSE)
p3=FeaturePlot(subclust, gene, order=T, raster = FALSE, split.by="stim")
DefaultAssay(subclust) = "SCT"
p2=FeaturePlot(subclust, gene, order=T, raster = FALSE)
(p1 | p2) / p3



DotPlot(subclust, features=c("alra_AT2G45050", "sct_AT2G45050"))


DefaultAssay(subclust) = "chromvar"

pdf("GA/GA_subclust26_mutants_chromatin_dotplot_accessibility_expression_motifactivity.pdf")
DotPlot(subclust, features=c("peaks_Chr2-18580112-18580513", "alra_AT2G45050", gata10.motif), col.min=0.2) + RotatedAxis()
DotPlot(subclust, features=c("peaks_Chr2-18580112-18580513", "alra_AT2G45050", "sct_AT2G45050", gata10.motif), col.min=0.2) + RotatedAxis()
FeaturePlot(subclust, gata10.motif)

CoveragePlot(subclust, region = "AT2G45050", features = "AT2G45050", extend.upstream = 2000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE)


dev.off()







c26_cells = colnames(subclust)
spatial_c26 = rownames(merged@meta.data)[rownames(merged@meta.data) %in% c26_cells]


merged@meta.data$c26 = "no"
merged@meta.data[spatial_c26,]$c26 = "c26"

merged@meta.data$c26_sub = "na"
merged@meta.data[spatial_c26,]$c26_sub = subclust$seurat_clusters




DefaultAssay(merged) = "alra"
Idents(merged) = "c26"
Idents(merged) = "seurat_clusters"

fovs_seg = names(merged)[27:42]
fovs_pt = names(merged)[14:26]

pdf("GA/GA_c26_spatial_location_segs.pdf")
for(i in fovs_seg)
{
print(ImageDimPlot(merged, fov = i, group.by = "c26", axes = TRUE, size = 0.7, border.color = "white", cols = c("gray", "magenta"),
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
print(ImageDimPlot(merged, fov = i, group.by = "c26_sub", axes = TRUE, size = 0.7, border.color = "white", cols = c("gray", rumi),
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}
dev.off()


pdf("GA/GA_c26_spatial_location_pts.pdf")
for(i in fovs_pt)
{
print(ImageDimPlot(merged, fov = i, group.by = "c26", axes = TRUE, size = 0.6, border.size = 0.01, border.color = "white", cols = c("gray27", rumi),
    coord.fixed = TRUE, crop=TRUE, alpha = 0.8) + ggtitle(i))
print(ImageDimPlot(merged, fov = i, group.by = "c26_sub", axes = TRUE, size = 0.6, border.size = 0.01, border.color = "white", cols = c("gray27", rumi),
    coord.fixed = TRUE, crop=TRUE, alpha = 0.8) + ggtitle(i))
}
dev.off()



color.palette  <- colorRampPalette(rumi[1:4]))(9)










######
##  GA subcluster module expression in spatial
######

GA_subclust = read.delim("GA/GA_c26_subcluster_markers.xls", header=T, row.names=1)

GA_sub5_markers = GA_subclust[GA_subclust$cluster == 5,]$gene
GA_sub4_markers = GA_subclust[GA_subclust$cluster == 4,]$gene

merged <- AddModuleScore(object = merged, features = list(GA_sub4_markers), name = "GA_sub4_SCT", assay="SCT")
merged <- AddModuleScore(object = merged, features = list(GA_sub5_markers), name = "GA_sub5_SCT", assay="SCT")

merged <- AddModuleScore(object = merged, features = list(GA_sub4_markers), name = "GA_sub4_alra", assay="alra")
merged <- AddModuleScore(object = merged, features = list(GA_sub5_markers), name = "GA_sub5_alra", assay="alra")


pdf("GA/GA_subclust_sub4_module_spatial.pdf")
for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "GA_sub4_alra1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}
dev.off()

pdf("GA/GA_subclust_sub5_module_spatial.pdf")
for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "GA_sub5_alra1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}
dev.off()


min_value_alra = min(c(merged$GA_sub4_alra1, merged$GA_sub5_alra1))
min_value_sct = min(c(merged$GA_sub4_SCT1, merged$GA_sub5_SCT1))

merged$delta_sub4_sub5_alra = (merged$GA_sub4_alra1 + min_value_alra) - (merged$GA_sub5_alra1 + min_value_alra)
merged$delta_sub5_sub4_alra = (merged$GA_sub5_alra1 + min_value_alra) - (merged$GA_sub4_alra1 + min_value_alra)

merged$delta_sub4_sub5_sct = (merged$GA_sub4_SCT1 + min_value_sct) - (merged$GA_sub5_SCT1 + min_value_sct)
merged$delta_sub5_sub4_sct = (merged$GA_sub5_SCT1 + min_value_sct) - (merged$GA_sub4_SCT1 + min_value_sct)

pal = rev(brewer.pal(11, "RdBu"))

pdf("GA/GA_subclust_DELTA_sub5_sub4_module_spatial_alra.pdf")
for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "delta_sub4_sub5_alra", axes = TRUE, size = 0.7, border.color = "white", cols = pal,
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}

for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "delta_sub5_sub4_alra", axes = TRUE, size = 0.7, border.color = "white", cols = pal,
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}
dev.off()



pdf("GA/GA_subclust_DELTA_sub5_sub4_module_spatial_sct.pdf")
for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "delta_sub4_sub5_sct", axes = TRUE, size = 0.7, border.color = "white", cols = pal,
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}

for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "delta_sub5_sub4_sct", axes = TRUE, size = 0.7, border.color = "white", cols = pal,
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}
dev.off()





for(i in fovs_pt)
{
print(ImageFeaturePlot(merged, fov = i, features= "delta_sub4_sub5_alra", axes = TRUE, size = 0.6, border.size = 0.01, border.color = "white", cols = pal,
    coord.fixed = TRUE, crop=TRUE, alpha = 0.8) + ggtitle(i))
}

for(i in fovs_pt)
{
print(ImageFeaturePlot(merged, fov = i, features= "delta_sub4_sub5_sct", axes = TRUE, size = 0.6, border.size = 0.01, border.color = "white", cols = pal,
    coord.fixed = TRUE, crop=TRUE, alpha = 0.8) + ggtitle(i))
}
dev.off()









merged <- AddModuleScore(object = merged, features = list(GA_biosyn), name = "GA_biosynthesis_alra", assay="alra")
merged <- AddModuleScore(object = merged, features = list(GA_biosyn), name = "GA_biosynthesis_sct", assay="SCT")

merged <- AddModuleScore(object = merged, features = list(GA_cat), name = "GA_catabolism_alra", assay="alra")
merged <- AddModuleScore(object = merged, features = list(GA_cat), name = "GA_catabolism_sct", assay="SCT")

merged <- AddModuleScore(object = merged, features = list(GA_sig), name = "GA_signaling_alra", assay="alra")
merged <- AddModuleScore(object = merged, features = list(GA_sig), name = "GA_signaling_sct", assay="SCT")


pdf("GA/GA_subclust26_GAsignaling_module_spatial_alra.pdf")
for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "GA_biosynthesis_alra1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))

print(ImageFeaturePlot(merged, fov = i, features= "GA_catabolism_alra1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))

print(ImageFeaturePlot(merged, fov = i, features= "GA_signaling_alra1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))

}
dev.off()

pdf("GA/GA_subclust26_GAsignaling_module_spatial_sct.pdf")
for(i in fovs_seg)
{
print(ImageFeaturePlot(merged, fov = i, features= "GA_biosynthesis_sct1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))

print(ImageFeaturePlot(merged, fov = i, features= "GA_catabolism_sct1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))

print(ImageFeaturePlot(merged, fov = i, features= "GA_signaling_sct1", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, crop=TRUE) + ggtitle(i))
}
dev.off()


