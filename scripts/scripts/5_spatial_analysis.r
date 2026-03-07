##### integration with apical seedling snRNA-seq data

library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
library(Matrix)
library(RColorBrewer)
library(paletteer)
library(sf) 
library(geojsonsf)
library(ggplot2)
library(stringr)
library(grid)
library(spdep)
library(pals)

setwd("")


hook_10x = readRDS("results/harmony/apical_subcluster/apical_clusters.rds")
hook_10x$assay = "10x"

spatial = readRDS(file="spatial_only_integration_250205.rds")

DefaultAssay(spatial) = "RNA"
spatial = JoinLayers(spatial)
merged = merge(spatial, hook_10x)

DefaultAssay(merged) = "RNA"
merged <- SCTransform(merged, residual.features=row.names(spatial), return.only.var.genes=TRUE, clip.range = c(-10, 10))
merged <- RunPCA(merged, assay="SCT", ndims=50, features = rownames(spatial))
merged = RunHarmony(merged, "assay", plot_convergence = FALSE, verbose=TRUE)


merged <- merged %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

merged = FindClusters(merged, resolution = 0.6) %>% 
    identity()


a = 0.6

pdf("merfish_integration_250205.pdf")
DimPlot(merged, reduction = "umap", label=T)
DimPlot(merged, split.by = "orig.ident")
DimPlot(merged, group.by = "orig.ident")
DimPlot(merged, group.by = "assay")
DimPlot(merged, split.by = "assay")
DimPlot(merged, split.by = "geno")

dev.off()


merged = PrepSCTFindMarkers(merged)



saveRDS(merged, "merfish_droplet_integration_250205")


merged = readRDS("merfish_droplet_integration_250205")



#### add cropped areas to the seurat object

### SAM and leaf primorida
cropped.coords <- Crop(merged[["hls1_air_sect2_r1"]], x = c(4175, 4310), y = c(2660, 2820), coords = "plot")
merged[["samhls1airsect2r1"]] <- cropped.coords
DefaultBoundary(merged[["samhls1airsect2r1"]]) <- "segmentation"

### cotyledon 1
cropped.coords <- Crop(merged[["Col_ET_sect2_r3"]], x = c(0, 1350), y = c(6800, 7300), coords = "plot")
merged[["cotyledonColETsect2r3"]] <- cropped.coords
DefaultBoundary(merged[["cotyledonColETsect2r3"]]) <- "segmentation"

### cotyledon 2
cropped.coords <- Crop(merged[["Col_ET_sect2_r2"]], x = c(3700, 4350), y = c(3150, 3900), coords = "plot")
merged[["cotyledon2ColETsect2r2"]] <- cropped.coords
DefaultBoundary(merged[["cotyledon2ColETsect2r2"]]) <- "segmentation"

### cotyledon 3
cropped.coords <- Crop(merged[["Col_air_sect1_r0"]], x = c(2500, 3300), y = c(5800, 6800), coords = "plot")
merged[["cotyledon3Colairsect1r0"]] <- cropped.coords
DefaultBoundary(merged[["cotyledon3Colairsect1r0"]]) <- "segmentation"

### cotyledon crop 4
cropped.coords <- Crop(merged[["Col_air_sect2_r0"]], x = c(1400, 1800), y = c(4400, 5600), coords = "plot")
merged[["cotyledon4Colairsect2r0"]] <- cropped.coords
DefaultBoundary(merged[["cotyledon4Colairsect2r0"]]) <- "segmentation"

# hypocotyl crop
cropped.coords <- Crop(merged[["Col_ET_sect1_r0"]], x = c(5500, 6000), y = c(8000, 8300), coords = "plot")
merged[["hypocotylColETsect1r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotylColETsect1r0"]]) <- "segmentation"


# hypocotyl 2 crop
cropped.coords <- Crop(merged[["Col_ET_sect1_r0"]], x = c(3200, 3600), y = c(6300, 6500), coords = "plot")
merged[["hypocotyl2ColETsect1r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl2ColETsect1r0"]]) <- "segmentation"


##### hypocotyl cross section
cropped.coords <- Crop(merged[["Col_ET_sect1_r0"]], x = c(4250, 4500), y = c(8300, 8700), coords = "plot")
merged[["hypocotyl3ColETsect1r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl3ColETsect1r0"]]) <- "segmentation"


##### hypocotyl long section
cropped.coords <- Crop(merged[["Col_air_sect2_r0"]], x = c(1200, 2100), y = c(3000, 3700), coords = "plot")
merged[["hypocotylColairsect2r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotylColairsect2r0"]]) <- "segmentation"


##### hypocotyl long section 2
cropped.coords <- Crop(merged[["Col_air_sect2_r0"]], x = c(2500, 3300), y = c(3500, 4700), coords = "plot")
merged[["hypocotyl2Colairsect2r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl2Colairsect2r0"]]) <- "segmentation"

##### hypocotyl long section 3
cropped.coords <- Crop(merged[["Col_air_sect2_r0"]], x = c(4500, 5200), y = c(2800, 3900), coords = "plot")
merged[["hypocotyl3Colairsect2r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl3Colairsect2r0"]]) <- "segmentation"


##### hypocotyl long section 4
cropped.coords <- Crop(merged[["Col_ET_sect1_r0"]], x = c(2500, 3700), y = c(9250, 9800), coords = "plot")
merged[["hypocotyl4ColETsect1r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl4ColETsect1r0"]]) <- "segmentation"


##### hypocotyl long section 5
cropped.coords <- Crop(merged[["Col_air_sect2_r0"]], x = c(1400, 1800), y = c(4300, 5600), coords = "plot")
merged[["hypocotyl5Colairsect2r0"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl5Colairsect2r0"]]) <- "segmentation"

##### hypocotyl long section 6
cropped.coords <- Crop(merged[["Col_air_sect1_r1"]], x = c(4800, 5500), y = c(3100, 4000), coords = "plot")
merged[["hypocotyl6Colairsect1r1"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl6Colairsect1r1"]]) <- "segmentation"

##### hypocotyl long section 7
cropped.coords <- Crop(merged[["Col_air_sect1_r1"]], x = c(5300, 6200), y = c(1000, 1700), coords = "plot")
merged[["hypocotyl7Colairsect1r1"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl7Colairsect1r1"]]) <- "segmentation"


##### hypocotyl long section 8
cropped.coords <- Crop(merged[["Col_ET_sect2_r1"]], x = c(3600, 4400), y = c(1500, 2200), coords = "plot")
merged[["hypocotyl8ColETsect2r1"]] <- cropped.coords
DefaultBoundary(merged[["hypocotyl8ColETsect2r1"]]) <- "segmentation"



###########
#### save merged object
###########

saveRDS(merged, "spatial_droplet_integrated_250205.rds")


###### merge hypocotyl ROIs

hypocotyl_merge = subset(merged, cells=c(Cells(merged[["hypocotylColETsect1r0"]]), Cells(merged[["hypocotyl3ColETsect1r0"]]),
Cells(merged[["hypocotylColairsect2r0"]]),
Cells(merged[["hypocotyl2Colairsect2r0"]]),
Cells(merged[["hypocotyl3Colairsect2r0"]]),
Cells(merged[["hypocotyl4ColETsect1r0"]]),
Cells(merged[["hypocotyl5Colairsect2r0"]]),
Cells(merged[["hypocotyl6Colairsect1r1"]]),
Cells(merged[["hypocotyl7Colairsect1r1"]]),
Cells(merged[["hypocotyl8ColETsect2r1"]])
))


hypocotyl_merge <- SCTransform(hypocotyl_merge, assay = "RNA", clip.range = c(-10, 10))
hypocotyl_merge <- RunPCA(hypocotyl_merge)
hypocotyl_merge <- RunUMAP(hypocotyl_merge, dims = 1:50)
hypocotyl_merge <- FindNeighbors(hypocotyl_merge, reduction = "pca", dims = 1:50)
hypocotyl_merge <- FindClusters(hypocotyl_merge, resolution = 0.8)


num_clust = max(as.numeric(levels(hypocotyl_merge$seurat_clusters)))
polychrome_pal = polychrome(num_clust + 1)

names(polychrome_pal) = c(0:num_clust)


pdf("fig2_B_hypocotyl_merge_spatial.pdf")
DimPlot(hypocotyl_merge, label=T, label.color = "white", cols = polychrome_pal) + DarkTheme()
ImageDimPlot(hypocotyl_merge, fov = "hypocotylColETsect1r0", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, cols = polychrome_pal)
dev.off()


saveRDS(hypocotyl_merge, "hypocotyl_spatialONLY_integration_250205.rds")

######

hypocotyl_merge = readRDS("hypocotyl_spatialONLY_integration_250205.rds")

hypocotyl_merge = PrepSCTFindMarkers(hypocotyl_merge)
cluster.markers <- FindAllMarkers(object = hypocotyl_merge, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
#cluster.markers <- FindAllMarkers(object = hypocotyl_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

markers.df = cluster.markers

desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
#desc[,1] <- gsub("\\..*", "", desc[,1]); desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
markers.df <- cbind(Desc=descv[gsub("_.*", "", markers.df$gene)], markers.df)

write.table(markers.df, "hypocotyl_spatialONLY_integration_markers_250825.xls", sep="\t", row.names=T, col.names=NA)


### top marker
top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = p_val_adj)



####### spatial analysis of convex / concave clusters

### get the hypocotyl subset
hypocotyl_merge = readRDS(file="hypocotyl_spatialONLY_integration_250205.rds")

hypocotyl = subset(hypocotyl_merge, cells=Cells(hypocotyl_merge[["hypocotylColETsect1r0"]]))

DefaultBoundary(hypocotyl[["hypocotylColETsect1r0"]]) <- "centroid"
coords = GetTissueCoordinates(hypocotyl)
coords[c("x","y")] = Rotation(coords[,1:2], 45)

coords$color = "gray"
coords$location = NA

convex_cortex = colnames(subset(hypocotyl, idents=14))
epi = colnames(subset(hypocotyl, idents=2))
endodermis = colnames(subset(hypocotyl, idents=4))
vascular = colnames(subset(hypocotyl, idents=13))


general_cortex_1 = colnames(subset(hypocotyl, idents=c(0)))
general_cortex_2 = colnames(subset(hypocotyl, idents=c(1)))
unannotated = colnames(subset(hypocotyl, idents=c(5,6,9,11,15)))
concave_cells = colnames(subset(hypocotyl, idents=7))



coords_vasc = coords
coords_convex = coords

coords[coords$cell %in% convex_cortex,]$location = "convex_cortex"
coords[coords$cell %in% epi,]$location = "epidermis"
coords[coords$cell %in% endodermis,]$location = "endodermis"
coords[coords$cell %in% vascular,]$location = "vascular"

coords[coords$cell %in% general_cortex_1,]$location = "general_cortex_1"
coords[coords$cell %in% general_cortex_2,]$location = "general_cortex_2"
coords[coords$cell %in% unannotated,]$location = "unannotated"
coords[coords$cell %in% concave_cells,]$location = "concave_cells"


coords_convex[coords_convex$cell %in% convex_cortex,]$location = "convex_cortex"
coords_convex[coords_convex$cell %in% epi,]$location = "epidermis"
coords_vasc[coords_vasc$cell %in% endodermis,]$location = "endodermis"
coords_vasc[coords_vasc$cell %in% vascular,]$location = "vascular"


cols = c("red", "#16ff32", "#f8a19f", "#325a9b", "#e4e1e3", "#5a5156", "#b00068", "white")
names(cols) = c("epidermis", "endodermis", "vascular", "convex_cortex", "general_cortex_1", "general_cortex_2", "concave_cells", "unannotated")

cols_convex = c("#325a9b", "red")
names(cols_convex) = c("convex_cortex", "epidermis")

cols_vasc = c("#16ff32", "#f8a19f")
names(cols_vasc) = c("endodermis", "vascular")

num_clust = max(as.numeric(levels(hypocotyl_merge$seurat_clusters)))
polychrome_pal = polychrome(num_clust + 1)
names(polychrome_pal) = c(0:num_clust)

coords$location = factor(coords$location, levels = rev(c("convex_cortex", "concave_cells", "general_cortex_1", "general_cortex_2", "endodermis", "vascular",  "unannotated", "epidermis")))
coords$location = factor(coords$location, levels = rev(c("epidermis", "convex_cortex", "concave_cells", "endodermis", "vascular", "general_cortex_1", "general_cortex_2", "unannotated")))

pdf("fig2B_spatial_hypocotyl_cells.pdf")
angle=45
plot.new()
print(ImageDimPlot(hypocotyl, fov = "hypocotylColETsect1r0", axes = TRUE, size = 3, border.color = "white",
    coord.fixed = TRUE, cols = polychrome_pal), vp = viewport(angle=angle))
plot.new()
print(ImageDimPlot(hypocotyl, fov = "hypocotylColETsect1r0", axes = TRUE, size = 0.7, border.color = "white",
    coord.fixed = TRUE, cols = polychrome_pal, boundaries = "segmentation"), vp = viewport(angle=angle))

pdf("fig2C_figS6D_E_spatial_hypocotyl_cell_quantificaiton_.pdf")
ggplot(as.data.frame(coords)) + 
	geom_histogram(binwidth = 10, aes(x=y, fill = location))  +
	scale_fill_manual(values = cols)
ggplot(as.data.frame(coords_convex)) + 
	geom_histogram(binwidth = 10, aes(x=y, fill = location))  +
	scale_fill_manual(values = cols_convex)
ggplot(as.data.frame(coords_vasc)) + 
	geom_histogram(binwidth = 10, aes(x=y, fill = location))  +
	scale_fill_manual(values = cols_vasc)
dev.off()





##### convex/concave BBM HDG12 transcript quantification

hypocotyl_merge = readRDS(file="/ceph/ethylene/merscope/analysis_250205/hypocotyl_spatialONLY_integration_250205.rds")

hypocotyl = subset(hypocotyl_merge, cells=Cells(hypocotyl_merge[["hypocotylColairsect2r0"]]))

HDG12_mols = as.data.frame(hypocotyl[["hypocotylColairsect2r0"]][["molecules"]][]$AT1G17920)
BBM_mols = as.data.frame(hypocotyl[["hypocotylColairsect2r0"]][["molecules"]][]$AT5G17430)


# #### check rotation angle
# DefaultBoundary(hypocotyl[["hypocotylColairsect2r0"]]) <- "segmentation"
# 
angle=-115
# print(ImageDimPlot(hypocotyl, fov = "hypocotylColairsect2r0", axes = TRUE, size = 3, border.color = "white",
#     coord.fixed = TRUE, cols = polychrome_palette), vp = viewport(angle=angle))
# dev.off()



HDG12_mols_rotated = as.data.frame(Rotation(HDG12_mols, angle*pi/180))
colnames(HDG12_mols_rotated) = c("x", "y")
HDG12_mols_rotated$gene = "HDG12"
HDG12_mols_rotated$color = "#ff00ff"

BBM_mols_rotated = as.data.frame(Rotation(BBM_mols, angle*pi/180))
colnames(BBM_mols_rotated) = c("x", "y")
BBM_mols_rotated$gene = "BBM"
BBM_mols_rotated$color = "#ffff00"

molecules_location = rbind(HDG12_mols_rotated, BBM_mols_rotated)


cols = c(BBM = "#ffff00", HDG12 = "#ff00ff")

pdf("fig2K_hook_BBM_HDG12_spatial_location.pdf")
p = ggplot(molecules_location, aes(x, y))
p + geom_point(aes(color=gene)) + 
	scale_color_manual(values = cols) +
	theme_dark()

ggplot(molecules_location) + 
	geom_freqpoly(binwidth = 90, aes(x=x, color = gene)) + 
	scale_color_manual(values = cols) +
	theme_dark()
dev.off()




