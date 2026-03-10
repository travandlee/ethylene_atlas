#### spatial multi-modal integration

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Signac)
library(GenomeInfoDb)
library(GenomicFeatures)
library(rtracklayer)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BiocParallel)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(gprofiler2)
library(reshape2)
library(ComplexUpset)
library(cowplot)
library(magrittr)
library(harmony)
library(magrittr)
library(RSQLite)
library(openxlsx)

register(SerialParam())
options(future.globals.maxSize = 3e+10)

system("mkdir /ceph/ethylene/multiome/spatial_integration_250825")
setwd("/ceph/ethylene/multiome/spatial_integration_250825")

huntrix = c("#252157", "#6ea0b4", "#aa82c8", "#8bd282", "#b4fac8", "#91beaa", "#8ca0a0", "#c8415a", "#d29f41", "#fadc82")
huntrix_pal = data.frame(row.names = c(0:36), cols = colorRampPalette(huntrix)(n = 37))

rumi = c("#ca6dca", "#5826a3", "#ffa22b", "#b6a364", "#f50501")


gene_list = read.delim("/ceph/ethylene/merscope/FinalGeneList_DM1695_Arabidopsis_thaliana.csv", row.names=1, header=TRUE, sep = ",")
blacklist = import("/gale/netapp/home/trlee/general/Tn5_rRNA_blacklist.bed")


###### updated subset from https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R
subset_opt <- function(
    object = NULL, 
    subset = NULL, 
    cells = NULL, 
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    ...)
{
  
  if (Update.slots) { 
    message("Updating object slots..")
    object %<>% UpdateSlots()
  }
  
  message("Cloing object..")
  obj_subset <- object
  
  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) { 
    cells <- Cells(obj_subset)[cells]
  }
  
  if (!missing(subset) || !is.null(idents)) {
    message("Extracting cells matched to `subset` and/or `idents`")
  }
  
  if (class(obj_subset) == "FOV") {
    message("object class is `FOV` ")
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }
  
  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    message("Matching cells for object class `FOV`..")
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else { 
    # check if cells are present in all FOV
    message("Matching cells in FOVs..")
    cells_check <-
      lapply(Images(obj_subset) %>% seq, 
             function(i) { 
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells) 
             }) %>% unlist
  }
  
  if (all(cells_check)) { 
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..")
    obj_subset %<>% base::subset(cells = cells, 
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    message("Subsetting FOVs..")
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells, 
                       idents = idents, 
                       features = features, 
                       ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }
    
  } else { 
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
          message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
          message("Subsetting Centroids..")
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells, 
                       idents = idents, 
                       features = features, 
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("Removing FOVs where cells are NOT found: ", "\n", 
            paste0(Images(object)[which(!cells_check == TRUE)], "\n"))
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] } 
    
    # subset final object
    message("..subset final object")
    obj_subset %<>% 
      base::subset(cells = cells,
                   idents = idents,
                   features = features, 
                   ...)
  }
  
  if (Update.object && !class(obj_subset) == "FOV") { 
    message("Updating object..")
    obj_subset %<>% UpdateSeuratObject() }
  
  message("Object is ready!")
  return(obj_subset)
  
}




#### load multiome and spatial + snRNA-seq datsaets
multiome = readRDS("/ceph/ethylene/multiome/integration_250805/multiome_integrated_timecourse_250825.rds")
spatial = readRDS("/ceph/ethylene/merscope/analysis_250205/spatial_droplet_integrated_250205.rds")

DefaultAssay(spatial) = "RNA"
DefaultAssay(multiome) = "RNA"

spatial[["SCT"]] = NULL
multiome[["SCT"]] = NULL

multiome$assay = "multiome"
multiome$geno = "Col"

spatial = JoinLayers(spatial)
multiome = JoinLayers(multiome)

merged = merge(spatial, multiome)


saveRDS(merged, "spatial_multiome_merged_250825.rds")


DefaultAssay(merged) = "RNA"


spatial_genes = Features(spatial[["Vizgen"]])

merged = JoinLayers(merged, layers = "counts")

merged[["RNA"]] <- split(merged[["RNA"]], f = merged$assay, layers = "counts")

#### integrate based on expression of spatial genes
merged <- SCTransform(merged, residual.features=spatial_genes, return.only.var.genes=TRUE, clip.range = c(-10, 10))

merged <- RunPCA(merged, assay="SCT", ndims=50, features = spatial_genes)

merged <- IntegrateLayers(
  object = merged, method = HarmonyIntegration,
  orig.reduction = "pca", normalization.method = "SCT", new.reduction = "harmony", features=spatial_genes,
  verbose = TRUE
)

merged <- merged %>% 
    RunUMAP(reduction = "harmony", dims = 1:50) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
    FindClusters(resolution = 0.2) %>% 
    identity()

merged = FindClusters(merged, resolution = 0.25) %>% 
    identity()


pdf("merfish_multiome_integration_250825.pdf")
DimPlot(merged, reduction = "umap", label=T, raster=FALSE, cols = huntrix_pal$cols) + theme(rect = element_rect(fill = "transparent"))
DimPlot(merged, split.by = "orig.ident", raster=FALSE)
DimPlot(merged, group.by = "orig.ident", raster=FALSE)
DimPlot(merged, group.by = "assay", raster=FALSE)
DimPlot(merged, split.by = "assay", raster=FALSE)
DimPlot(merged, split.by = "geno", raster=FALSE)
dev.off()

pdf("merfish_multiome_integration_250825_split_assay.pdf", width=16)
DimPlot(merged, split.by = "assay", raster=FALSE)
dev.off()



saveRDS(merged, "spatial_multiome_integrated_250825.rds")



################
## DE analysis
################


merged = PrepSCTFindMarkers(merged)
cluster.markers <- FindAllMarkers(object = merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
#cluster.markers <- FindAllMarkers(object = merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)



markers.df = cluster.markers

desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
#desc[,1] <- gsub("\\..*", "", desc[,1]); desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
markers.df <- cbind(Desc=descv[gsub("_.*", "", markers.df$gene)], markers.df)

write.table(markers.df, "spatial_multiome_integrated_markers_250825.xls", sep="\t", row.names=T, col.names=NA)


### top markers
top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)




pdf("spatial_multiome_integration_dotplot_top_marker.pdf")
DotPlot(merged, features = unique(top1$gene), dot.min=.03) + RotatedAxis()
DotPlot(merged, features = unique(top10$gene), dot.min=.03) + RotatedAxis()
dev.off()



############
#### Spatial and droplet markers
############

setwd("/ceph/ethylene/multiome/spatial_integration_250825")

merged = readRDS("spatial_multiome_integrated_250825.rds")


########
### MERFISH markers only
########


droplet = subset_opt(merged, assay==c("10x"))
multi = subset_opt(merged, assay==c("multiome"))
merfish = subset(merged, assay=="merfish")
merfish.cluster.markers <- FindAllMarkers(object = merfish, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, features = gene_list$geneId)
merfish.cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


merfish.markers.df = merfish.cluster.markers
merfish.markers.df$gene_name = gene_list$Vizgen.Gene[match(merfish.markers.df$gene, gene_list$geneId)]


desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
merfish.markers.df <- cbind(Desc=descv[gsub("_.*", "", merfish.markers.df$gene)], merfish.markers.df)

write.table(merfish.markers.df, "spatial_multiome_integration_markers_spatialonly_250825.xls", sep="\t", row.names=T, col.names=NA)


### top marker
merfish.top1 <- merfish.cluster.markers %>% group_by(cluster) %>% top_n(n = -1, wt = p_val_adj) %>% slice_max(n = 1, order_by=p_val_adj, with_ties = FALSE)
merfish.top10 <- merfish.cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% slice_max(n = 10, order_by=p_val_adj, with_ties = FALSE)



pdf("spatial_multiome_integration_spatialonly_dotplot_top_marker.pdf")
DotPlot(merged, features = unique(merfish.top1$gene), dot.min=.03)
DotPlot(merged, features = unique(merfish.top10$gene), dot.min=.03)

DotPlot(merfish, features = unique(merfish.top1$gene), dot.min=.03)
DotPlot(merfish, features = unique(merfish.top10$gene), dot.min=.03)

DotPlot(droplet, features = unique(merfish.top1$gene), dot.min=.03)
DotPlot(droplet, features = unique(merfish.top10$gene), dot.min=.03)

DotPlot(multi, features = unique(merfish.top1$gene), dot.min=.03)
DotPlot(multi, features = unique(merfish.top10$gene), dot.min=.03)

dev.off()




########
### markers from droplet data only
########


droplet.cluster.markers <- FindAllMarkers(object = droplet, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
droplet.cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


droplet.markers.df = droplet.cluster.markers


desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
droplet.markers.df <- cbind(Desc=descv[gsub("_.*", "", droplet.markers.df$gene)], droplet.markers.df)

write.table(droplet.markers.df, "spatial_multiome_integration_markers_dropletonly_250825.xls", sep="\t", row.names=T, col.names=NA)


### top marker
droplet.top1 <- droplet.cluster.markers %>% group_by(cluster) %>% top_n(n = -1, wt = p_val_adj) %>% slice_max(n = 1, order_by=p_val_adj, with_ties = FALSE)
droplet.top10 <- droplet.cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% slice_max(n = 10, order_by=p_val_adj, with_ties = FALSE)


pdf("spatial_multiome_integration_dropletonly_dotplot_top_marker.pdf")
DotPlot(merged, features = unique(droplet.top1$gene), dot.min=.03)
DotPlot(merged, features = unique(droplet.top10$gene), dot.min=.03)

DotPlot(merfish, features = unique(droplet.top1$gene), dot.min=.03)
DotPlot(merfish, features = unique(droplet.top10$gene), dot.min=.03)

DotPlot(droplet, features = unique(droplet.top1$gene), dot.min=.03)
DotPlot(droplet, features = unique(droplet.top10$gene), dot.min=.03)

DotPlot(multi, features = unique(droplet.top1$gene), dot.min=.03)
DotPlot(multi, features = unique(droplet.top10$gene), dot.min=.03)

dev.off()


##########
#### Call cluster specific peaks from integrated dataset
##########

DefaultAssay(merged) = "ATAC"

multiome_clusters = table(subset_opt(merged, assay="multiome")$seurat_clusters)

peaks_cluster <- CallPeaks(merged, group.by = "seurat_clusters", idents=c(0:32), effective.genome.size = 119481543, assay = "ATAC")
# peaks_cluster <- CallPeaks(merged, group.by = "seurat_clusters", idents=c(0:36), effective.genome.size = 119481543, assay = "ATAC")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
# peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks_cluster <- dropSeqlevels(peaks_cluster, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster <- subsetByOverlaps(x = peaks_cluster, ranges = blacklist, invert = TRUE)

# quantify counts in each peak
DefaultAssay(merged) = "ATAC"

macs2_counts <- FeatureMatrix(
  fragments = Fragments(merged),
  features = peaks_cluster,
  cells = colnames(merged)
)


# create a new assay using the MACS2 peak set and add it to the Seurat object
merged[["peaks_cluster"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(merged),
  annotation = annotations
)

DefaultAssay(merged) <- "peaks_cluster"


#########
# add motif information to peak sets
#########


DefaultAssay(merged) = "peaks"

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'plants', species = "3702", all_versions = FALSE)
)

merged <- AddMotifs(
  object = merged,
  genome = BSgenome.Athaliana.TAIR.TAIR9,
  pfm = pfm
)

DefaultAssay(merged) = "peaks_cluster"

merged <- AddMotifs(
  object = merged,
  genome = BSgenome.Athaliana.TAIR.TAIR9,
  pfm = pfm
)

saveRDS(merged, "spatial_multiome_integrated_250825.rds")



##### find cluster specific peaks with differential accessibility between clusters
DefaultAssay(merged) = "peaks_cluster"

da_peaks <- FindAllMarkers(
  object = merged,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks',
  assay = "peaks_cluster"
)

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])


write.table(da_peaks, "spatial_multiome_integration_cluster_diffaccess_markers_250825.xls", sep="\t", row.names=T, col.names=NA)

da_peaks = read.table( "spatial_multiome_integration_cluster_diffaccess_markers_250825.xls", sep="\t", row.names=1, header=T)


#### plot top motifs for differentially accessible cluster specific peaks

pdf("cluster_da_peaks_motif.pdf", width = 10, height = 6)
for(i in seq(1:length(levels(merged$seurat_clusters))))
{
	i=i-1
	
	try(print(MotifPlot(
	  object = merged,
	  assay = "peaks_cluster",
	  motifs = head(rownames(get(paste0("cluster.motifs.", i))))
		)
		+ ggtitle(paste0("cluster_", i, "_differntially_accesible_peak_motifs")))
		)
}
dev.off()





#######
#### Use Chromvar to calculate motif activities within accessibile chromatin
#######

pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 3702, all_versions = FALSE))

DefaultAssay(merged) = "peaks"
motif.matrix <- CreateMotifMatrix(features = granges(merged), pwm = pwm_set, genome = 'BSgenome.Athaliana.TAIR.TAIR9', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
merged <- SetAssayData(merged, assay = 'peaks', slot = 'motifs', new.data = motif.object)


DefaultAssay(merged) = "peaks_cluster"
motif.matrix <- CreateMotifMatrix(features = granges(merged), pwm = pwm_set, genome = 'BSgenome.Athaliana.TAIR.TAIR9', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
merged <- SetAssayData(merged, assay = 'peaks_cluster', slot = 'motifs', new.data = motif.object)

saveRDS(merged, "spatial_multiome_integrated_250825.rds")


DefaultAssay(merged) = 'peaks_cluster'
merged <- RunChromVAR(
  object = merged,
  genome = BSgenome.Athaliana.TAIR.TAIR9,
  verbose = TRUE,
  new.assay.name = "chromvar_cluster"
)


saveRDS(merged, "spatial_multiome_integrated_250825.rds")





########
## Cluster motif activities
########

da_peaks = read.delim("spatial_multiome_integration_cluster_diffaccess_markers_annotated_250825.xls", header=TRUE, row.names=1)

DefaultAssay(merged) = "peaks_cluster"



##### identify cluster specific differential motif activity 
for(i in levels(factor(merged$seurat_clusters)))
{

cluster_peaks = da_peaks[da_peaks$p_val < 0.05,]
cluster_peaks = cluster_peaks[cluster_peaks$cluster == i,]$peak_id

try(assign(paste0("cluster.motifs.", i), FindMotifs(
  object = merged,
  features = cluster_peaks,
  background = 20000,
  assay = "peaks_cluster"
)))
}


enriched.motifs.list = list()

for(i in levels(factor(merged$seurat_clusters)))
{
	enriched.motifs = paste0("cluster.motifs.", i)
	temp.list = NULL
	motifs.list = list(get(enriched.motifs))
	names(motifs.list) = enriched.motifs
	enriched.motifs.list = c(enriched.motifs.list, motifs.list)
}


write.xlsx(setNames(as.list(lapply(enriched.motifs.list, data.frame)), names(enriched.motifs.list)), file="spatial_multiome_integration_cluster_diffaccess_enriched_motifs_annotated_250825.xlsx")





####### ARF analysis
DefaultAssay(merged) = "ATAC"
DefaultAssay(subclust) = "ATAC"

multiome_data = subset_opt(subclust, assay=="multiome")
peaks = rownames(subclust[["peaks"]])
peaks_cluster = rownames(merged[["peaks_cluster"]])


### ARF6

arf2.motif <- ConvertMotifID(merged, name = 'ARF2', assay="peaks")

gene = "AT1G30330"

## peaks[grep("Chr1-1069", peaks)]
## peaks[grep("Chr1-1069", peaks_cluster)]

region1 = "Chr1-10690680-10691024"
region2 = "Chr1-10691662-10692099"
region3 = "Chr1-10692803-10693317"
region4 = "Chr1-10695430-10696573"
region1_grange = StringToGRanges(region1)
region2_grange = StringToGRanges(region2)
region3_grange = StringToGRanges(region3)
region4_grange = StringToGRanges(region4)

regions = c(region1, region2, region3, region4)

regions_grange = c(region1_grange, region2_grange, region3_grange, region4_grange)

multiome_data <- LinkPeaks(
  object = multiome_data,
  peak.assay = "peaks",
  expression.assay = "alra",
  genes.use = gene,
  gene.id = TRUE,
  verbose=T,
  gene.coords = annotations
)


pdf("mechanism_ARF6_regulation.pdf")
DotPlot(merged, features=c(paste0("peaks_", regions), paste0("alra_", gene), paste0("sct_", gene), arf2.motif), col.min=0.2) + RotatedAxis()
DotPlot(subclust, features=c(paste0("peaks_", regions), paste0("alra_", gene), paste0("sct_", gene), arf2.motif), col.min=0.2) + RotatedAxis()

CoveragePlot(merged, region = gene, features = gene, region.highlight = c(regions_grange), extend.upstream = 5000, extend.downstream = 5000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE) & scale_fill_manual(values = huntrix_pal$cols)
CoveragePlot(subclust, region = gene, features = gene, region.highlight = c(regions_grange), extend.upstream = 5000, extend.downstream = 5000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE) & scale_fill_manual(values = rumi)
dev.off()



### ARF2

gene = "AT5G62000"

## peaks[grep("Chr5-249", peaks)]

region1 = "Chr5-24906138-24906592"
region2 = "Chr5-24908005-24908398"
region3 = "Chr5-24909409-24909816"
region1_grange = StringToGRanges(region1)
region2_grange = StringToGRanges(region2)
region3_grange = StringToGRanges(region3)

regions = c(region1, region2, region3)

regions_grange = c(region1_grange, region2_grange, region3_grange)

multiome_data <- LinkPeaks(
  object = multiome_data,
  peak.assay = "peaks",
  expression.assay = "alra",
  genes.use = gene,
  gene.id = TRUE,
  verbose=T,
  gene.coords = annotations
)


pdf("mechanism_ARF2_regulation.pdf")
DotPlot(merged, features=c(paste0("peaks_", regions), paste0("alra_", gene), paste0("sct_", gene), arf2.motif), col.min=0.2) + RotatedAxis()
DotPlot(subclust, features=c(paste0("peaks_", regions), paste0("alra_", gene), paste0("sct_", gene), arf2.motif), col.min=0.2) + RotatedAxis()

CoveragePlot(merged, region = gene, features = gene, region.highlight = c(regions_grange), extend.upstream = 5000, extend.downstream = 5000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE) & scale_fill_manual(values = huntrix_pal$cols)
CoveragePlot(subclust, region = gene, features = gene, region.highlight = c(regions_grange), extend.upstream = 5000, extend.downstream = 5000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE) & scale_fill_manual(values = rumi)
dev.off()


VlnPlot(merged, features=gene, pt.size=0, assay="alra")
VlnPlot(subclust, features=gene, pt.size=0, assay="alra")






### GATA2

gene = "AT2G45050"

## peaks[grep("Chr2-185", peaks)]
## peaks_cluster[grep("Chr2-185", peaks_cluster)]

region1 = "Chr2-18579993-18581056"
region2 = "Chr2-18581278-18581953"
region3 = "Chr2-18582689-18582903"
region1_grange = StringToGRanges(region1)
region2_grange = StringToGRanges(region2)
region3_grange = StringToGRanges(region3)

regions = c(region1, region2, region3)

regions_grange = c(region1_grange, region2_grange, region3_grange)

multiome_data <- LinkPeaks(
  object = multiome_data,
  peak.assay = "peaks_cluster",
  expression.assay = "alra",
  genes.use = gene,
  gene.id = TRUE,
  verbose=T,
  gene.coords = annotations
)

DefaultAssay(merged) = "peaks_cluster"
DefaultAssay(subclust) = "peaks_cluster"

pdf("mechanism_GATA2_regulation.pdf")
DotPlot(merged, features=c(regions, paste0("alra_", gene), paste0("sct_", gene), gata10.motif), col.min=0.2) + RotatedAxis()
DotPlot(subclust, features=c(regions, paste0("alra_", gene), paste0("sct_", gene), gata10.motif), col.min=0.2) + RotatedAxis()

CoveragePlot(merged, region = gene, features = gene, region.highlight = c(regions_grange), extend.upstream = 1000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE) & scale_fill_manual(values = huntrix_pal$cols)
CoveragePlot(subclust, region = gene, features = gene, region.highlight = c(regions_grange), extend.upstream = 1000, extend.downstream = 2000, assay = 'ATAC', expression.assay = "alra", peaks = TRUE, annotation = TRUE, group.by = "seurat_clusters", show.bulk=TRUE) & scale_fill_manual(values = rumi)
dev.off()


VlnPlot(merged, features=gene, pt.size=0, assay="alra")
VlnPlot(subclust, features=gene, pt.size=0, assay="alra")



#######












cell_type = "convex"

gene.list = list(c("AT1G72290", "AT5G17430", "AT3G61830"))


merged <- AddModuleScore(object = merged, features = gene.list, name = cell_type, nbin=15)

pdf(width=12)
print(FeaturePlot(object = merged, features = paste0(cell_type, "1")) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))))
print(FeaturePlot(object = merged, features = paste0(cell_type, "1"), min.cutoff = "q5", max.cutoff = "q95") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))))
print(FeaturePlot(object = merged, split.by="stim", features = paste0(cell_type, "1"), min.cutoff = "q5", max.cutoff = "q95", raster=F) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))))
VlnPlot(merged, features = paste0(cell_type, "1"), pt.size=0)
VlnPlot(merged, features = paste0(cell_type, "1"), pt.size=0, split.by="stim")

dev.off()



cell_type = "concave"

gene.list = list(c("AT5G67100", "AT1G57820", "AT1G47210"))


merged <- AddModuleScore(object = merged, features = gene.list, name = cell_type, nbin=15)

pdf(width=12)
print(FeaturePlot(object = merged, features = paste0(cell_type, "1")) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))))
print(FeaturePlot(object = merged, features = paste0(cell_type, "1"), min.cutoff = "q5", max.cutoff = "q95") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))))
print(FeaturePlot(object = merged, split.by="stim", features = paste0(cell_type, "1"), min.cutoff = "q5", max.cutoff = "q95", raster=F) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))))
VlnPlot(merged, features = paste0(cell_type, "1"), pt.size=0)
VlnPlot(merged, features = paste0(cell_type, "1"), pt.size=0, split.by="stim")

dev.off()


pdf("example_genes/spatial_convex_concave_markers_featureplot.pdf")
FeaturePlot(object = merged, features = paste0("convex", "1"), min.cutoff = "q5", max.cutoff = "q95", raster=FALSE, order=T) & scale_colour_gradientn(colours = c("#cccccc", "#00ff15", "#008906"))
FeaturePlot(object = merged, features = paste0("concave", "1"), min.cutoff = "q5", max.cutoff = "q95", raster=FALSE, order=T) & scale_colour_gradientn(colours = c("#cccccc", "#ffaaff", "#ff00ff", "#a300a3"))
dev.off()






##### Fig4 spatial markers


region = "hypocotylColairsect2r0"


DefaultAssay(merged) = "RNA"
roi = subset_opt(merged, cells = Cells(merged[[region]]$centroid))

DefaultBoundary(roi[[region]]) <- "segmentation"

dataset_name = names(merged)[13:25]
polychrome_palette = DiscretePalette(n = length(levels(merged$seurat_clusters)), palette = "polychrome")
names(polychrome_palette) = as.integer(levels(merged$seurat_clusters))


mols = c("AT4G23750", "AT5G27120")
mol_size = .1
mols_col = c("magenta", "yellow")
mol_alpha = 0.8
alpha=0.5


pdf("GA/GA_c26_basal_hook_markers_spatial.pdf")
ImageFeaturePlot(roi, features = mols, fov = region, axes=FALSE, max.cutoff="q95", border.color="white")
ImageDimPlot(roi, molecules = mols, mols.size=mol_size, mols.cols=mols_col, mols.alpha = mol_alpha, fov = region, cols = polychrome_palette, axes=FALSE, alpha=alpha)
ImageDimPlot(roi, molecules = mols, mols.size=mol_size, mols.cols=mols_col, mols.alpha = mol_alpha, fov = region, cols = rep("gray", length(levels(factor(roi$seurat_clusters)))), axes=FALSE, border.size=0.1)

ImageDimPlot(roi, molecules = mols, mols.size=mol_size, mols.cols=mols_col, mols.alpha = mol_alpha, fov = region, axes=FALSE, alpha=alpha)
ImageDimPlot(roi, fov = region, cols = polychrome_palette, axes=FALSE, alpha=alpha)


dev.off()

