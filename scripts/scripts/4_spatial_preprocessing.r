##### spatial preprocessing and integration with apical seedling snRNA-seq data

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

setwd("")

gene_list = read.delim("/ceph/ethylene/merscope/FinalGeneList_DM1695_Arabidopsis_thaliana.csv", row.names=1, header=TRUE, sep = ",")

polygon_files = c(
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section1_region3_segmentation_noZ_250205/baysor_hls1_air_section1_region3_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section1_region2_segmentation_noZ_250205/baysor_hls1_air_section1_region2_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section1_region1_segmentation_noZ_250205/baysor_hls1_air_section1_region1_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section2_region3_segmentation_noZ_250205/baysor_hls1_air_section2_region3_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section2_region2_segmentation_noZ_250205/baysor_hls1_air_section2_region2_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section2_region1_segmentation_noZ_250205/baysor_hls1_air_section2_region1_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section2_region3_segmentation_noZ_250205/baysor_Col_ET_section2_region3_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section2_region2_segmentation_noZ_250205/baysor_Col_ET_section2_region2_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section2_region1_segmentation_noZ_250205/baysor_Col_ET_section2_region1_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section1_region0_segmentation_noZ_250205/baysor_Col_ET_section1_region0_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_air_section1_region0_segmentation_noZ_250205/baysor_Col_air_section1_region0_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_air_section1_region1_segmentation_noZ_250205/baysor_Col_air_section1_region1_segmentation_noZ_250205_polygons_2d.json", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_air_section2_region0_segmentation_noZ_250205/baysor_Col_air_section2_region0_segmentation_noZ_250205_polygons_2d.json"
)

dataset_name = c(
"hls1_air_sect1_r3", 
"hls1_air_sect1_r2", 
"hls1_air_sect1_r1", 
"hls1_air_sect2_r3", 
"hls1_air_sect2_r2", 
"hls1_air_sect2_r1", 
"Col_ET_sect2_r3", 
"Col_ET_sect2_r2", 
"Col_ET_sect2_r1", 
"Col_ET_sect1_r0", 
"Col_air_sect1_r0", 
"Col_air_sect1_r1",
"Col_air_sect2_r0"
)

#### load in polygon files
for(i in seq(length(polygon_files)))
{
assign(dataset_name[i], read_sf(polygon_files[i]))
}

#### remove crs system
st_crs(hls1_air_sect1_r3) = NA
st_crs(hls1_air_sect1_r2) = NA
st_crs(hls1_air_sect1_r1) = NA
st_crs(hls1_air_sect2_r3) = NA
st_crs(hls1_air_sect2_r2) = NA
st_crs(hls1_air_sect2_r1) = NA
st_crs(Col_ET_sect2_r3) = NA
st_crs(Col_ET_sect2_r2) = NA
st_crs(Col_ET_sect2_r1) = NA
st_crs(Col_ET_sect1_r0) = NA
st_crs(Col_air_sect1_r0) = NA
st_crs(Col_air_sect1_r1) = NA
st_crs(Col_air_sect2_r0) = NA


#### polygon coordinates sanity check 
pdf("all_polygons_spatial_datasets.pdf")
for(i in seq(length(dataset_name)))
{
	print(plot(get(dataset_name[i])))
}
dev.off()


##### remove all invalid polygons
valid_polygons = as.character()

for(i in seq(length(dataset_name)))
{
	assign(paste0(dataset_name[i], "_valid"), subset(get(dataset_name[i]), st_is_valid(get(dataset_name[i]))))
	valid_polygons = c(valid_polygons, paste0(dataset_name[i], "_valid"))
}

#### polygon coordinates sanity check 
pdf("valid_polygons_spatial_datasets.pdf")
for(i in seq(length(valid_polygons)))
{
	p = ggplot(get(valid_polygons[i])) +
	  geom_sf(aes(fill=id), show.legend=F)
	print(p)
}
dev.off()


##### get polygon centroids
valid_centroids = as.character()
cent_coords = as.character()

for(i in seq(length(valid_polygons)))
{
	assign(paste0(dataset_name[i], "_centroids"),st_centroid(get(valid_polygons[i])))
	valid_centroids = c(valid_centroids, paste0(dataset_name[i], "_centroids"))
	assign(paste0(dataset_name[i], "_cent_coords"), st_coordinates(get(paste0(dataset_name[i], "_centroids"))) %>% data.frame() %>% setNames(c("x", "y")))
	cent_coords = c(cent_coords, paste0(dataset_name[i], "_cent_coords"))
	assign(cent_coords[i], cbind(get(cent_coords[i]), cell=get(valid_centroids[i])$id))
}


##### get polygon coordinates
dataset_polygons = as.character()
polygon_coords = as.character()

for(i in seq(length(valid_polygons)))
{
	cell_id = get(valid_polygons[i])$id
	cell_id = as.data.frame(cbind(cell_num=1:length(cell_id), cell=cell_id))
	
	assign(paste0(dataset_name[i], "_polygons"),st_coordinates(get(valid_polygons[i]))[,c("X", "Y", "L2")] %>% data.frame() %>% setNames(c("x", "y", "L2")))
	dataset_polygons = c(dataset_polygons, paste0(dataset_name[i], "_polygons"))
	
	polygon_id = get(dataset_polygons[i])$L2 %>% data.frame()
	polygon_id$cell <- cell_id$cell[match(polygon_id[,1], cell_id$cell_num)]
	assign(dataset_polygons[i], cbind(get(dataset_polygons[i]), cell=polygon_id$cell))
	assign(dataset_polygons[i], get(dataset_polygons[i])[,c("x", "y", "cell")])
}

### load spatial molecule data
molecule_files = c(
"/ceph/ethylene/merscope/baysor/hls1_section1/hls1_section1_region3/detected_transcripts_hls1_r3_240912.csv",
"/ceph/ethylene/merscope/baysor/hls1_section1/hls1_section1_region2/detected_transcripts_hls1_region2_240912.csv", 
"/ceph/ethylene/merscope/baysor/hls1_section1/hls1_section1_region1/detected_transcripts_hls1_region1_240912.csv", 
"/ceph/ethylene/merscope/baysor/hls1_air_section2_region3/detected_transcripts_hls1_air_section2_region3.csv", 
"/ceph/ethylene/merscope/baysor/hls1_air_section2_region2/detected_transcripts_hls1_air_section2_region2.csv", 
"/ceph/ethylene/merscope/baysor/hls1_air_section2_region1/detected_transcripts_hls1_air_section2_region1.csv", 
"/ceph/ethylene/merscope/baysor/Col_ET_section2_region3/detected_transcripts_Col_ET_section2_region3.csv", 
"/ceph/ethylene/merscope/baysor/Col_ET_section2_region2/detected_transcripts_Col_ET_section2_region2.csv", 
"/ceph/ethylene/merscope/baysor/Col_ET_section2_region1/detected_transcripts_Col_ET_section2_region1.csv", 
"/ceph/ethylene/merscope/baysor/Col_ET_section1_region0/detected_transcripts_Col_ET_section1_region0.csv", 
"/ceph/ethylene/merscope/baysor/exp1/section1_region0/detected_transcripts_section1_region0.csv", 
"/ceph/ethylene/merscope/baysor/exp1/detected_transcripts_section1.csv", 
"/ceph/ethylene/merscope/baysor/exp2/Col_air_section2_region0_detected_transcripts.csv"
)


molecules_name = c(
"hls1_air_sect1_r3_molecules", 
"hls1_air_sect1_r2_molecules", 
"hls1_air_sect1_r1_molecules", 
"hls1_air_sect2_r3_molecules", 
"hls1_air_sect2_r2_molecules", 
"hls1_air_sect2_r1_molecules", 
"Col_ET_sect2_r3_molecules", 
"Col_ET_sect2_r2_molecules", 
"Col_ET_sect2_r1_molecules", 
"Col_ET_sect1_r0_molecules", 
"Col_air_sect1_r0_molecules", 
"Col_air_sect1_r1_molecules", 
"Col_air_sect2_r0_molecules"
)

#### reassign variable names

for(i in seq(length(molecule_files)))
{
	assign(molecules_name[i], read.delim(molecule_files[i], header=TRUE, sep = ","))
	assign(molecules_name[i], get(molecules_name[i])[,c("global_x", "global_y", "gene")] %>% `colnames<-`(c("x", "y", "gene")) %>% filter(!str_detect(gene,"Blank*")))
}

for(i in seq(length(molecules_name)))
{
	temp = get(molecules_name[i])
	temp$gene = gene_list$geneId[match(temp$gene, gene_list$Vizgen.Gene)]
	assign(molecules_name[i], temp)
}

gcm_files = c(
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section1_region3_segmentation_noZ_250205/baysor_hls1_air_section1_region3_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section1_region2_segmentation_noZ_250205/baysor_hls1_air_section1_region2_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section1_region1_segmentation_noZ_250205/baysor_hls1_air_section1_region1_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section2_region3_segmentation_noZ_250205/baysor_hls1_air_section2_region3_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section2_region2_segmentation_noZ_250205/baysor_hls1_air_section2_region2_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/hls1_air_section2_region1_segmentation_noZ_250205/baysor_hls1_air_section2_region1_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section2_region3_segmentation_noZ_250205/baysor_Col_ET_section2_region3_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section2_region2_segmentation_noZ_250205/baysor_Col_ET_section2_region2_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section2_region1_segmentation_noZ_250205/baysor_Col_ET_section2_region1_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_ET_section1_region0_segmentation_noZ_250205/baysor_Col_ET_section1_region0_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_air_section1_region0_segmentation_noZ_250205/baysor_Col_air_section1_region0_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_air_section1_region1_segmentation_noZ_250205/baysor_Col_air_section1_region1_segmentation_noZ_250205_counts.tsv", 
"/ceph/ethylene/merscope/Baysor_segmentation_noZ_250205/Col_air_section2_region0_segmentation_noZ_250205/baysor_Col_air_section2_region0_segmentation_noZ_250205_counts.tsv"
)

gcm_name = c(
"hls1_air_sect1_r3_gcm", 
"hls1_air_sect1_r2_gcm", 
"hls1_air_sect1_r1_gcm", 
"hls1_air_sect2_r3_gcm", 
"hls1_air_sect2_r2_gcm", 
"hls1_air_sect2_r1_gcm", 
"Col_ET_sect2_r3_gcm", 
"Col_ET_sect2_r2_gcm", 
"Col_ET_sect2_r1_gcm", 
"Col_ET_sect1_r0_gcm", 
"Col_air_sect1_r0_gcm", 
"Col_air_sect1_r1_gcm", 
"Col_air_sect2_r0_gcm"
)


for(i in seq(length(gcm_files)))
{
	assign(gcm_name[i], read.delim(gcm_files[i], header=TRUE, row.names=1))
}


for(i in seq(length(gcm_name)))
{
	hold = get(gcm_name[i])
	hold = hold[!grepl("Blank*", row.names(hold)),]
	row.names(hold) = gene_list$geneId[match(row.names(hold), gene_list$Vizgen.Gene)]
	colnames(hold) = gsub("\\.", "-", colnames(hold))
	
	cell_id = get(valid_polygons[i])$id
	cell_id = as.data.frame(cbind(cell_num=1:length(cell_id), cell=cell_id))
	
	colnames(hold) %in% cell_id$cell
	hold = hold[colnames(hold) %in% cell_id$cell]
	assign(gcm_name[i], hold)
}


spatial_dir = c(
"/ceph/ethylene/merscope/baysor/hls1_section1/hls1_section1_region3/",
"/ceph/ethylene/merscope/baysor/hls1_section1/hls1_section1_region2/", 
"/ceph/ethylene/merscope/baysor/hls1_section1/hls1_section1_region1/", 
"/ceph/ethylene/merscope/baysor/hls1_air_section2_region3/", 
"/ceph/ethylene/merscope/baysor/hls1_air_section2_region2/", 
"/ceph/ethylene/merscope/baysor/hls1_air_section2_region1/", 
"/ceph/ethylene/merscope/baysor/Col_ET_section2_region3/", 
"/ceph/ethylene/merscope/baysor/Col_ET_section2_region2/", 
"/ceph/ethylene/merscope/baysor/Col_ET_section2_region1/", 
"/ceph/ethylene/merscope/baysor/Col_ET_section1_region0/", 
"/ceph/ethylene/merscope/baysor/exp1/section1_region0/", 
"/ceph/ethylene/merscope/baysor/exp1/", 
"/ceph/ethylene/merscope/baysor/exp2/"
)


for(i in seq(length(spatial_dir)))
{
	assign(paste0(dataset_name[i], "_list"), list(ReadVizgen(spatial_dir[1]), transcripts = get(gcm_name[i]), centroids=get(cent_coords[i]), segmentations=get(dataset_polygons[i]), molecules=get(molecules_name[i])))
}



### sanity check of centroid locations
pdf("centroids_spatial_datasets.pdf")
for(i in seq(length(spatial_dir)))
{
	print(ggplot(get(paste0(dataset_name[i], "_list"))$centroids, aes(x= x, y = y))+
			geom_point(size = 0.1, color = "grey") +
			theme_classic())
}
dev.off()

### sanity check of spatial molecule locations
pdf("molecules_spatial_datasets.pdf")
for(i in seq(length(spatial_dir)))
{
	print(ggplot(get(paste0(dataset_name[i], "_list"))$molecules, aes(x= x, y = y))+
			geom_point(size = 0.1, color = "grey") +
			theme_classic())
}
dev.off()


for(i in seq(length(dataset_name)))
{
transcripts = get(paste0(dataset_name[i], "_list"))$transcripts
assign(paste0(dataset_name[i], "_seurat"), CreateSeuratObject(counts = transcripts, assay = "Vizgen"))
}

for(i in seq(length(dataset_name)))
{
assign(paste0(dataset_name[i], "_seurat_centroids"), CreateCentroids(get(cent_coords[i])))
assign(paste0(dataset_name[i], "_seurat_segmentations"), CreateSegmentation(get(dataset_polygons[i])))
}

for(i in seq(length(dataset_name)))
{
assign(paste0(dataset_name[i], "_segmentation_data"),list(
    "centroids" = get(paste0(dataset_name[i], "_seurat_centroids")),
    "segmentation" = get(paste0(dataset_name[i], "_seurat_segmentations"))))
}


for(i in seq(length(dataset_name)))
{
assign(paste0(dataset_name[i], "_coordinates"), CreateFOV(
												coords = get(paste0(dataset_name[i], "_segmentation_data")),
												type = c("segmentation", "centroids"),
												molecules = get(molecules_name[i]),
												mol.type = "microns",
												assay = "Vizgen"
  ) %>% 
    subset(x = .,
           cells = intersect(x = Cells(x = .[["segmentation"]]),
                               y = Cells(x = paste0(dataset_name[i], "_seurat")))))
}


seurat_datasets = NULL
seurat_datastes = character()

for(i in seq(length(dataset_name)))
{
	seurat_obj = get(paste0(dataset_name[i], "_seurat"))
	seurat_obj[[dataset_name[i]]] = get(paste0(dataset_name[i], "_coordinates"))
	seurat_obj = subset(seurat_obj, nCount_Vizgen > 10)
	assign(paste0(dataset_name[i], "_seurat"), seurat_obj)
	seurat_datasets = c(seurat_datasets, paste0(dataset_name[i], "_seurat"))
}

#### sanity check of spatial loadings into Seurat
pdf("seurat_spatial_datasets.pdf")
for(i in seq(length(dataset_name)))
{
print(ImageDimPlot(get(paste0(dataset_name[i], "_seurat")), fov = dataset_name[i], cols = "polychrome", axes = TRUE) + ggtitle(dataset_name[i]))
}
dev.off()



for(dataset in seurat_datasets)
{
	assign(dataset, SCTransform(get(dataset), assay = "Vizgen", clip.range = c(-10, 10), variable.features.n = nrow(get(dataset))) %>%
			RunPCA(npcs = 30) %>%
			RunUMAP(dims = 1:30) %>%
			FindNeighbors(reduction = "pca", dims = 1:30) %>%
			FindClusters(resolution = 0.3))
}



pdf("segmented_clusters.pdf")
for(i in seq(length(dataset_name)))
{
	print(DimPlot(get(seurat_datasets[i]), reduction = "umap", label=T))
	print(ImageDimPlot(get(seurat_datasets[i]), fov = dataset_name[i], cols = "polychrome", axes = TRUE, border.color = "NA", alpha = 0.8) + ggtitle(dataset_name[i]))
}
dev.off()


hls1_air_sect1_r3_seurat$assay = "merfish"
hls1_air_sect1_r2_seurat$assay = "merfish"
hls1_air_sect1_r1_seurat$assay = "merfish"
hls1_air_sect2_r3_seurat$assay = "merfish"
hls1_air_sect2_r2_seurat$assay = "merfish"
hls1_air_sect2_r1_seurat$assay = "merfish"
Col_ET_sect2_r3_seurat$assay = "merfish"
Col_ET_sect2_r2_seurat$assay = "merfish"
Col_ET_sect2_r1_seurat$assay = "merfish"
\Col_ET_sect1_r0_seurat$assay = "merfish"
Col_air_sect1_r0_seurat$assay = "merfish"
Col_air_sect1_r1_seurat$assay = "merfish"
Col_air_sect2_r0_seurat$assay = "merfish"

hls1_air_sect1_r3_seurat$orig.ident = "merfish"
hls1_air_sect1_r2_seurat$orig.ident = "merfish"
hls1_air_sect1_r1_seurat$orig.ident = "merfish"
hls1_air_sect2_r3_seurat$orig.ident = "merfish"
hls1_air_sect2_r2_seurat$orig.ident = "merfish"
hls1_air_sect2_r1_seurat$orig.ident = "merfish"
\Col_ET_sect2_r3_seurat$orig.ident = "merfish"
Col_ET_sect2_r2_seurat$orig.ident = "merfish"
Col_ET_sect2_r1_seurat$orig.ident = "merfish"
\Col_ET_sect1_r0_seurat$orig.ident = "merfish"
Col_air_sect1_r0_seurat$orig.ident = "merfish"
Col_air_sect1_r1_seurat$orig.ident = "merfish"
Col_air_sect2_r0_seurat$orig.ident = "merfish"

hls1_air_sect1_r3_seurat$geno = "hls1"
hls1_air_sect1_r2_seurat$geno = "hls1"
hls1_air_sect1_r1_seurat$geno = "hls1"
hls1_air_sect2_r3_seurat$geno = "hls1"
hls1_air_sect2_r2_seurat$geno = "hls1"
hls1_air_sect2_r1_seurat$geno = "hls1"
Col_ET_sect2_r3_seurat$geno = "Col"
Col_ET_sect2_r2_seurat$geno = "Col"
Col_ET_sect2_r1_seurat$geno = "Col"
Col_ET_sect1_r0_seurat$geno = "Col"
Col_air_sect1_r0_seurat$geno = "Col"
Col_air_sect1_r1_seurat$geno = "Col"
Col_air_sect2_r0_seurat$geno = "Col"

hls1_air_sect1_r3_seurat$stim = "CTRL"
hls1_air_sect1_r2_seurat$stim = "CTRL"
hls1_air_sect1_r1_seurat$stim = "CTRL"
hls1_air_sect2_r3_seurat$stim = "CTRL"
hls1_air_sect2_r2_seurat$stim = "CTRL"
hls1_air_sect2_r1_seurat$stim = "CTRL"
Col_ET_sect2_r3_seurat$stim = "ET"
Col_ET_sect2_r2_seurat$stim = "ET"
Col_ET_sect2_r1_seurat$stim = "ET"
Col_ET_sect1_r0_seurat$stim = "ET"
Col_air_sect1_r0_seurat$stim = "CTRL"
Col_air_sect1_r1_seurat$stim = "CTRL"
Col_air_sect2_r0_seurat$stim = "CTRL"


### save RDS file of Seurat spatial datasets
for(i in seq(length(dataset_name)))
{
saveRDS(get(seurat_datasets[i]), file=paste0("spatial_dataset_", dataset_name[i], "_250205.rds"))
}


####### integrate spatial datasets

merged = merge(hls1_air_sect1_r3_seurat, c(
hls1_air_sect1_r2_seurat,
hls1_air_sect1_r1_seurat,
hls1_air_sect2_r3_seurat,
hls1_air_sect2_r2_seurat,
hls1_air_sect2_r1_seurat,
Col_ET_sect2_r3_seurat,
Col_ET_sect2_r2_seurat,
Col_ET_sect2_r1_seurat,
Col_ET_sect1_r0_seurat,
Col_air_sect1_r0_seurat,
Col_air_sect1_r1_seurat,
Col_air_sect2_r0_seurat
))


merged[["RNA"]] = merged[["Vizgen"]]
DefaultAssay(merged) = "RNA"

#### spatial only clustering
merged <- SCTransform(merged, assay = "RNA", clip.range = c(-10, 10))
merged <- RunPCA(merged, npcs = 30, features = rownames(merged))
merged <- RunUMAP(merged, dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)


datasets = Images(merged)

pdf("spatial_merged_clusters.pdf")
print(DimPlot(merged, reduction = "umap", label=T))
for(i in seq(length(datasets)))
{
	print(ImageDimPlot(merged, fov = datasets[i], cols = "polychrome", axes = TRUE, border.color = "NA", alpha = 0.8) + ggtitle(dataset_name[i]))
}
dev.off()



pdf("fig2A_spatial_integration.pdf")
DimPlot(merged, cols="polychrome")
dev.off()

saveRDS(merged, "spatial_only_integration_250205.rds")
