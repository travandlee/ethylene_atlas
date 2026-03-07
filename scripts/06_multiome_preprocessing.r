#### Multiome dataset preprocessing

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


system("mkdir /ceph/ethylene/multiome/integration_241121")

setwd("/ceph/ethylene/multiome/integration_241121")

###  Tn5 blacklist regions
gr_obj = import("/gale/netapp/home/trlee/general/Tn5_rRNA_blacklist.bed")
gr_list = split(gr_obj, gr_obj$name)

### load annotations file
annotations <- rtracklayer::import("/gale/netapp/home/trlee/reference/multiome/Araport11_GTF_genes_transposons.current.gtf")
genome(annotations) <- "TAIR10"
mcols(annotations)[,2] = as.character(mcols(annotations)[,2])
names(mcols(annotations))[5] = "tx_id"
annotations$gene_name = annotations$gene_id
annotations$gene_biotype = annotations$type
annotations = annotations[annotations@seqnames != "ChrC" & annotations@seqnames != "ChrM"]

#######
## Seedling 72h sample 1
#######



inputdata.10x <- Read10X_h5("/ceph/ethylene/multiome/processed/seedling_72h_seedling_1_241007_reseq/outs/filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/ceph/ethylene/multiome/processed/seedling_72h_seedling_1_241007_reseq/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )

multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(multiome_dataset, file="/ceph/ethylene/multiome/processed/seedling_72h_1_filtered_241121.rds")

#######
## Seedling 72h sample 2
#######

inputdata.10x <- Read10X_h5("/ceph/ethylene/multiome/processed/seedling_72h_seedling_2_241007_reseq/outs/raw_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/ceph/ethylene/multiome/processed/seedling_72h_seedling_2_241007_reseq/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )

multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(multiome_dataset, file="/ceph/ethylene/multiome/processed/seedling_72h_2_filtered_241121.rds")

#######
## Seedling 72h + ET sample 2
#######

setwd("/ceph/ethylene/multiome/processed")

inputdata.10x <- Read10X_h5("/ceph/ethylene/multiome/seedling_ET_seedling_250805/outs/raw_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/ceph/ethylene/multiome/seedling_ET_seedling_250805/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )


multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(multiome_dataset, file="/ceph/ethylene/multiome/processed/seedling_ET_2_filtered_241121.rds")


######
inputdata.10x <- Read10X_h5("/gale/ddn/ddn_neomorph/trlee/et/hook/seedling_2do_seedling_240111/outs/filtered_feature_bc_matrix.h5")


rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/gale/ddn/ddn_neomorph/trlee/et/hook/seedling_2do_seedling_240111/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )


multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(multiome_dataset, file="seedling_48h_raw_filtered_240708.rds")




#####

inputdata.10x <- Read10X_h5("/gale/ddn/ddn_neomorph/trlee/et/hook/seedling_2do_seedling_240111/outs/filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/gale/ddn/ddn_neomorph/trlee/et/hook/seedling_2do_seedling_240111/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )


multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(pbmc, file="seedling_48h_raw_filtered_240708.rds")


######
inputdata.10x <- Read10X_h5("/ceph/ethylene/multiome/seedling_28h_seedling_240206/outs/filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/ceph/ethylene/multiome/seedling_28h_seedling_240206/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )

multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(pbmc, file="seedling_28h_240603.rds")




#####


inputdata.10x <- Read10X_h5("/gale/ddn/ddn_neomorph/trlee/et/multiome/seedling_4ET_230404/outs/filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

multiome_dataset <- CreateSeuratObject(counts = rna_counts, project = "seedling_72h")
multiome_dataset[["percent.mt"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATMG")
multiome_dataset[["percent.cp"]] <- PercentageFeatureSet(multiome_dataset, pattern = "ATCG")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)[1:5]
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "/gale/ddn/ddn_neomorph/trlee/et/multiome/seedling_4ET_230404/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = "TAIR10",
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations,
   validate.fragments = FALSE
 )

multiome_dataset[["ATAC"]] <- chrom_assay
DefaultAssay(multiome_dataset) = "ATAC"
multiome_dataset = TSSEnrichment(multiome_dataset, tss.positions = StringToGRanges(rownames(multiome_dataset), sep = c("-", "-")))
multiome_dataset <- NucleosomeSignal(multiome_dataset)

multiome_dataset$blacklist_fraction <- FractionCountsInRegion(multiome_dataset, regions = gr_list)

multiome_dataset <- subset(x = multiome_dataset, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nCount_RNA < 25000 & nCount_RNA > 500)

DefaultAssay(multiome_dataset) <- "RNA"
multiome_dataset <- SCTransform(multiome_dataset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(pbmc, file="ET_4_230411.rds")


