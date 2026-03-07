#### multiome integration
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
library(harmony)


system("mkdir /ceph/ethylene/multiome/peaks")
system("mkdir /ceph/ethylene/multiome/integration_250805")
setwd("/ceph/ethylene/multiome/integration_250805")

# load datasets
air28 = readRDS(file="/ceph/ethylene/multiome/seedling_28h_raw_filtered_240705.rds")
air48 = readRDS(file="/ceph/ethylene/multiome/seedling_48h_raw_filtered_240708.rds")
air0 = readRDS(file="/ceph/ethylene/multiome/Air0_raw_filtered_240708.rds")
ET4 = readRDS(file="/ceph/ethylene/multiome/ET4_raw_filtered_240708.rds")
ET4_2 = readRDS(file="/ceph/ethylene/multiome/processed/seedling_ET_2_filtered_241121.rds")
air72_1 = readRDS(file="/ceph/ethylene/multiome/processed/seedling_72h_1_filtered_241121.rds")
air72_2 = readRDS(file="/ceph/ethylene/multiome/processed/seedling_72h_2_filtered_241121.rds")

air0 = RenameCells(air0, add.cell.id="air0")
air28 = RenameCells(air28, add.cell.id="air28")
air48 = RenameCells(air48, add.cell.id="air48")
ET4 = RenameCells(ET4, add.cell.id="ET4")
ET4_2 = RenameCells(ET4, add.cell.id="ET4")
air72_1 = RenameCells(air72_1, add.cell.id="air72")
air72_2 = RenameCells(air72_2, add.cell.id="air72")

DefaultAssay(air0) = "ATAC"
DefaultAssay(air28) = "ATAC"
DefaultAssay(air48) = "ATAC"
DefaultAssay(ET4) = "ATAC"
DefaultAssay(ET4_2) = "ATAC"
DefaultAssay(air72_1) = "ATAC"
DefaultAssay(air72_2) = "ATAC"

###call peaks and export to .bed file for all datasets independently
blacklist = import("/gale/netapp/home/trlee/general/Tn5_rRNA_blacklist.bed")

peaks_cluster <- CallPeaks(air28, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster <- dropSeqlevels(peaks_cluster, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster <- subsetByOverlaps(x = peaks_cluster, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster, "/ceph/ethylene/multiome/peaks/seedling_28h_peaks.bed")

peaks_cluster_air48 <- CallPeaks(air48, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster_air48 <- dropSeqlevels(peaks_cluster_air48, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster_air48 <- subsetByOverlaps(x = peaks_cluster_air48, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster_air48, "/ceph/ethylene/multiome/peaks/air48_peaks.bed")

peaks_cluster_ET <- CallPeaks(ET4, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster_ET <- dropSeqlevels(peaks_cluster_ET, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster_ET <- subsetByOverlaps(x = peaks_cluster_ET, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster_ET, "/ceph/ethylene/multiome/peaks/ET4_peaks.bed")

peaks_cluster_ET_2 <- CallPeaks(ET4_2, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster_ET_2 <- dropSeqlevels(peaks_cluster_ET_2, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster_ET_2 <- subsetByOverlaps(x = peaks_cluster_ET_2, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster_ET_2, "/ceph/ethylene/multiome/peaks/ET4_2_peaks.bed")

peaks_cluster_air0 <- CallPeaks(air0, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster_air0 <- dropSeqlevels(peaks_cluster_air0, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster_air0 <- subsetByOverlaps(x = peaks_cluster_air0, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster_air0, "/ceph/ethylene/multiome/peaks/Air0_peaks.bed")

peaks_cluster_air72_1 <- CallPeaks(air72_1, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster_air72_1 <- dropSeqlevels(peaks_cluster_air72_1, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster_air72_1 <- subsetByOverlaps(x = peaks_cluster_air72_1, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster_air72_1, "/ceph/ethylene/multiome/peaks/Air72_1_peaks_241121.bed")

peaks_cluster_air72_2 <- CallPeaks(air72_2, group.by = "seurat_clusters", effective.genome.size = 119481543)
peaks_cluster_air72_2 <- dropSeqlevels(peaks_cluster_air72_2, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks_cluster_air72_2 <- subsetByOverlaps(x = peaks_cluster_air72_2, ranges = blacklist, invert = TRUE)
rtracklayer::export.bed(peaks_cluster_air72_2, "/ceph/ethylene/multiome/peaks/Air72_2_peaks_241121.bed")

peaks.28h <- read.table(
  file = "/ceph/ethylene/multiome/peaks/seedling_28h_peaks.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)

peaks.air48 <- read.table(
  file = "/ceph/ethylene/multiome/peaks/air48_peaks.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)

peaks.ET4 <- read.table(
  file = "/ceph/ethylene/multiome/peaks/ET4_peaks.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)

peaks.ET4_2 <- read.table(
  file = "/ceph/ethylene/multiome/peaks/ET4_2_peaks.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)
peaks.air0 <- read.table(
  file = "/ceph/ethylene/multiome/peaks/Air0_peaks.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)

peaks.air72_1 <- read.table(
  file = "/ceph/ethylene/multiome/peaks/Air72_1_peaks_241121.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)

peaks.air72_2 <- read.table(
  file = "/ceph/ethylene/multiome/peaks/Air72_2_peaks_241121.bed",
  col.names = c("Chr", "start", "end", "name", "score", "strand")
)

gr.28h <- makeGRangesFromDataFrame(peaks.28h)
gr.48h <- makeGRangesFromDataFrame(peaks.air48)

gr.ET4 <- makeGRangesFromDataFrame(peaks.ET4)
gr.ET4_2 <- makeGRangesFromDataFrame(peaks.ET4_2)

gr.air0 <- makeGRangesFromDataFrame(peaks.air0)

gr.72h_1 <- makeGRangesFromDataFrame(peaks.air72_1)
gr.72h_2 <- makeGRangesFromDataFrame(peaks.air72_2)


combined.peaks <- reduce(x = c(gr.28h, gr.48h, gr.ET4, gr.ET4_2, gr.air0, gr.72h_1, gr.72h_2))

### collapse overlapping peaks
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


frags.seedling_28h = Fragments(air28)

seedling_28h.counts <- FeatureMatrix(
  fragments = frags.seedling_28h,
  features = combined.peaks
  )


frags.air48 = Fragments(air48)

air48.counts <- FeatureMatrix(
  fragments = frags.air48,
  features = combined.peaks
  )


frags.ET4 = Fragments(ET4)

ET4.counts <- FeatureMatrix(
  fragments = frags.ET4,
  features = combined.peaks
  )

frags.ET4_2 = Fragments(ET4_2)

ET4_2.counts <- FeatureMatrix(
  fragments = frags.ET4_2,
  features = combined.peaks
  )


frags.air0 = Fragments(air0)

air0.counts <- FeatureMatrix(
  fragments = frags.air0,
  features = combined.peaks
  )


frags.air72_1 = Fragments(air72_1)

air72_1.counts <- FeatureMatrix(
  fragments = frags.air72_1,
  features = combined.peaks
  )

frags.air72_2 = Fragments(air72_2)

air72_2.counts <- FeatureMatrix(
  fragments = frags.air72_2,
  features = combined.peaks
  )


seedling_28h_assay <- CreateChromatinAssay(seedling_28h.counts, fragments = frags.seedling_28h)
air28[["ATAC"]] <- seedling_28h_assay

seedling_48h_assay <- CreateChromatinAssay(air48.counts, fragments = frags.air48)
air48[["ATAC"]] <- seedling_48h_assay

ET4_assay <- CreateChromatinAssay(ET4.counts, fragments = frags.ET4)
ET4[["ATAC"]] <- ET4_assay

ET4_2_assay <- CreateChromatinAssay(ET4.counts, fragments = frags.ET4_2)
ET4_2[["ATAC"]] <- ET4_2_assay

air0_assay <- CreateChromatinAssay(air0.counts, fragments = frags.air0)
air0[["ATAC"]] <- air0_assay

air72_1_assay <- CreateChromatinAssay(air72_1.counts, fragments = frags.air72_1)
air72_1[["ATAC"]] <- air72_1_assay

air72_2_assay <- CreateChromatinAssay(air72_2.counts, fragments = frags.air72_2)
air72_2[["ATAC"]] <- air72_2_assay


combined <- merge(
  x = air28,
  y = list(air48, ET4, ET4_2, air0, air72_1, air72_2)
)


combined[["ATAC"]]
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q5')
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


#### Add annotation information to ATAC assay
annotations <- rtracklayer::import("/gale/netapp/home/trlee/reference/multiome/Araport11_GTF_genes_transposons.current.gtf")

genome(annotations) <- "TAIR10"
mcols(annotations)[,2] = as.character(mcols(annotations)[,2])
names(mcols(annotations))[5] = "tx_id"
annotations$gene_name = annotations$gene_id
annotations$gene_biotype = annotations$type
annotations = annotations[annotations@seqnames != "ChrC" & annotations@seqnames != "ChrM"]
Annotation(combined) = annotations



DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined, verbose = FALSE) %>% RunPCA()
combined = RunHarmony(combined, "orig.ident")
combined = RunUMAP(combined, dims = 1:50, reduction = "harmony", reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


DimPlot(combined, reduction = "umap.rna", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")



combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

combined$multi_clusters = combined$seurat_clusters


combined = FindNeighbors(combined, dims = 1:30, reduction = "harmony", assay = "SCT", verbose = FALSE) %>% FindClusters(verbose = FALSE, resolution = .4)


#############
#### call peaks of merged dataset
#############

blacklist_gr = import("/gale/netapp/home/trlee/general/Tn5_rRNA_blacklist.bed")

blacklist = split(blacklist_gr, blacklist_gr$name)


# call peaks using MACS2
DefaultAssay(combined) = "ATAC"

peaks <- CallPeaks(combined, assay = "ATAC")
peaks <- dropSeqlevels(peaks, value = c("ChrM", "ChrC"), pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)

# quantify counts in each peak
DefaultAssay(combined) = "ATAC"

macs2_counts <- FeatureMatrix(
  fragments = Fragments(combined),
  features = peaks,
  cells = colnames(combined)
)


fragments_ET4 <- CreateFragmentObject(
   path = "/gale/ddn/ddn_neomorph/trlee/et/multiome/seedling_4ET_230404/outs/atac_fragments.tsv.gz",
   cells = colnames(ET4),
   validate.fragments = FALSE
 )


fragments_28h <- CreateFragmentObject(
   path = "/ceph/ethylene/multiome/seedling_28h_seedling_240206/outs/atac_fragments.tsv.gz",
   cells = colnames(air28),
   validate.fragments = FALSE
 )


fragments_48h <- CreateFragmentObject(
   path = "/gale/ddn/ddn_neomorph/trlee/et/hook/seedling_2do_seedling_240111/outs/atac_fragments.tsv.gz",
   cells = colnames(air48),
   validate.fragments = FALSE
 )


fragments_0h <- CreateFragmentObject(
   path = "/gale/ddn/ddn_neomorph/trlee/et/multiome/seedling_0ET_230404/outs/atac_fragments.tsv.gz",
   cells = colnames(air0),
   validate.fragments = FALSE
 )


fragments_72h_1 <- CreateFragmentObject(
   path = "/ceph/ethylene/multiome/processed/seedling_72h_seedling_1_241007_reseq/outs/atac_fragments.tsv.gz",
   cells = colnames(air72_1),
   validate.fragments = FALSE
 )


fragments_72h_2 <- CreateFragmentObject(
   path = "/ceph/ethylene/multiome/processed/seedling_72h_seedling_2_241007_reseq/outs/atac_fragments.tsv.gz",
   cells = colnames(air72_2),
   validate.fragments = FALSE
 )


# create a new assay using the MACS2 peak set and add it to the Seurat object

combined[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = c(fragments_ET4, fragments_28h, fragments_48h, fragments_0h, fragments_72h_1, fragments_72h_2),
  annotation = annotations
)

DefaultAssay(combined) <- "peaks"

combined <- RegionStats(combined, genome = BSgenome.Athaliana.TAIR.TAIR9)


saveRDS(combined, file="multiome_integrated_timecourse_250825.rds")
