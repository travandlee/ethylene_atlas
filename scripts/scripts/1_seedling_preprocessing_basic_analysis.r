### Analysis of all single-nucleus RNA-seq datasets

library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(reshape2)
library(dplyr)

#### set working directory
setwd("")

system("mkdir results")
system("mkdir results/harmony")


ctrl.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/10x_rna_combined/control_premrna_transcriptome_reseq/outs/filtered_feature_bc_matrix")
ctrl_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/042919_10x_rna_ET/data/control_rna_042919_hiseq/outs/filtered_feature_bc_matrix")
ctrl_3.data <- Read10X(data.dir = '/gale/raidix/rdx-7/trlee/080719_10x_rna/control_rna_080719_hiseq/outs/filtered_feature_bc_matrix')
ctrl_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200718_col_hls_rna/Col_24air_rna_200709_novaseq/outs/filtered_feature_bc_matrix")
ctrl_5.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200724_col_hls_rna/Col_24air_rna_200716_novaseq/outs/filtered_feature_bc_matrix")
ctrl_6.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210426_col_air_rna/Col_24air_1_rna_210427_novaseq/outs/filtered_feature_bc_matrix")
ctrl_7.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210426_col_air_rna/Col_24air_2_rna_210427_novaseq/outs/filtered_feature_bc_matrix")
ctrl_8.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210426_col_air_rna/Col_24air_3_rna_210427_novaseq/outs/filtered_feature_bc_matrix")
ctrl_9.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210426_col_air_rna/Col_24air_4_rna_210427_novaseq/outs/filtered_feature_bc_matrix")

### 24h ET treatment
et_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/10x_rna_combined/ethylene_premrna_transcriptome_reseq/outs/filtered_feature_bc_matrix")
et_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/042919_10x_rna_ET/data/ethylene_rna_042919_hiseq/outs/filtered_feature_bc_matrix")
et_3.data <- Read10X(data.dir = '/gale/raidix/rdx-7/trlee/080719_10x_rna/ethylene_rna_080719_hiseq/outs/filtered_feature_bc_matrix')
ET_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200718_col_hls_rna/Col_24ET_rna_200709_novaseq/outs/filtered_feature_bc_matrix")
ET_5.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200724_col_hls_rna/Col_24ET_rna_200716_novaseq/outs/filtered_feature_bc_matrix")

### 72h ET treatment
col_72ET_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210713_col_3dET/Col_72ET_1_rna_210713_novaseq/outs/filtered_feature_bc_matrix")
col_72ET_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210713_col_3dET/Col_72ET_2_rna_210713_novaseq/outs/filtered_feature_bc_matrix")
col_72ET_3.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210713_col_3dET/Col_72ET_3_rna_210713_novaseq/outs/filtered_feature_bc_matrix")

### ET timecourse datasets
ET0.data <- Read10X(data.dir = '/gale/raidix/rdx-7/trlee/080719_10x_rna/0hr_rna_080719_hiseq/outs/filtered_feature_bc_matrix')
ET30.data <- Read10X(data.dir = '/gale/raidix/rdx-7/trlee/080719_10x_rna/30min_rna_080719_hiseq/outs/filtered_feature_bc_matrix')
ET4.data <- Read10X(data.dir = '/gale/raidix/rdx-7/trlee/080719_10x_rna/4hr_rna_080719_hiseq/outs/filtered_feature_bc_matrix')
ET24.data <- Read10X(data.dir = '/gale/raidix/rdx-7/trlee/080719_10x_rna/24hr_rna_080719_hiseq/outs/filtered_feature_bc_matrix')

ET30_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200116_10x_ET/200116_ET_rna/30minET_hook_rna_200116_novaseq/outs/filtered_feature_bc_matrix")
ET4_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200116_10x_ET/200116_ET_rna/4ET_hook_rna_200116_novaseq/outs/filtered_feature_bc_matrix")
ET8.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200116_10x_ET/200116_ET_rna/8ET_hook_rna_200116_novaseq/outs/filtered_feature_bc_matrix")

col_2ET_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210713_col_ETcourse/Col_2ET_1_rna_210713_novaseq/outs/filtered_feature_bc_matrix")
col_4ET_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210713_col_ETcourse/Col_4ET_rna_210713_novaseq/outs/filtered_feature_bc_matrix")
col_4ET_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/220106_col_ET/Col_4ET_rna_220106_novaseq/outs/filtered_feature_bc_matrix")
col_8ET_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210713_col_ETcourse/Col_8ET_rna_210713_novaseq/outs/filtered_feature_bc_matrix")
col_72ET_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210827_col_3dET/Col_72ET_4_rna_210827_novaseq/outs/filtered_feature_bc_matrix")
col_72ET_5.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210827_col_3dET/Col_72ET_5_rna_210827_novaseq/outs/filtered_feature_bc_matrix")
col_8ET_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210827_col_3dET/Col_8ET_rna_210827_novaseq/outs/filtered_feature_bc_matrix")

### apical hook dissection datasets
air.data <- Read10X(data.dir =  "/ceph/ethylene/TL_10X_RNA_apical_hook_230918/cellranger3/24air_hook_rna_231220_reseq_novaseq/outs/filtered_feature_bc_matrix")
et.data <- Read10X(data.dir = "/ceph/ethylene/TL_10X_RNA_apical_hook_230918/cellranger3/24et_hook_rna_231220_reseq_novaseq/outs/filtered_feature_bc_matrix")

### hls1 datasets
hls1_air.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200303_hls1_rna/hls1_24air_rna_200303_novaseq/outs/filtered_feature_bc_matrix")
hls1_ET.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200303_hls1_rna/hls1_24ET_rna_200303_novaseq/outs/filtered_feature_bc_matrix")

hls1_air_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200323_hls1_rna/hls1_24air_rna_200323_novaseq/outs/filtered_feature_bc_matrix")
hls1_ET_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200323_hls1_rna/hls1_24ET_rna_200323_novaseq/outs/filtered_feature_bc_matrix")

hls1_air_3.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200718_col_hls_rna/hls1_24air_rna_200709_novaseq/outs/filtered_feature_bc_matrix")
hls1_ET_3.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200718_col_hls_rna/hls1_24ET_rna_200709_novaseq/outs/filtered_feature_bc_matrix")

hls1_air_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200724_col_hls_rna/hls1_24air_rna_200716_novaseq/outs/filtered_feature_bc_matrix")
hls1_ET_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200724_col_hls_rna/hls1_24ET_rna_200716_novaseq/outs/filtered_feature_bc_matrix")

### ein2-5 datasets
ein2_air_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200816_ein2_rna/ein2_24air_rna_200814_novaseq/outs/filtered_feature_bc_matrix")
ein2_ET_1.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/200816_ein2_rna/ein2_24ET_rna_200814_novaseq/outs/filtered_feature_bc_matrix")

ein2_air_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/201122_ein2_rna/ein2_24air_rna_201103_novaseq/outs/filtered_feature_bc_matrix")
ein2_ET_2.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/201122_ein2_rna/ein2_24ET_rna_201103_novaseq/outs/filtered_feature_bc_matrix")

ein2_air_3.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/201122_ein2_rna/ein2_24air_rna_201030_novaseq/outs/filtered_feature_bc_matrix")
ein2_ET_3.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/201122_ein2_rna/ein2_24ET_rna_201030_novaseq/outs/filtered_feature_bc_matrix")

ein2_air_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210309_ein2_rna/ein2_24air_rna_210131_novaseq/outs/filtered_feature_bc_matrix")
ein2_ET_4.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210309_ein2_rna/ein2_24ET_rna_210131_novaseq/outs/filtered_feature_bc_matrix")

ein2_air_5.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210324_ein2_rna/resequence/ein2_24air_1_rna_resequence_210427_novaseq/outs/filtered_feature_bc_matrix")
ein2_ET_5.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210324_ein2_rna/resequence/ein2_24ET_1_rna_resequence_210427_novaseq/outs/filtered_feature_bc_matrix")

ein2_air_6.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210324_ein2_rna/resequence/ein2_24air_2_rna_resequence_210427_novaseq/outs/filtered_feature_bc_matrix")
ein2_ET_6.data <- Read10X(data.dir = "/gale/raidix/rdx-7/trlee/210324_ein2_rna/resequence/ein2_24ET_2_rna_resequence_210427_novaseq/outs/filtered_feature_bc_matrix")



seedling <- CreateSeuratObject(counts = cbind(ctrl.data, et_1.data, ctrl_2.data, et_2.data, ctrl_3.data, et_3.data, 
	ET0.data, ET30.data, ET4.data, ET24.data, 
	air.data, et.data, 
	ET30_2.data, ET4_2.data, ET8.data, 
	hls1_air.data, hls1_ET.data, hls1_air_2.data, hls1_ET_2.data, 
	ctrl_4.data, ET_4.data, 
	hls1_air_3.data, hls1_ET_3.data, 
	ctrl_5.data, ET_5.data, 
	hls1_air_4.data, hls1_ET_4.data, 
	ein2_air_1.data, ein2_ET_1.data, ein2_air_2.data, ein2_ET_2.data, ein2_air_4.data, ein2_ET_4.data, ein2_air_5.data, ein2_ET_5.data, ein2_air_6.data, ein2_ET_6.data, 
	ctrl_6.data, ctrl_7.data, ctrl_8.data, ctrl_9.data,
	col_72ET_1.data, col_72ET_2.data, col_72ET_3.data,
	col_2ET_1.data, col_4ET_1.data, col_8ET_1.data,
	col_72ET_4.data, col_72ET_5.data, col_8ET_2.data,
	col_4ET_2.data
	), project = "ethylene") 



### identify low quality and doublet cell cutoff
pdf("results/nFeature_dropout.pdf")
ggplot(as.data.frame(seedling@meta.data), aes(x=nFeature_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 350)

ggplot(as.data.frame(seedling@meta.data), aes(x=nFeature_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 350) +
	xlim(0, 1000)

ggplot(as.data.frame(seedling@meta.data), aes(x=nFeature_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 275)

ggplot(as.data.frame(seedling@meta.data), aes(x=nFeature_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 275) +
	xlim(0, 1000)

ggplot(as.data.frame(seedling@meta.data), aes(x=nFeature_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 200)

ggplot(as.data.frame(seedling@meta.data), aes(x=nFeature_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 200) +
	xlim(0, 1000)
dev.off()



pdf("results/nCount_dropout.pdf")
ggplot(as.data.frame(seedling@meta.data), aes(x=nCount_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 400)

ggplot(as.data.frame(seedling@meta.data), aes(x=nCount_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 400) +
	xlim(0, 1000)

ggplot(as.data.frame(seedling@meta.data), aes(x=nCount_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 325)

ggplot(as.data.frame(seedling@meta.data), aes(x=nCount_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 325) +
	xlim(0, 1000)

ggplot(as.data.frame(seedling@meta.data), aes(x=nCount_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 200)

ggplot(as.data.frame(seedling@meta.data), aes(x=nCount_RNA)) + 
	geom_histogram(binwidth = 5) + 
	scale_x_continuous(limits = c(0,10000)) +
	geom_vline(color = "red", linetype = "solid", xintercept = 200) +
	xlim(0, 1000)

dev.off()



seedling@meta.data$stim <- c(rep("CTRL", ncol(ctrl.data)), rep("ET", ncol(et_1.data)), rep("CTRL", ncol(ctrl_2.data)), rep("ET", ncol(et_2.data)), rep("CTRL", ncol(ctrl_3.data)), rep("ET", ncol(et_3.data)), rep("ET", ncol(ET0.data)), rep("ET", ncol(ET30.data)), rep("ET", ncol(ET4.data)), rep("ET", ncol(ET24.data)), 
							rep("CTRL", ncol(air.data)), rep("ET", ncol(et.data)), 
							rep("ET", ncol(ET30_2.data)), rep("ET", ncol(ET4_2.data)), rep("ET", ncol(ET8.data)), rep("CTRL", ncol(hls1_air.data)), rep("ET", ncol(hls1_ET.data)), rep("CTRL", ncol(hls1_air_2.data)), rep("ET", ncol(hls1_ET_2.data)), 
							rep("CTRL", ncol(ctrl_4.data)), rep("ET", ncol(ET_4.data)), 
							rep("CTRL", ncol(hls1_air_3.data)), rep("ET", ncol(hls1_ET_3.data)),
							rep("CTRL", ncol(ctrl_5.data)), rep("ET", ncol(ET_5.data)), 
							rep("CTRL", ncol(hls1_air_4.data)), rep("ET", ncol(hls1_ET_4.data)),
							rep("CTRL", ncol(ein2_air_1.data)), rep("ET", ncol(ein2_ET_1.data)),
							rep("CTRL", ncol(ein2_air_2.data)), rep("ET", ncol(ein2_ET_2.data)),
							rep("CTRL", ncol(ein2_air_4.data)), rep("ET", ncol(ein2_ET_4.data)),
							rep("CTRL", ncol(ein2_air_5.data)), rep("ET", ncol(ein2_ET_5.data)),
							rep("CTRL", ncol(ein2_air_6.data)), rep("ET", ncol(ein2_ET_6.data)),
							rep("CTRL", ncol(ctrl_6.data)), rep("CTRL", ncol(ctrl_7.data)),
							rep("CTRL", ncol(ctrl_8.data)), rep("CTRL", ncol(ctrl_9.data)),
							rep("ET", ncol(col_72ET_1.data)), rep("ET", ncol(col_72ET_2.data)),
							rep("ET", ncol(col_72ET_3.data)), rep("ET", ncol(col_2ET_1.data)),
							rep("ET", ncol(col_4ET_1.data)), rep("ET", ncol(col_8ET_1.data)),
							rep("ET", ncol(col_72ET_4.data)),  rep("ET", ncol(col_72ET_5.data)),
							rep("ET", ncol(col_8ET_2.data)),
							rep("ET", ncol(col_4ET_2.data))
							)

seedling@meta.data$tissue <- c(rep("seedling", ncol(ctrl.data)), rep("seedling", ncol(et_1.data)), rep("seedling", ncol(ctrl_2.data)), rep("seedling", ncol(et_2.data)), rep("seedling", ncol(ctrl_3.data)), rep("seedling", ncol(et_3.data)), rep("seedling", ncol(ET0.data)), rep("seedling", ncol(ET30.data)), rep("seedling", ncol(ET4.data)), rep("seedling", ncol(ET24.data)), 
							rep("hook", ncol(air.data)), rep("hook", ncol(et.data)), 
							rep("seedling", ncol(ET30_2.data)), rep("seedling", ncol(ET4_2.data)), rep("seedling", ncol(ET8.data)), rep("seedling", ncol(hls1_air.data)), rep("seedling", ncol(hls1_ET.data)), rep("seedling", ncol(hls1_air_2.data)), rep("seedling", ncol(hls1_ET_2.data)),
							rep("seedling", ncol(ctrl_4.data)), rep("seedling", ncol(ET_4.data)), 
							rep("seedling", ncol(hls1_air_3.data)), rep("seedling", ncol(hls1_ET_3.data)),
							rep("seedling", ncol(ctrl_5.data)), rep("seedling", ncol(ET_5.data)), 
							rep("seedling", ncol(hls1_air_4.data)), rep("seedling", ncol(hls1_ET_4.data)),
							rep("seedling", ncol(ein2_air_1.data)), rep("seedling", ncol(ein2_ET_1.data)),
							rep("seedling", ncol(ein2_air_2.data)), rep("seedling", ncol(ein2_ET_2.data)),
							rep("seedling", ncol(ein2_air_4.data)), rep("seedling", ncol(ein2_ET_4.data)),
							rep("seedling", ncol(ein2_air_5.data)), rep("seedling", ncol(ein2_ET_5.data)),
							rep("seedling", ncol(ein2_air_6.data)), rep("seedling", ncol(ein2_ET_6.data)),
							rep("seedling", ncol(ctrl_6.data)), rep("seedling", ncol(ctrl_7.data)),
							rep("seedling", ncol(ctrl_8.data)), rep("seedling", ncol(ctrl_9.data)),
							rep("seedling", ncol(col_72ET_1.data)), rep("seedling", ncol(col_72ET_2.data)),
							rep("seedling", ncol(col_72ET_3.data)), rep("seedling", ncol(col_2ET_1.data)),
							rep("seedling", ncol(col_4ET_1.data)), rep("seedling", ncol(col_8ET_1.data)),
							rep("seedling", ncol(col_72ET_4.data)), rep("seedling", ncol(col_72ET_5.data)),
							rep("seedling", ncol(col_8ET_2.data)),
							rep("seedling", ncol(col_4ET_2.data))
							)

seedling@meta.data$cond <- c(rep("AIR_4", ncol(ctrl.data)), rep("ET_4", ncol(et_1.data)), rep("AIR_4", ncol(ctrl_2.data)), rep("ET_4", ncol(et_2.data)), rep("AIR_4", ncol(ctrl_3.data)), rep("ET_4", ncol(et_3.data)), rep("ET_0", ncol(ET0.data)), rep("ET_30m", ncol(ET30.data)), rep("ET_4", ncol(ET4.data)), rep("ET_24", ncol(ET24.data)), 
							rep("AIR_24_hook", ncol(air.data)), rep("ET_24_hook", ncol(et.data)), 
							rep("ET_30m", ncol(ET30_2.data)), rep("ET_4", ncol(ET4_2.data)), rep("ET_8", ncol(ET8.data)), rep("AIR_24", ncol(hls1_air.data)), rep("ET_24", ncol(hls1_ET.data)), rep("AIR_24", ncol(hls1_air_2.data)), rep("ET_24", ncol(hls1_ET_2.data)),
							rep("AIR_24", ncol(ctrl_4.data)), rep("ET_24", ncol(ET_4.data)), 
							rep("AIR_24", ncol(hls1_air_3.data)), rep("ET_24", ncol(hls1_ET_3.data)),
							rep("AIR_24", ncol(ctrl_5.data)), rep("ET_24", ncol(ET_5.data)), 
							rep("AIR_24", ncol(hls1_air_4.data)), rep("ET_24", ncol(hls1_ET_4.data)),
							rep("AIR_24", ncol(ein2_air_1.data)), rep("ET_24", ncol(ein2_ET_1.data)),
							rep("AIR_24", ncol(ein2_air_2.data)), rep("ET_24", ncol(ein2_ET_2.data)),
							rep("AIR_24", ncol(ein2_air_4.data)), rep("ET_24", ncol(ein2_ET_4.data)),
							rep("AIR_24", ncol(ein2_air_5.data)), rep("ET_24", ncol(ein2_ET_5.data)),
							rep("AIR_24", ncol(ein2_air_6.data)), rep("ET_24", ncol(ein2_ET_6.data)),
							rep("AIR_24", ncol(ctrl_6.data)), rep("AIR_24", ncol(ctrl_7.data)),
							rep("AIR_24", ncol(ctrl_8.data)), rep("AIR_24", ncol(ctrl_9.data)),
							rep("ET_72", ncol(col_72ET_1.data)), rep("ET_72", ncol(col_72ET_2.data)),
							rep("ET_72", ncol(col_72ET_3.data)), rep("ET_2", ncol(col_2ET_1.data)),
							rep("ET_4", ncol(col_4ET_1.data)), rep("ET_8", ncol(col_8ET_1.data)),
							rep("ET_72", ncol(col_72ET_4.data)), rep("ET_72", ncol(col_72ET_5.data)),
							rep("ET_8", ncol(col_8ET_2.data)),
							rep("ET_4", ncol(col_4ET_2.data))
							)

seedling@meta.data$geno <- c(rep("Col", ncol(ctrl.data)), rep("Col", ncol(et_1.data)), rep("Col", ncol(ctrl_2.data)), rep("Col", ncol(et_2.data)), rep("Col", ncol(ctrl_3.data)), rep("Col", ncol(et_3.data)), rep("Col", ncol(ET0.data)), rep("Col", ncol(ET30.data)), rep("Col", ncol(ET4.data)), rep("Col", ncol(ET24.data)), 
							rep("Col", ncol(air.data)), rep("Col", ncol(et.data)), 
							rep("Col", ncol(ET30_2.data)), rep("Col", ncol(ET4_2.data)), rep("Col", ncol(ET8.data)), rep("hls1", ncol(hls1_air.data)), rep("hls1", ncol(hls1_ET.data)), rep("hls1", ncol(hls1_air_2.data)), rep("hls1", ncol(hls1_ET_2.data)),
							rep("Col", ncol(ctrl_4.data)), rep("Col", ncol(ET_4.data)), 
							rep("hls1", ncol(hls1_air_3.data)), rep("hls1", ncol(hls1_ET_3.data)),
							rep("Col", ncol(ctrl_5.data)), rep("Col", ncol(ET_5.data)), 
							rep("hls1", ncol(hls1_air_4.data)), rep("hls1", ncol(hls1_ET_4.data)),
							rep("ein2", ncol(ein2_air_1.data)), rep("ein2", ncol(ein2_ET_1.data)),
							rep("ein2", ncol(ein2_air_2.data)), rep("ein2", ncol(ein2_ET_2.data)),
							rep("ein2", ncol(ein2_air_4.data)), rep("ein2", ncol(ein2_ET_4.data)),
							rep("ein2", ncol(ein2_air_5.data)), rep("ein2", ncol(ein2_ET_5.data)),
							rep("ein2", ncol(ein2_air_6.data)), rep("ein2", ncol(ein2_ET_6.data)),
							rep("Col", ncol(ctrl_6.data)), rep("Col", ncol(ctrl_7.data)),
							rep("Col", ncol(ctrl_8.data)), rep("Col", ncol(ctrl_9.data)),
							rep("Col", ncol(col_72ET_1.data)), rep("Col", ncol(col_72ET_2.data)),
							rep("Col", ncol(col_72ET_3.data)), rep("Col", ncol(col_2ET_1.data)),
							rep("Col", ncol(col_4ET_1.data)), rep("Col", ncol(col_8ET_1.data)),
							rep("Col", ncol(col_72ET_4.data)), rep("Col", ncol(col_72ET_5.data)),
							rep("Col", ncol(col_8ET_2.data)),
							rep("Col", ncol(col_4ET_2.data))
							)


seedling[["percent.mt"]] <- PercentageFeatureSet(object = seedling, pattern = "^^ATMG")
seedling[["percent.cp"]] <- PercentageFeatureSet(object = seedling, pattern = "^^ATCG")

seedling$percent.cp = round(seedling$percent.cp, digits = 2)
seedling$percent.mt = round(seedling$percent.mt, digits = 2)


seedling <- subset(x = seedling, subset = nFeature_RNA > 275 & nFeature_RNA < 5000 & nCount_RNA > 325 & percent.mt < 5 & percent.cp < 10)

seedling = seedling[row.names(seedling)[!grepl("ATCG", row.names(seedling))]]
seedling = seedling[row.names(seedling)[!grepl("ATMG", row.names(seedling))]]


seedling = Seurat::NormalizeData(seedling, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 10000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = seedling@var.genes, npcs = 50, verbose = FALSE)


seedling <- seedling %>% 
    RunHarmony(c("stim", "percent.cp", "tissue"), plot_convergence = FALSE, verbose=TRUE)


seedling <- seedling %>% 
    RunUMAP(reduction = "harmony", dims = 1:20, umap.method = 'umap-learn') %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()


derpy = c("#6b6b6b", "#00afed", "#294b99", "#9f3b5c", "#5f4dd4", "#edbfc4", "#e98737", "#dbcf21")
derpy_pal = data.frame(row.names = c(0:23), cols = colorRampPalette(derpy)(n = 24))

pdf("fig1C_seedling_UMAP.pdf")
DimPlot(seedling, label=T, cols = derpy_pal[,1])
DimPlot(seedling, label=T, cols = derpy_pal[,1], raster=TRUE)
dev.off()

pdf("results/harmony/harmony_umap.pdf")
DimPlot(object = seedling, reduction = "umap", group.by = "stim")
DimPlot(object = seedling, reduction = "umap", group.by = "tissue")
DimPlot(object = seedling, reduction = "umap", group.by = "cond")
DimPlot(object = seedling, reduction = "umap", group.by = "geno")
DimPlot(object = seedling, reduction = "umap", label = TRUE, split.by = "geno")

DimPlot(object = subset(seedling, geno == "Col"), reduction = "umap", label=T) + ggtitle("Col")
DimPlot(object = subset(seedling, geno == "hls1"), reduction = "umap", label=T) + ggtitle("hls1")
DimPlot(object = subset(seedling, geno == "ein2"), reduction = "umap", label=T) + ggtitle("ein2")
DimPlot(object = subset(seedling, tissue == "seedling"), reduction = "umap", label=T) + ggtitle("seedling")
DimPlot(object = subset(seedling, tissue == "hook"), reduction = "umap", label=T) + ggtitle("hook") + xlim(range(seedling@reductions$umap[[]][,1])) + ylim(range(seedling@reductions$umap[[]][,2]))
DimPlot(object = subset(seedling, stim == "CTRL"), reduction = "umap", label=T) + ggtitle("Air")
DimPlot(object = subset(seedling, stim == "ET"), reduction = "umap", label=T) + ggtitle("ET")
dev.off()


saveRDS(seedling, "results/harmony/harmony.rds")



###### Identify all cluster markers

cluster.markers <- FindAllMarkers(object = seedling, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

markers.df = cluster.markers

desc <- read.delim("/gale/netapp/home/trlee/general/Araport11_annotation.txt")[, c(1,2)]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
markers.df <- cbind(Desc=descv[gsub("_.*", "", markers.df$gene)], markers.df)

write.table(markers.df, "results/harmony/cluster_markers.xls", sep="\t", row.names=T, col.names=NA)


### top marker
top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


DefaultAssay(seedling) <- "RNA"
seedling <- NormalizeData(seedling, verbose = FALSE)
seedling = ScaleData(seedling)

pdf("results/harmony/dotplot_clusters.pdf")
DotPlot(seedling, features = unique(top1$gene), dot.scale = 8,
    dot.min=.15) + RotatedAxis()
dev.off()


saveRDS(seedling, "results/harmony/harmony.rds")



##### expression of literature marker genes
pdf("figS1_epidermis_markers.pdf")
epi = c("AT1G68530", "AT4G21750", "AT5G57800")
DotPlot(seedling, features=epi)
FeaturePlot(seedling, features=epi, order=T, min.cutoff="q5", max.cutoff="q95")
dev.off()

pdf("figS1_mesophyll_markers.pdf")
meso = c("AT1G29910", "AT5G38420", "AT2G34430", "AT3G27690")
DotPlot(seedling, features=meso)
FeaturePlot(seedling, features=meso, order=T, min.cutoff="q5", max.cutoff="q95")
dev.off()

pdf("figS1_xylem_markers.pdf")
xylem = c("AT1G71930", "AT1G62700", "AT4G35350", "AT1G20850")
FeaturePlot(seedling, features=xylem, order=T, min.cutoff="q5", max.cutoff="q95")
dev.off()

pdf("fig1_endodermis_markers.pdf")
endodermis = c("AT5G15290", "AT3G11550", "AT2G27370", "AT1G61590")
DotPlot(seedling, features=endodermis)
FeaturePlot(seedling, features=endodermis, order=T, min.cutoff="q5", max.cutoff="q95")
dev.off()

pdf("fig1_root_hair_markers.pdf")
hair = c("AT5G51060",  "AT3G51460", "AT5G19560", "AT3G03050")
DotPlot(seedling, features=hair)
FeaturePlot(seedling, features=hair, order=T, min.cutoff="q5", max.cutoff="q95")
dev.off()

pdf("fig1_meristem_markers.pdf")
meristem = c("AT4G33270",  "AT4G33260", "AT4G35620", "AT1G20610")
DotPlot(seedling, features=meristem)
FeaturePlot(seedling, features=meristem, order=T, min.cutoff="q5", max.cutoff="q95")
dev.off()


pdf("fig_S2_ET_induced_examples.pdf", width=12)
FeaturePlot(seedling, features=c("AT1G49570"), order=T, min.cutoff="q5", max.cutoff="q95", split.by="stim", cols=c("gray", "magenta"))
FeaturePlot(seedling, features=c("AT2G47160"), order=T, min.cutoff="q5", max.cutoff="q95", split.by="stim", cols=c("gray", "magenta"))
FeaturePlot(seedling, features=c("AT5G19890"), order=T, min.cutoff="q5", max.cutoff="q95", split.by="stim", cols=c("gray", "magenta"))
FeaturePlot(seedling, features=c("AT1G62380"), order=T, min.cutoff="q5", max.cutoff="q95", split.by="stim", cols=c("gray", "magenta"))
FeaturePlot(seedling, features=c("AT1G62380"), order=T, min.cutoff="q5", max.cutoff="q95", cols=c("gray", "magenta"))
dev.off()

pdf("fig_S3_roothair_ET_markers.pdf")
FeaturePlot(seedling, features=c("AT1G66470", "AT5G19560"), order=T, min.cutoff="q5", max.cutoff="q95", split.by="stim", cols=c("gray", "magenta"))
FeaturePlot(seedling, features=c("AT1G66470", "AT5G19560"), order=T, min.cutoff="q5", max.cutoff="q95", cols=c("gray", "magenta"))
dev.off()





###### quantify number of cells per sample

num_air = ncol(subset(seedling, stim=="CTRL"))
num_ET = ncol(subset(seedling, stim=="ET"))

num_Col = ncol(subset(seedling, geno=="Col"))
num_hls1 = ncol(subset(seedling, geno=="hls1"))
num_ein2 = ncol(subset(seedling, geno=="ein2"))

num_seedling = ncol(subset(seedling, tissue=="seedling"))
num_hook = ncol(subset(seedling, tissue=="hook"))

quantify_cells = as.data.frame(t(data.frame(num_air, num_ET, num_Col, num_hls1, num_ein2, num_seedling, num_hook)))
quantify_cells$sample = row.names(quantify_cells)

quantify_cells$sample <- factor(quantify_cells$sample, levels = quantify_cells$sample[order(quantify_cells$V1)])
quantify_cells$sample <- factor(quantify_cells$sample, levels = quantify_cells$sample)

quantify_cells$color = c("1", "2", "1", "2", "3", "1", "2")

pdf("results/harmony/cell_distribution.pdf")
p = ggplot(data=quantify_cells, aes(x = sample, y = V1, fill = color))
p + geom_bar(stat="identity") +
	ylab("number of cells") +
	theme_classic() +
	theme(legend.position = "none", axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dev.off()


##### number of cells per cluster

clust_num = length(levels(seedling@meta.data$seurat_clusters))

clusters = as.numeric(as.character(seedling@meta.data$seurat_clusters))
clusters = cbind(as.numeric(clusters), seedling@meta.data$tissue, seedling@meta.data$stim, seedling@meta.data$geno)

cluster_cells = data.frame()
seedling.cluster_cells = data.frame()
hook.cluster_cells = data.frame()
air.cluster_cells = data.frame()
ET.cluster_cells = data.frame()
hls.cluster_cells = data.frame()
col.cluster_cells = data.frame()
ein2.cluster_cells = data.frame()


for(i in seq(1:clust_num))
{
num_cells = data.frame(i-1, length(clusters[,1][clusters == (i-1)]))
seedling.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,2] == "seedling"]))
hook.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,2] == "hook"]))
air.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,3] == "CTRL"]))
ET.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,3] == "ET"]))
col.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,4] == "Col"]))
hls.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,4] == "hls1"]))
ein2.num_cells = data.frame(i-1, length(clusters[clusters == (i-1) & clusters[,4] == "ein2"]))

cluster_cells = rbind(cluster_cells, num_cells)

seedling.cluster_cells = rbind(seedling.cluster_cells, seedling.num_cells)
hook.cluster_cells = rbind(hook.cluster_cells, hook.num_cells)
air.cluster_cells = rbind(air.cluster_cells, air.num_cells)
ET.cluster_cells = rbind(ET.cluster_cells, ET.num_cells)
col.cluster_cells = rbind(col.cluster_cells, col.num_cells)
hls.cluster_cells = rbind(hls.cluster_cells, hls.num_cells)
ein2.cluster_cells = rbind(ein2.cluster_cells, ein2.num_cells)

}



colnames(cluster_cells) = c("Cluster_number", "Cell_number")
row.names(cluster_cells) = levels(Idents(seedling))
cluster_cells$Cluster_ident = levels(Idents(seedling))
cluster_cells$Cluster_ident <- factor(cluster_cells$Cluster_ident, levels = cluster_cells$Cluster_ident)


cluster_cells$seedling = seedling.cluster_cells[,2]
cluster_cells$hook = hook.cluster_cells[,2]
cluster_cells$air = air.cluster_cells[,2]
cluster_cells$et = ET.cluster_cells[,2]
cluster_cells$col = col.cluster_cells[,2]
cluster_cells$hls = hls.cluster_cells[,2]
cluster_cells$ein2 = ein2.cluster_cells[,2]


###### Plot number of cells per cluster


cluster_cells$percent = cluster_cells$Cell_number / sum(cluster_cells$Cell_number)
cluster_cells$rep = 1

cluster_cells$percent_seedling = cluster_cells$seedling / cluster_cells$Cell_number
cluster_cells$percent_hook = cluster_cells$hook / cluster_cells$Cell_number

cluster_cells$proportion_seedling = cluster_cells$seedling / sum(cluster_cells$seedling)
cluster_cells$proportion_hook = cluster_cells$hook / sum(cluster_cells$hook)

subset = cluster_cells[,c("Cluster_number","proportion_hook", "proportion_seedling")]

df.melt = melt(subset, id.vars=c("Cluster_number"))

pdf("results/harmony/Cluster_proportion_enrichment_cells.pdf")

p = ggplot(data=df.melt, aes(fill=variable))
p + geom_bar(aes(x=Cluster_number, y=value, fill=variable), stat = "identity", width=.8, color="black", size=.3) +
	ylab("Cluster percentage") +
	xlab("Cluster number") +
	guides(fill=guide_legend(title="Cluster")) +
	theme_classic()
	#theme(aspect.ratio = 3/1)

dev.off()


##### Proportion of air vs. ET treated cells

cluster_cells$percent_air = cluster_cells$air / cluster_cells$Cell_number
cluster_cells$percent_ET = cluster_cells$et / cluster_cells$Cell_number

subset = cluster_cells[,c("Cluster_number","percent_air","percent_ET")]
df.melt = melt(subset, id.vars=c("Cluster_number"))

pdf("results/harmony/Cluster_proportion_treatment.pdf")

p = ggplot(data=df.melt, aes(fill=variable))
p + geom_bar(aes(x=Cluster_number, y=value, fill=variable), stat = "identity", width=.8, color="black", size=.3) +
	ylab("Cluster percentage") +
	xlab("Cluster number") +
	guides(fill=guide_legend(title="Cluster")) +
	theme_classic()
dev.off()


##### Proportion of cells by genotype

cluster_cells$percent_col = cluster_cells$col / cluster_cells$Cell_number
cluster_cells$percent_hls = cluster_cells$hls / cluster_cells$Cell_number
cluster_cells$percent_ein2 = cluster_cells$ein2 / cluster_cells$Cell_number

subset = cluster_cells[,c("Cluster_number","percent_col","percent_hls","percent_ein2")]
df.melt = melt(subset, id.vars=c("Cluster_number"))

pdf("results/harmony/Cluster_proportion_genotype.pdf")

p = ggplot(data=df.melt, aes(fill=variable))
p + geom_bar(aes(x=Cluster_number, y=value, fill=variable), stat = "identity", width=.8, color="black", size=.3) +
	ylab("Cluster percentage") +
	xlab("Cluster number") +
	guides(fill=guide_legend(title="Cluster")) +
	theme_classic()
	#theme(aspect.ratio = 3/1)

dev.off()


cluster_cells$proportion_air = cluster_cells$air / sum(cluster_cells$air)
cluster_cells$proportion_ET = cluster_cells$et / sum(cluster_cells$et)

cluster_cells$proportion_col = cluster_cells$col / sum(cluster_cells$col)
cluster_cells$proportion_hls = cluster_cells$hls / sum(cluster_cells$hls)
cluster_cells$proportion_ein2 = cluster_cells$ein2 / sum(cluster_cells$ein2)

cluster_cells$proportion_col_norm = (1 / (cluster_cells$proportion_col + cluster_cells$proportion_hls + cluster_cells$proportion_ein2)) * cluster_cells$col / sum(cluster_cells$col)
cluster_cells$proportion_hls_norm = (1 / (cluster_cells$proportion_col + cluster_cells$proportion_hls + cluster_cells$proportion_ein2)) * cluster_cells$hls / sum(cluster_cells$hls)
cluster_cells$proportion_ein2_norm = (1 / (cluster_cells$proportion_col + cluster_cells$proportion_hls + cluster_cells$proportion_ein2)) * cluster_cells$ein2 / sum(cluster_cells$ein2)

subset = cluster_cells[,c("Cluster_number","proportion_col_norm","proportion_hls_norm","proportion_ein2_norm")]
df.melt = melt(subset, id.vars=c("Cluster_number"))

pdf("results/harmony/Cluster_relative_percentage_cells_genotype.pdf")

p = ggplot()
p + geom_bar(data = df.melt, aes(x=Cluster_number, y=value, fill=variable), stat = "identity", width=.8, color="black", size=.3) +
	geom_rect(aes(xmin = min(df.melt$Cluster_number) - 0.5, ymin = 0.33, xmax= max(df.melt$Cluster_number) + 0.5, ymax = 0.66, alpha=0.1), fill="gray", color="blue") +
	ylab("Cluster percentage") +
	xlab("Cluster number") +
	guides(fill=guide_legend(title="Cluster")) +
	theme_classic()
dev.off()

