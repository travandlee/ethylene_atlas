########
##### Imputation
########

library(Seurat)
library(SeuratWrappers)
library(dplyr)


system("ulimit -u 16384")
Sys.setenv(R_MAX_VSIZE = 8e10)


setwd("/ceph/ethylene/multiome/spatial_integration_250825")

merged = readRDS("spatial_multiome_integrated_250825.rds")

merged_RNA = CreateSeuratObject(merged[["RNA"]])
merged=NULL
gc()

merged_RNA[["RNA"]]$data.ethylene.1 = NULL
merged_RNA[["RNA"]]$scale.data.ethylene.1 = NULL
merged_RNA[["RNA"]]$data.1 = NULL
merged_RNA[["RNA"]]$scale.data.1 = NULL
merged_RNA[["RNA"]]$data.1.2 = NULL
merged_RNA[["RNA"]]$scale.data.1.2 = NULL
merged_RNA[["RNA"]]$data.2 = NULL
merged_RNA[["RNA"]]$scale.data.2 = NULL

merged_RNA = JoinLayers(merged_RNA)

merged_RNA = FindVariableFeatures(merged_RNA)
features = VariableFeatures(merged_RNA)
all_genes = row.names(merged_RNA[["RNA"]])

all_genes = {set.seed(3386); sample(all_genes)}

merged_RNA = NormalizeData(merged_RNA)
alra_mat = matrix(ncol=ncol(merged_RNA)) %>% as("sparseMatrix")
alra_mat = as(alra_mat, "sparseMatrix")


##### run imputation in chunks due to memory limitations
for(i in c(1:10))
{

assign(paste0("alra_subset_1_", i), matrix(ncol=ncol(merged_RNA)) %>% as("sparseMatrix"))


genes_start = 1+(1000*(i-1))
genes_end = (1000*i)
genes_use = all_genes[genes_start:genes_end]
assign("alra_subset", RunALRA(merged_RNA, k = 41, genes.use = genes_use))
assign(paste0("alra_mat_1_", i), alra_subset[["alra"]]$data)
rm(alra_subset)
gc(reset = TRUE)
}


alra_matrix_1 = rbind(alra_mat_1_1,
alra_mat_1_2,
alra_mat_1_3,
alra_mat_1_4,
alra_mat_1_5,
alra_mat_1_6,
alra_mat_1_7,
alra_mat_1_8,
alra_mat_1_9,
alra_mat_1_10
)

rm(alra_mat_1_1,
alra_mat_1_2,
alra_mat_1_3,
alra_mat_1_4,
alra_mat_1_5,
alra_mat_1_6,
alra_mat_1_7,
alra_mat_1_8,
alra_mat_1_9,
alra_mat_1_10
)



for(i in c(11:20))
{

genes_start = 1+(1000*(i-1))
genes_end = (1000*i)
genes_use = all_genes[genes_start:genes_end]
assign("alra_subset", RunALRA(merged_RNA, k = 41, genes.use = genes_use))
assign(paste0("alra_mat_1_", i), alra_subset[["alra"]]$data)
rm(alra_subset)
gc(reset = TRUE)
}

alra_matrix_2 = rbind(alra_mat_1_11,
alra_mat_1_12,
alra_mat_1_13,
alra_mat_1_14,
alra_mat_1_15,
alra_mat_1_16,
alra_mat_1_17,
alra_mat_1_18,
alra_mat_1_19,
alra_mat_1_20
)

rm(alra_mat_1_11,
alra_mat_1_12,
alra_mat_1_13,
alra_mat_1_14,
alra_mat_1_15,
alra_mat_1_16,
alra_mat_1_17,
alra_mat_1_18,
alra_mat_1_19,
alra_mat_1_20
)


alra_matrix = rbind(alra_matrix_1, alra_matrix_2)

writeMM(obj = alra_matrix, file = "alra_20k.mtx")


rm(alra_matrix_1, alra_matrix_2)

for(i in c(21:30))
{

genes_start = 1+(1000*(i-1))
genes_end = (1000*i)
genes_use = all_genes[genes_start:genes_end]
assign("alra_subset", RunALRA(merged_RNA, k = 41, genes.use = genes_use))
assign(paste0("alra_mat_1_", i), alra_subset[["alra"]]$data)
rm(alra_subset)
gc(reset = TRUE)
}


i=31

genes_start = 1+(1000*(i-1))
genes_end = (1000*i)
genes_use = all_genes[genes_start:genes_end]
assign("alra_subset", RunALRA(merged_RNA, k = 41, genes.use = genes_use))
assign(paste0("alra_mat_1_", i), alra_subset[["alra"]]$data)
rm(alra_subset)
gc(reset = TRUE)


i=32

genes_start = 1+(1000*(i-1))
genes_end = length(all_genes)
genes_use = all_genes[genes_start:genes_end]
assign("alra_subset", RunALRA(merged_RNA, k = 41, genes.use = genes_use))
assign(paste0("alra_mat_1_", i), alra_subset[["alra"]]$data)
rm(alra_subset)
gc(reset = TRUE)



alra_matrix_3 = rbind(alra_mat_1_21,
alra_mat_1_22,
alra_mat_1_23,
alra_mat_1_24,
alra_mat_1_25,
alra_mat_1_26,
alra_mat_1_27,
alra_mat_1_28,
alra_mat_1_29,
alra_mat_1_30,
alra_mat_1_31,
alra_mat_1_32
)

rm(alra_mat_1_21,
alra_mat_1_22,
alra_mat_1_23,
alra_mat_1_24,
alra_mat_1_25,
alra_mat_1_26,
alra_mat_1_27,
alra_mat_1_28,
alra_mat_1_29,
alra_mat_1_30,
alra_mat_1_31,
alra_mat_1_32
)

subset_alra = readRDS("alra_20k_features_250917.rds")
alra_matrix_1 = subset_alra[["alra"]]@data

alra = CreateAssayObject(data = rbind(alra_matrix_1, alra_matrix_3))

saveRDS(alra, "alra_all_features_250917.rds")

merged[["alra"]] <-  alra

saveRDS(merged, "spatial_multiome_imputed_250825.rds")
