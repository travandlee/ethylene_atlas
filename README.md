# ethylene_atlas

This repository contains scripts used in the manuscript "Spatially Defined Transient Cell States Regulate Differential Growth in the Apical Hook"

[preprint link](biorxiv.org)

Browse the single-cell datasets at our web resource [https://arabidopsisdevatlas.salk.edu/](https://arabidopsisdevatlas.salk.edu/)

Processed files can be found at [GEO accession # XXXXXX](https://www.ncbi.nlm.nih.gov/)

<img width="1650" height="1050" alt="multimodal_integration" src="https://github.com/user-attachments/assets/731f338e-8148-41e1-b7a7-a59f2fdabb03" />

![spatial_clusters_resize](https://github.com/user-attachments/assets/5c0415de-830f-408e-95c6-065df2a159fc)


## Arabidopsis seedling analysis

Scripts used to preprocess, cluster, and analyze the snRNA-seq datasets of whole seedlings are included below.

[01_seedling_preprocessing_basic_analysis](scripts/01_seedling_preprocessing_basic_analysis.r)

[02_ethylene_analysis.r](scripts/02_ethylene_analysis.r)

[03_apical_seedling_analysis.r](scripts/03_apical_seedling_analysis.r)

## Spatial analysis

Scripts used to preprocess, cluster, and analyze the MERFISH spatial transcripomic datasets are included below.

Baysor was used to perform cell-segmentation. Processed cell segmented polygon files can be found at our [web resource](https://neomorph.salk.edu/tlee_public/apical_hook/)

Spatial analyses requires Seurat v5 for spatial dataset handling.

[04_spatial_preprocessing.r](scripts/04_spatial_preprocessing.r)

[05_spatial_analysis.r](scripts/05_spatial_analysis.r)

## Multiome analysis

Scripts used to preprocess and prepare the multiome datasets for integration are included below.

[06_multiome_preprocessing.r](scripts/06_multiome_preprocessing.r)

[07_multiome_integration.r](scripts/07_multiome_integration.r)

## Spatial multi-modal analyses

Scripts used to integrate the spatial and single-nucleus transcriptome and multiome datasets are included below.

Spatial analyses requires Seurat v5 for spatial dataset handling.

[08_spatial_integration_analysis.r](scripts/08_spatial_integration_analysis.r)

[09_expression_imputation.r](scripts/09_expression_imputation.r)

## AHR cell population analysis

Scripts used to analyze the AHR cells and construct the multiomic gene regulatory network are included below.

[10_AHR_subset_analysis.r](scripts/10_AHR_subset_analysis.r)

[11_scMEGA_GRN_analysis.r](scripts/11_scMEGA_GRN_analysis.r)
