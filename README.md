# Human diabetic pancreas single-cell and spatial transcriptomics analysis

This repository contains generalized analysis code and notebook workflows accompanying a manuscript on epithelial remodeling, ductal–endocrine transitional states, and spatial cell-state organization in human diabetic pancreatic tissue.

The repository is organized as a reusable analysis template. Local file paths, raw data files, large intermediate objects, and protected donor-level data are not included. Users should update the sample sheets and input paths to match their own data organization.

## Analysis overview

The study integrates longitudinal single-cell RNA-seq, Visium HD spatial transcriptomics, Xenium in situ transcriptomics, spatial trajectory inference, ligand–receptor analysis, RNA velocity, and spatial PAGA-based transition mapping.

Major analysis modules include:

1. **Longitudinal scRNA-seq analysis** of human T1D pancreatic slice cultures across baseline, control, and BMP7-treated conditions.
2. **Visium HD spatial transcriptomics**, including 8-µm binned outputs, segmentation-aware analysis, module scoring, and analysis of islet-associated and duct-associated INS+ cells.
3. **Xenium in situ spatial transcriptomics**, including cell segmentation, ROI analysis, marker heatmaps, dot plots, and scRNA-seq-based label transfer.
4. **stLearn spatial trajectory and ligand–receptor analysis** of ducto-endocrine Xenium regions.
5. **RNA velocity, scVelo, SIRV, and PAGA analysis** to infer candidate transcriptional transitions and project scRNA-seq dynamics into spatial tissue context.

## Repository contents

```text
R/
  Generalized R scripts for scRNA-seq preprocessing, QC, DoubletFinder, integration,
  annotation, marker discovery, hybrid-state analysis, and scVelo export.

notebooks/
  Jupyter notebook templates for stLearn spatial trajectory/CCI analysis and
  SIRV/scVelo/PAGA spatial RNA velocity analysis. Replace placeholders with local paths.

config/
  Example sample sheets, cluster annotation templates, and marker gene set tables.

docs/
  GitHub Pages landing page.

figures/
  Placeholder folder for exported manuscript-style plots. Raw or large images should not be committed.

environment/
  Example dependency files for R and Python environments.
```

## 1. Longitudinal scRNA-seq workflow

The scRNA-seq workflow begins with 10x Genomics filtered feature-barcode matrices and performs:

- import of 10x matrices using Seurat
- Seurat object creation
- sample metadata assignment
- mitochondrial RNA fraction calculation
- quality-control filtering
- DoubletFinder-based doublet detection
- replicate-level integration within each experimental condition
- PCA, UMAP, graph-based clustering, and marker discovery
- curated cell-state annotation
- hybrid-state extraction and comparative marker analysis
- export of counts, metadata, PCA, UMAP, and gene names for downstream scVelo analysis

The workflow is designed for longitudinal pancreatic slice datasets with conditions such as baseline, untreated/control culture, and BMP7-treated culture, but the scripts can be adapted to other experimental designs by editing the sample sheet.

## 2. Visium HD workflow

The Visium HD workflow supports analysis of Space Ranger outputs, including:

- 8-µm binned spatial expression matrices
- segmentation-aware `filtered_feature_cell_matrix.h5` outputs
- barcode-to-segmentation mapping files
- spatial marker visualization
- module scoring for endocrine, beta-cell, ductal, acinar, progenitor, and maturity programs
- comparison of islet-associated and duct-associated INS+ cells
- scRNA-seq-guided spatial annotation and deconvolution

This repository provides generalized code structure and marker templates. Raw Space Ranger outputs are not included.

## 3. Xenium workflow

The Xenium workflow supports analysis of cell-segmented in situ transcriptomic outputs, including:

- import of Xenium cell-feature matrices and spatial coordinates
- QC, normalization, dimensionality reduction, clustering, and spatial visualization
- import of Xenium Analyzer cluster annotations
- ROI-based analysis using Xenium Explorer exports
- heatmap and dot plot analysis of curated marker programs
- scRNA-seq-to-Xenium label transfer
- prediction-score filtering and marker-guided refinement of selected epithelial states

Custom or study-specific cell-state labels should be updated in the annotation template rather than hard-coded into the scripts.

## 4. stLearn spatial trajectory notebook

The `notebooks/stLearn_spatial_trajectory_CCI.ipynb` notebook provides a structured workflow for:

- loading Xenium expression data and spatial coordinates
- importing ROI cell lists or polygons from Xenium Explorer
- Louvain clustering within focused ducto-endocrine regions
- mapping original Xenium clusters to ROI-specific Louvain states
- stLearn pseudotime spatial trajectory analysis
- branch-specific transition-marker discovery
- ligand–receptor hotspot analysis using a curated human ligand–receptor database
- exporting publication-ready plots and tables

## 5. SIRV, scVelo, and PAGA notebook

The `notebooks/SIRV_scVelo_PAGA_spatial_velocity.ipynb` notebook provides a structured workflow for:

- loading a processed scRNA-seq reference
- integrating velocyto loom files containing spliced and unspliced RNA layers
- constructing an AnnData object with expression, metadata, and embeddings
- running scVelo dynamical RNA velocity analysis
- computing velocity confidence, velocity length, and velocity pseudotime
- applying PAGA to summarize velocity-informed state relationships
- loading Xenium spatial cells and coordinates
- projecting scRNA-seq velocity information onto spatial cells using SIRV-style integration
- overlaying inferred cell states and transition arrows on Xenium coordinates or matched H&E images

## Data availability

This repository contains analysis code and notebook templates only. Raw FASTQ files, BAM files, loom files, Xenium output folders, Space Ranger output folders, Seurat RDS objects, AnnData H5AD files, and large image files are not included.

## How to run

1. Clone or download this repository.
2. Install the required R and Python packages.
3. Edit `config/sample_sheet_template.csv` to point to local 10x Genomics filtered matrix folders.
4. Edit notebook path variables such as `DATA_DIR`, `OUT_DIR`, and `PROJECT_DIR`.
5. Run the R scripts in numerical order or use `R/run_all.R`.
6. Run the notebooks after generating the required intermediate objects.

## Session Information

sessionInfo()
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] purrr_1.1.0        tibble_3.3.0       stringr_1.5.1      tidyr_1.3.1        readr_2.1.5        Matrix_1.7-3       patchwork_1.3.2    ggplot2_3.5.2     
 [9] dplyr_1.1.4        Seurat_5.4.0       SeuratObject_5.4.0 sp_2.2-0          

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0         magrittr_2.0.3         ggbeeswarm_0.7.2       spatstat.utils_3.1-5  
  [7] farver_2.1.2           ragg_1.5.0             fs_1.6.6               vctrs_0.6.5            ROCR_1.0-11            memoise_2.0.1         
 [13] spatstat.explore_3.5-2 rstatix_0.7.2          htmltools_0.5.8.1      forcats_1.0.1          usethis_3.2.1          broom_1.0.10          
 [19] cellranger_1.1.0       Formula_1.2-5          sctransform_0.4.2      parallelly_1.45.1      KernSmooth_2.23-24     htmlwidgets_1.6.4     
 [25] ica_1.0-3              plyr_1.8.9             plotly_4.11.0          zoo_1.8-14             cachem_1.1.0           igraph_2.1.4          
 [31] mime_0.13              lifecycle_1.0.4        pkgconfig_2.0.3        R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-4    
 [37] future_1.67.0          shiny_1.11.1           digest_0.6.37          tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.5.1         
 [43] pkgload_1.4.1          textshaping_1.0.3      ggpubr_0.6.1           labeling_0.4.3         progressr_0.17.0       spatstat.sparse_3.1-0 
 [49] httr_1.4.7             polyclip_1.10-7        abind_1.4-8            compiler_4.4.1         remotes_2.5.0          bit64_4.6.0-1         
 [55] withr_3.0.2            backports_1.5.0        carData_3.0-5          fastDummies_1.7.5      pkgbuild_1.4.8         ggsignif_0.6.4        
 [61] MASS_7.3-60.2          sessioninfo_1.2.3      tools_4.4.1            vipor_0.4.7            lmtest_0.9-40          beeswarm_0.4.0        
 [67] httpuv_1.6.16          future.apply_1.20.0    goftest_1.2-3          glue_1.8.0             nlme_3.1-164           promises_1.3.3        
 [73] grid_4.4.1             Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.4         generics_0.1.4         hdf5r_1.3.12          
 [79] gtable_0.3.6           spatstat.data_3.1-8    tzdb_0.5.0             hms_1.1.3              data.table_1.17.8      utf8_1.2.6            
 [85] car_3.1-3              spatstat.geom_3.5-0    RcppAnnoy_0.0.22       ggrepel_0.9.6          RANN_2.6.2             pillar_1.11.1         
 [91] vroom_1.6.5            spam_2.11-1            RcppHNSW_0.6.0         later_1.4.4            splines_4.4.1          lattice_0.22-6        
 [97] bit_4.6.0              survival_3.6-4         deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
[103] gridExtra_2.3          scattermore_1.2        devtools_2.4.6         matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.2        
[109] codetools_0.2-20       BiocManager_1.30.26    cli_3.6.5              uwot_0.2.3             systemfonts_1.2.3      xtable_1.8-4          
[115] reticulate_1.43.0      dichromat_2.0-0.1      Rcpp_1.1.0             readxl_1.4.5           globals_0.18.0         spatstat.random_3.4-1 
[121] png_0.1-8              ggrastr_1.0.2          spatstat.univar_3.1-4  parallel_4.4.1         ellipsis_0.3.2         dotCall64_1.2         
[127] listenv_0.9.1          viridisLite_0.4.2      scales_1.4.0           ggridges_0.5.7         crayon_1.5.3           writexl_1.5.4         
[133] rlang_1.1.6            cowplot_1.2.0

## Reproducibility

The scripts are provided as generalized, reusable workflows. Exact thresholds, marker sets, and sample labels should be reviewed and adjusted for each dataset. Environment and package versions should be recorded at the time of execution using the included session-information script and environment templates.

## Citation

If you use this repository, please cite the associated manuscript once available.
