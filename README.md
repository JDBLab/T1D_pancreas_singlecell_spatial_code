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

## Reproducibility

The scripts are provided as generalized, reusable workflows. Exact thresholds, marker sets, and sample labels should be reviewed and adjusted for each dataset. Environment and package versions should be recorded at the time of execution using the included session-information script and environment templates.

## Citation

If you use this repository, please cite the associated manuscript once available.
