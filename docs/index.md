# Analysis code for human diabetic pancreas single-cell and spatial transcriptomics

This GitHub Pages site summarizes the analysis workflows used in the accompanying manuscript.

## Workflows

### Longitudinal scRNA-seq

R scripts process 10x Genomics scRNA-seq matrices from human pancreatic slice cultures. The workflow includes QC, DoubletFinder, Seurat integration, PCA/UMAP, clustering, marker discovery, annotation, hybrid-state analysis, and export for scVelo.

### Visium HD

The Visium HD workflow supports 8-µm binned and segmentation-aware Space Ranger outputs, including spatial marker visualization, module scoring, and comparison of islet-associated versus duct-associated INS+ cells.

### Xenium

The Xenium workflow supports cell-segmented in situ transcriptomic analysis, ROI-based transition-region analysis, heatmaps, dot plots, and scRNA-seq-guided label transfer.

### stLearn spatial trajectory

The stLearn notebook performs ROI-level spatial pseudotime analysis, branch marker discovery, and ligand–receptor hotspot analysis in ducto-endocrine regions.

### SIRV, scVelo, and PAGA spatial velocity

The SIRV/scVelo notebook integrates scRNA-seq RNA velocity with Xenium spatial cells and summarizes inferred transition relationships using PAGA.

## Repository organization

- `R/`: R scripts for scRNA-seq and related analysis workflows
- `notebooks/`: Jupyter notebooks for stLearn and SIRV/scVelo/PAGA analyses
- `config/`: templates for sample sheets, annotations, and marker sets
- `environment/`: example environment files
- `figures/`: placeholder for exported figures

## Note

Raw sequencing data, large intermediate objects, and protected donor-level data are not included. Users should update local paths and sample sheets before running the workflows.
