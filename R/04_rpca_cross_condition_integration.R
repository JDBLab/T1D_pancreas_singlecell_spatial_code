# ============================================================
# 04_rpca_cross_condition_integration.R
# Generalized reciprocal PCA integration for cross-condition comparisons.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

source("R/functions_scRNAseq.R")

objdir <- "objects/annotated"
outdir <- "results/rpca_integrations"
safe_dir_create(outdir)
safe_dir_create(file.path(outdir, "plots"))
safe_dir_create(file.path(outdir, "tables"))
safe_dir_create("objects/rpca_integrated")

# Define generalized comparisons. Edit labels to match your study.
comparisons <- list(
  Control_D5_vs_Treated_D5 = c("CTD5", "BTD5"),
  Control_D10_vs_Treated_D10 = c("CTD10", "BTD10"),
  Baseline_Control_Timecourse = c("T0", "CTD5", "CTD10"),
  Baseline_Treated_Timecourse = c("T0", "BTD5", "BTD10")
)

safe_join_layers <- function(x, assay = "RNA") {
  DefaultAssay(x) <- assay
  if (inherits(x[[assay]], "Assay5")) {
    x <- JoinLayers(x, assay = assay)
  }
  x
}

integrate_RPCA_general <- function(name, condition_names, nfeatures = 3000, npcs = 30, dims = 1:30, resolution = 0.4) {
  message("Running RPCA integration: ", name)

  objs <- lapply(condition_names, function(cond) {
    path <- file.path(objdir, paste0(cond, "_annotated.rds"))
    if (!file.exists(path)) path <- file.path("results/condition_integrations", paste0(cond, "_combined.rds"))
    obj <- readRDS(path)
    obj$dataset <- cond
    obj <- RenameCells(obj, add.cell.id = cond)
    obj <- safe_join_layers(obj, assay = "RNA")
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
    obj
  })
  names(objs) <- condition_names

  features <- SelectIntegrationFeatures(object.list = objs, nfeatures = nfeatures)
  objs <- lapply(objs, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, npcs = npcs, verbose = FALSE)
    x
  })

  anchors <- FindIntegrationAnchors(
    object.list = objs,
    anchor.features = features,
    reduction = "rpca",
    dims = dims
  )

  integ <- IntegrateData(anchorset = anchors, dims = dims)
  DefaultAssay(integ) <- "integrated"
  integ <- ScaleData(integ, verbose = FALSE)
  integ <- RunPCA(integ, npcs = npcs, verbose = FALSE)
  integ <- RunUMAP(integ, dims = dims, reduction.name = "umap_int", verbose = FALSE)
  integ <- FindNeighbors(integ, dims = dims, verbose = FALSE)
  integ <- FindClusters(integ, resolution = resolution, verbose = FALSE)

  p1 <- DimPlot(integ, reduction = "umap_int", group.by = "dataset") + ggtitle(paste0(name, " by dataset"))
  p2 <- DimPlot(integ, reduction = "umap_int", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle(paste0(name, " clusters"))
  ggsave(file.path(outdir, "plots", paste0(name, "_UMAP_dataset.png")), p1, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "plots", paste0(name, "_UMAP_clusters.png")), p2, width = 7, height = 6, dpi = 300, bg = "white")

  DefaultAssay(integ) <- "RNA"
  integ <- safe_join_layers(integ, assay = "RNA")
  Idents(integ) <- "seurat_clusters"
  markers <- FindAllMarkers(integ, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
  write_csv(markers, file.path(outdir, "tables", paste0(name, "_FindAllMarkers.csv")))

  saveRDS(integ, file.path("objects/rpca_integrated", paste0(name, "_integrated.rds")))
  invisible(integ)
}

for (nm in names(comparisons)) {
  integrate_RPCA_general(nm, comparisons[[nm]])
}
