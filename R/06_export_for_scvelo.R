# ============================================================
# 06_export_for_scvelo.R
# Export generalized Seurat objects for downstream scVelo/RNA velocity analysis.
# This exports counts.mtx, metadata.csv, pca.csv, umap.csv, and gene_names.csv.
# Spliced/unspliced loom files should be generated separately from BAM files using velocyto.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(dplyr)
})

source("R/functions_scRNAseq.R")

safe_dir_create("results/scvelo_export")

export_for_scvelo <- function(seurat_obj,
                              prefix,
                              out_dir = "results/scvelo_export",
                              assay = "RNA",
                              reduction_umap = NULL) {
  odir <- file.path(out_dir, prefix)
  safe_dir_create(odir)

  DefaultAssay(seurat_obj) <- assay
  if (inherits(seurat_obj[[assay]], "Assay5")) {
    seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  }

  seurat_obj$barcode <- colnames(seurat_obj)

  reductions <- Reductions(seurat_obj)
  if (is.null(reduction_umap)) {
    if ("umap" %in% reductions) reduction_umap <- "umap"
    else if ("umap_int" %in% reductions) reduction_umap <- "umap_int"
    else reduction_umap <- NA_character_
  }

  if (!is.na(reduction_umap)) {
    umap_emb <- Embeddings(seurat_obj, reduction = reduction_umap)
    write.csv(umap_emb, file = file.path(odir, "umap.csv"), quote = FALSE, row.names = TRUE)
    seurat_obj$UMAP_1 <- umap_emb[, 1]
    seurat_obj$UMAP_2 <- umap_emb[, 2]
  }

  meta_df <- seurat_obj@meta.data
  write.csv(meta_df, file = file.path(odir, "metadata.csv"), quote = FALSE, row.names = TRUE)

  counts_matrix <- if (inherits(seurat_obj[[assay]], "Assay5")) {
    GetAssayData(seurat_obj, assay = assay, layer = "counts")
  } else {
    GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }
  writeMM(counts_matrix, file = file.path(odir, "counts.mtx"))

  write.table(
    data.frame(gene = rownames(counts_matrix)),
    file = file.path(odir, "gene_names.csv"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  if ("pca" %in% reductions) {
    pca_emb <- Embeddings(seurat_obj, reduction = "pca")
    write.csv(pca_emb, file = file.path(odir, "pca.csv"), quote = FALSE, row.names = TRUE)
  }

  message("Exported scVelo input files to: ", normalizePath(odir))
}

# Example: export condition-level annotated objects.
# Edit object paths depending on which Seurat object is used as RNA velocity reference.
objects_to_export <- list.files("objects/annotated", pattern = "_annotated\\.rds$", full.names = TRUE)

for (path in objects_to_export) {
  obj <- readRDS(path)
  prefix <- sub("_annotated\\.rds$", "", basename(path))
  export_for_scvelo(obj, prefix = prefix)
}
