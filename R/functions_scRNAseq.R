# ============================================================
# Shared functions for longitudinal scRNA-seq analysis
# Generalized from the manuscript workflow.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

safe_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_10x_sample <- function(sample_id,
                            condition,
                            data_dir,
                            donor_id = NA_character_,
                            min_cells = 3,
                            min_features = 200) {
  message("Reading sample: ", sample_id)
  x <- Read10X(data.dir = data_dir)
  if (is.list(x)) {
    # Keep gene expression matrix if Read10X returns multiple modalities.
    if ("Gene Expression" %in% names(x)) x <- x[["Gene Expression"]] else x <- x[[1]]
  }
  obj <- CreateSeuratObject(
    counts = x,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )
  obj$sample_id <- sample_id
  obj$condition <- condition
  obj$donor_id <- donor_id
  obj
}

add_qc_metrics <- function(obj, mito_pattern = "^MT-") {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)
  obj
}

filter_cells <- function(obj,
                         min_features = 200,
                         max_features = 8000,
                         max_percent_mt = 20) {
  subset(
    obj,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )
}

preprocess_single_sample <- function(obj,
                                     nfeatures = 3000,
                                     npcs = 30,
                                     dims = 1:20,
                                     resolution = 0.4,
                                     verbose = FALSE) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = verbose)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = verbose)
  obj <- ScaleData(obj, verbose = verbose)
  obj <- RunPCA(obj, npcs = npcs, verbose = verbose)
  obj <- FindNeighbors(obj, dims = dims, verbose = verbose)
  obj <- FindClusters(obj, resolution = resolution, verbose = verbose)
  obj <- RunUMAP(obj, dims = dims, verbose = verbose)
  obj
}

# DoubletFinder wrapper. Works with common DoubletFinder versions.
run_doubletfinder <- function(obj,
                              dims = 1:20,
                              expected_doublet_rate = 0.075,
                              pN = 0.25,
                              pK = NULL,
                              resolution_for_homotypic = 0.4,
                              verbose = TRUE) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    warning("DoubletFinder is not installed. Returning object without doublet calls.")
    obj$DoubletFinder_class <- "NotRun"
    return(obj)
  }

  # Preliminary processing if needed.
  if (!"pca" %in% Reductions(obj)) {
    obj <- preprocess_single_sample(obj, dims = dims, resolution = resolution_for_homotypic, verbose = FALSE)
  }
  Idents(obj) <- "seurat_clusters"

  n_cells <- ncol(obj)
  n_exp <- round(expected_doublet_rate * n_cells)

  # Homotypic adjustment using preliminary clusters.
  homotypic_prop <- tryCatch({
    DoubletFinder::modelHomotypic(obj$seurat_clusters)
  }, error = function(e) 0)
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))

  if (is.null(pK)) {
    if (verbose) message("Estimating pK using DoubletFinder parameter sweep...")
    sweep_res <- tryCatch({
      if ("paramSweep_v3" %in% getNamespaceExports("DoubletFinder")) {
        DoubletFinder::paramSweep_v3(obj, PCs = dims, sct = FALSE)
      } else {
        DoubletFinder::paramSweep(obj, PCs = dims, sct = FALSE)
      }
    }, error = function(e) NULL)

    if (!is.null(sweep_res)) {
      sweep_stats <- tryCatch(DoubletFinder::summarizeSweep(sweep_res, GT = FALSE), error = function(e) NULL)
      bcmvn <- tryCatch(DoubletFinder::find.pK(sweep_stats), error = function(e) NULL)
      if (!is.null(bcmvn) && nrow(bcmvn) > 0) {
        pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      }
    }
  }

  if (is.null(pK) || is.na(pK)) {
    warning("Could not estimate pK; using pK = 0.09. Please inspect DoubletFinder sweep results for final analyses.")
    pK <- 0.09
  }

  if (verbose) {
    message("DoubletFinder: n_cells=", n_cells,
            "; expected=", n_exp,
            "; homotypic_adjusted=", n_exp_adj,
            "; pK=", pK)
  }

  obj <- if ("doubletFinder_v3" %in% getNamespaceExports("DoubletFinder")) {
    DoubletFinder::doubletFinder_v3(obj, PCs = dims, pN = pN, pK = pK, nExp = n_exp_adj, reuse.pANN = FALSE, sct = FALSE)
  } else {
    DoubletFinder::doubletFinder(obj, PCs = dims, pN = pN, pK = pK, nExp = n_exp_adj, reuse.pANN = FALSE, sct = FALSE)
  }

  df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(df_col) > 0) {
    obj$DoubletFinder_class <- obj@meta.data[[tail(df_col, 1)]]
  } else {
    obj$DoubletFinder_class <- NA_character_
  }
  obj
}

remove_doublets <- function(obj, class_col = "DoubletFinder_class") {
  if (!class_col %in% colnames(obj@meta.data)) {
    warning(class_col, " not found. Returning original object.")
    return(obj)
  }
  subset(obj, subset = .data[[class_col]] != "Doublet")
}

integrate_condition <- function(obj_list,
                                condition_name,
                                nfeatures = 3000,
                                npcs = 30,
                                dims = 1:20,
                                resolution = 0.4,
                                outdir = "results") {
  safe_dir_create(file.path(outdir, "plots"))
  safe_dir_create(file.path(outdir, "tables"))

  message("Integrating condition: ", condition_name)
  obj_list <- lapply(obj_list, function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    x
  })

  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
  anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, dims = dims)
  combined <- IntegrateData(anchorset = anchors, dims = dims)

  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = npcs, verbose = FALSE)
  combined <- FindNeighbors(combined, dims = dims, verbose = FALSE)
  combined <- FindClusters(combined, resolution = resolution, verbose = FALSE)
  combined <- RunUMAP(combined, dims = dims, verbose = FALSE)

  p_sample <- DimPlot(combined, reduction = "umap", group.by = "sample_id") + ggtitle(paste0(condition_name, " by sample"))
  p_cluster <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle(paste0(condition_name, " clusters"))
  ggsave(file.path(outdir, "plots", paste0(condition_name, "_UMAP_sample.png")), p_sample, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "plots", paste0(condition_name, "_UMAP_clusters.png")), p_cluster, width = 7, height = 6, dpi = 300, bg = "white")

  DefaultAssay(combined) <- "RNA"
  combined <- NormalizeData(combined, verbose = FALSE)
  Idents(combined) <- "seurat_clusters"
  markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
  write_csv(markers, file.path(outdir, "tables", paste0(condition_name, "_FindAllMarkers.csv")))

  saveRDS(combined, file.path(outdir, paste0(condition_name, "_combined.rds")))
  combined
}

apply_cluster_annotation <- function(obj, mapping_csv, cluster_col = "seurat_clusters", output_col = "Celltype_renamed") {
  map <- read_csv(mapping_csv, show_col_types = FALSE)
  stopifnot(all(c("cluster", "celltype") %in% colnames(map)))
  map$cluster <- as.character(map$cluster)
  map_vec <- setNames(map$celltype, map$cluster)
  clusters <- as.character(obj@meta.data[[cluster_col]])
  obj[[output_col]] <- ifelse(clusters %in% names(map_vec), map_vec[clusters], paste0("Cluster_", clusters))
  obj[[output_col]] <- factor(obj[[output_col]][, 1])
  Idents(obj) <- output_col
  obj
}

get_assay_data_safe <- function(obj, assay = "RNA", layer = "data") {
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) {
    lyr <- Layers(obj[[assay]])
    use_layer <- if (layer %in% lyr) layer else if ("data" %in% lyr) "data" else lyr[1]
    GetAssayData(obj, assay = assay, layer = use_layer)
  } else {
    GetAssayData(obj, assay = assay, slot = layer)
  }
}
