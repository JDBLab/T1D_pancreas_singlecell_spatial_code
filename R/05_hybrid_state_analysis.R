# ============================================================
# 05_hybrid_state_analysis.R
# Generalized analysis of hybrid/progenitor-like states across conditions.
# Does not include low-threshold deconvolution or sample-specific rescue steps.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

source("R/functions_scRNAseq.R")

objdir <- "objects/annotated"
outdir <- "results/hybrid_state_analysis"
safe_dir_create(outdir)
safe_dir_create(file.path(outdir, "plots"))
safe_dir_create(file.path(outdir, "tables"))
safe_dir_create("objects/hybrid")

# Edit this table to match the final annotations used in the manuscript.
# pattern is matched to the Celltype_renamed column.
hybrid_query <- tibble::tribble(
  ~condition, ~pattern, ~hybrid_family,
  "T0",    "^Hy1($|_)", "Hy1",
  "CTD5",  "^Hy2($|_)", "Hy2",
  "CTD10", "^Hy2($|_)", "Hy2",
  "BTD5",  "^Hy3($|_)", "Hy3",
  "BTD10", "^Hy4($|_)", "Hy4"
)

make_rna_only <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- DietSeurat(obj, assays = "RNA", dimreducs = NULL, graphs = NULL)
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
  obj
}

extract_hybrid <- function(condition, pattern, hybrid_family) {
  path <- file.path(objdir, paste0(condition, "_annotated.rds"))
  obj <- readRDS(path)
  stopifnot("Celltype_renamed" %in% colnames(obj@meta.data))
  labs <- as.character(obj$Celltype_renamed)
  keep <- grepl(pattern, labs)
  if (sum(keep) == 0) {
    warning("No cells found for ", hybrid_family, " in ", condition)
    return(NULL)
  }
  sub <- subset(obj, cells = Cells(obj)[keep])
  sub <- make_rna_only(sub)
  sub$HybridFamily <- hybrid_family
  sub$Condition <- condition
  sub$OriginalCelltype <- labs[keep]
  sub
}

hy_list <- pmap(hybrid_query, extract_hybrid)
hy_list <- hy_list[!sapply(hy_list, is.null)]

if (length(hy_list) < 2) stop("Need at least two hybrid groups for integrated analysis.")

Hybrid_all <- merge(hy_list[[1]], y = hy_list[-1], add.cell.ids = hybrid_query$condition[seq_along(hy_list)])
DefaultAssay(Hybrid_all) <- "RNA"
if (inherits(Hybrid_all[["RNA"]], "Assay5")) Hybrid_all <- JoinLayers(Hybrid_all, assay = "RNA")

Hybrid_all <- NormalizeData(Hybrid_all, verbose = FALSE)
Hybrid_all <- FindVariableFeatures(Hybrid_all, nfeatures = 3000, verbose = FALSE)
Hybrid_all <- ScaleData(Hybrid_all, verbose = FALSE)
Hybrid_all <- RunPCA(Hybrid_all, npcs = 30, verbose = FALSE)
Hybrid_all <- RunUMAP(Hybrid_all, dims = 1:20, verbose = FALSE)

Hybrid_all$HybridFamily <- factor(Hybrid_all$HybridFamily, levels = unique(hybrid_query$hybrid_family))

p1 <- DimPlot(Hybrid_all, reduction = "umap", group.by = "HybridFamily", label = TRUE, repel = TRUE) +
  ggtitle("Hybrid-state continuum")
ggsave(file.path(outdir, "plots", "Hybrid_UMAP_by_family.png"), p1, width = 7, height = 6, dpi = 300, bg = "white")

# General marker programs.
marker_programs <- list(
  Endocrine = c("CHGA", "SCG2", "SCG3", "PCSK1", "ISL1", "NEUROD1", "RFX6", "INSM1", "NKX6-1"),
  Ductal = c("KRT19", "KRT8", "KRT18", "SPP1", "MUC1", "KRT17", "TFF1", "TFF2"),
  Acinar = c("PTF1A", "CPA1", "CPA2", "PRSS1", "PRSS2", "CTRB1", "CTRB2"),
  Progenitor_Remodeling = c("SOX9", "SOX4", "KLF6", "JUN", "CLDN4", "CLU", "S100A6", "STMN1")
)

for (nm in names(marker_programs)) {
  feats <- intersect(marker_programs[[nm]], rownames(Hybrid_all))
  if (length(feats) >= 2) {
    Hybrid_all <- AddModuleScore(Hybrid_all, features = list(feats), name = paste0(nm, "Score"))
  }
}

score_cols <- grep("Score1$", colnames(Hybrid_all@meta.data), value = TRUE)
if (length(score_cols) > 0) {
  score_df <- Hybrid_all@meta.data %>%
    select(HybridFamily, all_of(score_cols)) %>%
    pivot_longer(cols = -HybridFamily, names_to = "Score", values_to = "Value")
  p_score <- ggplot(score_df, aes(x = HybridFamily, y = Value, fill = HybridFamily)) +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~Score, scales = "free_y", ncol = 2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Module score")
  ggsave(file.path(outdir, "plots", "Hybrid_module_scores.png"), p_score, width = 10, height = 6, dpi = 300, bg = "white")
}

panel_genes <- unique(unlist(marker_programs))
panel_genes <- intersect(panel_genes, rownames(Hybrid_all))
if (length(panel_genes) > 0) {
  p_dot <- DotPlot(Hybrid_all, features = panel_genes, group.by = "HybridFamily") + RotatedAxis() +
    ggtitle("Hybrid-state marker programs")
  ggsave(file.path(outdir, "plots", "Hybrid_DotPlot_marker_programs.png"), p_dot, width = 13, height = 5.5, dpi = 300, bg = "white")
}

Idents(Hybrid_all) <- "HybridFamily"
markers <- FindAllMarkers(Hybrid_all, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
write_csv(markers, file.path(outdir, "tables", "Hybrid_FindAllMarkers.csv"))
write_csv(Hybrid_all@meta.data, file.path(outdir, "tables", "Hybrid_metadata.csv"))
saveRDS(Hybrid_all, file.path("objects/hybrid", "Hybrid_all.rds"))
