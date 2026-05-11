# ============================================================
# 03_cluster_annotation_and_markers.R
# Apply cluster annotations from a user-provided mapping table and plot markers.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

source("R/functions_scRNAseq.R")

indir <- "results/condition_integrations"
outdir <- "results/annotation"
safe_dir_create(outdir)
safe_dir_create(file.path(outdir, "plots"))
safe_dir_create(file.path(outdir, "tables"))
safe_dir_create("objects/annotated")

# The mapping file should contain columns:
# condition,cluster,celltype
# Example:
# T0,0,Acinar
# T0,1,Ductal
mapping_file <- "config/cluster_annotation_template.csv"

if (!file.exists(mapping_file)) {
  message("No cluster annotation file found at ", mapping_file)
  message("Creating a template. Fill it and rerun this script.")
  template <- tibble(condition = character(), cluster = character(), celltype = character())
  write_csv(template, mapping_file)
  quit(save = "no")
}

map_all <- read_csv(mapping_file, show_col_types = FALSE)
conditions <- unique(map_all$condition)

# Marker panels can be adjusted for the final manuscript.
marker_panels <- list(
  Acinar = c("PRSS1", "PRSS2", "CPA1", "CPB1", "CTRB1", "REG1A"),
  Ductal = c("KRT7", "KRT8", "KRT18", "KRT19", "SOX9", "SPP1", "TFF1", "TFF2"),
  Endocrine = c("INS", "IAPP", "CHGA", "CHGB", "SCG3", "PCSK1", "PDX1", "NKX6-1", "MAFA"),
  Stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "VIM"),
  Immune = c("PTPRC", "CD3D", "CD3E", "CD14", "LYZ", "MS4A1")
)

for (cond in conditions) {
  obj_path <- file.path(indir, paste0(cond, "_combined.rds"))
  if (!file.exists(obj_path)) next

  obj <- readRDS(obj_path)
  cond_map <- map_all %>% filter(condition == cond) %>% select(cluster, celltype)
  tmp_map <- file.path(outdir, "tables", paste0(cond, "_cluster_annotation_map.csv"))
  write_csv(cond_map, tmp_map)
  obj <- apply_cluster_annotation(obj, tmp_map, output_col = "Celltype_renamed")

  # UMAP with annotations.
  p <- DimPlot(obj, reduction = "umap", group.by = "Celltype_renamed", label = TRUE, repel = TRUE, raster = FALSE) +
    ggtitle(paste0(cond, " annotated cell states"))
  ggsave(file.path(outdir, "plots", paste0(cond, "_UMAP_annotated.png")), p, width = 8, height = 6, dpi = 300, bg = "white")

  # DotPlot for lineage markers.
  genes <- unique(unlist(marker_panels))
  genes <- intersect(genes, rownames(obj))
  if (length(genes) > 0) {
    p_dot <- DotPlot(obj, features = genes, group.by = "Celltype_renamed") + RotatedAxis() +
      ggtitle(paste0(cond, " lineage marker dot plot"))
    ggsave(file.path(outdir, "plots", paste0(cond, "_DotPlot_lineage_markers.png")), p_dot, width = 14, height = 6, dpi = 300, bg = "white")
  }

  Idents(obj) <- "Celltype_renamed"
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
  write_csv(markers, file.path(outdir, "tables", paste0(cond, "_FindAllMarkers_by_celltype.csv")))
  saveRDS(obj, file.path("objects/annotated", paste0(cond, "_annotated.rds")))
}
