# ============================================================
# 01_qc_doubletfinder.R
# Import 10x matrices, perform QC, optional DoubletFinder, and save cleaned objects.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(DoubletFinder)
})

source("R/functions_scRNAseq.R")

sample_sheet <- read_csv("config/sample_sheet_template.csv", show_col_types = FALSE)
outdir <- "results/qc"
objdir <- "objects/qc"
safe_dir_create(outdir)
safe_dir_create(objdir)

# General QC settings. Adjust after inspecting each library.
qc_min_features <- 200
qc_max_features <- 8000
qc_max_percent_mt <- 20
expected_doublet_rate <- 0.075

qc_summary <- list()

for (i in seq_len(nrow(sample_sheet))) {
  row <- sample_sheet[i, ]

  obj <- read_10x_sample(
    sample_id = row$sample_id,
    condition = row$condition,
    data_dir = row$data_dir,
    donor_id = row$donor_id
  )

  obj <- add_qc_metrics(obj)

  n_initial <- ncol(obj)

  p_qc <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
  ggsave(file.path(outdir, paste0(row$sample_id, "_QC_violin_before_filtering.png")), p_qc, width = 10, height = 4, dpi = 300, bg = "white")

  obj <- filter_cells(
    obj,
    min_features = qc_min_features,
    max_features = qc_max_features,
    max_percent_mt = qc_max_percent_mt
  )
  n_after_qc <- ncol(obj)

  # Preliminary processing for DoubletFinder.
  obj <- preprocess_single_sample(obj, dims = 1:20, resolution = 0.4, verbose = FALSE)
  obj <- run_doubletfinder(
    obj,
    dims = 1:20,
    expected_doublet_rate = expected_doublet_rate,
    resolution_for_homotypic = 0.4
  )
  n_doublets <- sum(obj$DoubletFinder_class == "Doublet", na.rm = TRUE)
  obj <- remove_doublets(obj)
  n_final <- ncol(obj)

  saveRDS(obj, file.path(objdir, paste0(row$sample_id, "_qc_doublet_removed.rds")))

  qc_summary[[row$sample_id]] <- tibble(
    sample_id = row$sample_id,
    condition = row$condition,
    donor_id = row$donor_id,
    n_initial = n_initial,
    n_after_qc = n_after_qc,
    n_doublets = n_doublets,
    n_final = n_final
  )
}

qc_summary <- bind_rows(qc_summary)
write_csv(qc_summary, file.path(outdir, "QC_summary_by_sample.csv"))
print(qc_summary)
