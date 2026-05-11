# ============================================================
# 02_integrate_by_condition.R
# Integrate replicate libraries within each condition and generate UMAPs/markers.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

source("R/functions_scRNAseq.R")

sample_sheet <- read_csv("config/sample_sheet_template.csv", show_col_types = FALSE)
qc_objdir <- "objects/qc"
outdir <- "results/condition_integrations"
safe_dir_create(outdir)
safe_dir_create(file.path(outdir, "plots"))
safe_dir_create(file.path(outdir, "tables"))

conditions <- unique(sample_sheet$condition)
condition_objects <- list()

for (cond in conditions) {
  ss <- sample_sheet %>% filter(condition == cond)
  obj_list <- lapply(ss$sample_id, function(sid) {
    readRDS(file.path(qc_objdir, paste0(sid, "_qc_doublet_removed.rds")))
  })
  names(obj_list) <- ss$sample_id

  combined <- integrate_condition(
    obj_list = obj_list,
    condition_name = cond,
    nfeatures = 3000,
    npcs = 30,
    dims = 1:20,
    resolution = 0.4,
    outdir = outdir
  )
  condition_objects[[cond]] <- combined
}

# Summary of final cells per condition.
cell_summary <- tibble(
  condition = names(condition_objects),
  n_cells = sapply(condition_objects, ncol)
)
write_csv(cell_summary, file.path(outdir, "cell_counts_by_condition.csv"))
print(cell_summary)
