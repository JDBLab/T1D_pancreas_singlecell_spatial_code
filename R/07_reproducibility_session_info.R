# ============================================================
# 07_reproducibility_session_info.R
# Save package versions for reproducibility.
# ============================================================

safe_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_dir_create("results/reproducibility")

sink("results/reproducibility/sessionInfo.txt")
print(sessionInfo())
sink()

writeLines("Saved sessionInfo to results/reproducibility/sessionInfo.txt")
