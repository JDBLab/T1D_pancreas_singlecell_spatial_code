# ============================================================
# run_all.R
# Driver script for the generalized scRNA-seq workflow.
# Run from the repository root.
# ============================================================

source("R/01_qc_doubletfinder.R")
source("R/02_integrate_by_condition.R")
source("R/03_cluster_annotation_and_markers.R")
source("R/04_rpca_cross_condition_integration.R")
source("R/05_hybrid_state_analysis.R")
source("R/06_export_for_scvelo.R")
source("R/07_reproducibility_session_info.R")
