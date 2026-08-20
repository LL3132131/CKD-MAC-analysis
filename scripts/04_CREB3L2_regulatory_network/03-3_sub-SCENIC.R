# =============================================================================
# Title: 03-3_sub-SCENIC.R
# Purpose: Export the SMC-enriched subset expression matrix for SCENIC
# Manuscript Figure/Table: Figure 4G and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: Final SMC Seurat object
# Generated outputs: SMC-enriched subset `sce_exp.csv`
# Upstream dependencies: Final SMC Seurat object
# Downstream consumers: `create_smc_scenic_loom.py`
# Configuration keys: DATA_PROCESSED_DIR, SCENIC_DB_DIR, SCENIC_TF_LIST, SCENIC_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/04_CREB3L2_regulatory_network/03-3_sub-SCENIC.R
# Reproducibility notes: This is the manuscript SMC-enriched subset SCENIC
# input export. It is independent of the all-cell SCENIC run.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
output_file <- configured_path(
  PATHS,
  "SCENIC_RESULTS_DIR",
  "01.sub_allCell",
  "sce_exp.csv"
)
output_parent <- dirname(output_file)
dir.create(
  output_parent,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(output_parent)) {
  stop(
    "Failed to create output directory: ",
    output_parent
  )
}
sce <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
write.csv(t(as.matrix(sce@assays$RNA@counts)),file = output_file)
# Next: run create_smc_scenic_loom.py to create sce.loom.
