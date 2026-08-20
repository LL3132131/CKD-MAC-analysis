# =============================================================================
# Title: 03_SCENIC_mk_count.R
# Purpose: all-cell loom conversion
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data_umap.rds`
# Generated outputs: all-cell `sce_exp.csv`
# Upstream dependencies: `7data_umap.rds`
# Downstream consumers: all-cell loom conversion
# Configuration keys: DATA_PROCESSED_DIR, SCENIC_DB_DIR, SCENIC_TF_LIST, SCENIC_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/04_CREB3L2_regulatory_network/03_SCENIC_mk_count.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(viridis)
output_dir <- configured_path(PATHS, "SCENIC_RESULTS_DIR", "00.allCell")
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(output_dir)) {
  stop(
    "Failed to create output directory: ",
    output_dir
  )
}
scenic_export_dir <- configured_path(
  PATHS,
  "SCENIC_RESULTS_DIR",
  "00.allCell",
  "7scenic"
)
dir.create(
  scenic_export_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(scenic_export_dir)) {
  stop(
    "Failed to create output directory: ",
    scenic_export_dir
  )
}
setwd(output_dir)
obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
count <- obj@assays$RNA@counts

###for-loom-SCE
write.csv(t(as.matrix(obj@assays$RNA@counts)),file = configured_path(PATHS, "SCENIC_RESULTS_DIR", "00.allCell", "7scenic", "sce_exp.csv"))
# Next: run create_all_cell_scenic_loom.py to create rSCENIC.loom.
