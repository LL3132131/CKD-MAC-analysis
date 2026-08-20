# =============================================================================
# Title: 01-4-01_prepare_umap-gene.R
# Purpose: Dynamo
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: Final SMC Seurat object
# Generated outputs: Dynamo metadata, UMAP and top-gene tables
# Upstream dependencies: Final SMC Seurat object
# Downstream consumers: Dynamo
# Configuration keys: DATA_PROCESSED_DIR, TRAJECTORY_RESULTS_DIR, MONOCLE2_PRIVATE_RDS, MONOCLE2_PUBLIC_RDS
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/03_Trajectory/01-4-01_prepare_umap-gene.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(viridis)
library(SeuratWrappers)
library(monocle3)
output_dir <- configured_path(PATHS, "TRAJECTORY_RESULTS_DIR", "dynamo")
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
setwd(output_dir)

obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
if (!"Celltype" %in% colnames(obj@meta.data)) {
  stop("smc_subclusters_final_seurat.rds is missing required metadata: Celltype")
}
Idents(obj) <- "Celltype"
umap <- obj@reductions$umap@cell.embeddings
meta <- obj@meta.data

write.table(meta , file = 'sub-meta-dynamo.xls',sep = '\t')
write.table(umap,file = 'sub-umap-dynamo.xls',sep = '\t')

gene <- FindMarkers(obj, ident.1 = 'OFB-SMC',log2FC.threshold = 0.5,only.pos = T,test.use = "wilcox")
gene1 <- gene[gene$avg_log2FC > 1 & gene$p_val < 0.05, ]
gene_name <- rownames(gene1)
write.table(gene_name,file = 'topgene-for-recipe.xls',sep = '\t')
