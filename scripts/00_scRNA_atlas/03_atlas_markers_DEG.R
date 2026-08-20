# =============================================================================
# Title: Atlas cell-type marker and DEG source data
# Purpose: Produce the manuscript major-cell marker/DEG table and canonical
# marker source panels from the final eight-sample atlas.
# Manuscript Figure/Table: Figure S1D; Table S1
# Release role: SUPPLEMENTARY SOURCE DATA (not a composite final figure)
# Required inputs: DATA_PROCESSED_DIR/8data_cca_3000.rds
# Generated outputs: DATA_PROCESSED_DIR/8data_cca_3000_RunDEtest.rds;
# ATLAS_RESULTS_DIR/8data_DEGs.csv; canonical marker panels
# Upstream dependencies: 01_umap.R
# Downstream consumers: Table S1 and Figure S1D assembly
# Configuration keys: DATA_PROCESSED_DIR, RESULTS_DIR, ATLAS_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/00_scRNA_atlas/03_atlas_markers_DEG.R
# Reproducibility notes: SCP RunDEtest, fold-change and adjusted-P-value
# thresholds and the canonical marker list are unchanged. The DEG-enriched
# object is now produced explicitly by this script.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(scop)
library(SCP)

dir.create(PATHS$ATLAS_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
marker_output_dir <- configured_path(
  PATHS, "RESULTS_DIR", "figures", "Supplementary", "sup_fig1-gene"
)
dir.create(marker_output_dir, recursive = TRUE, showWarnings = FALSE)

final_celltype_levels <- c(
  "Neutrophil", "Fibroblast", "non-OFB-SMC", "Endothelial cell",
  "CD8 T cell", "OFB-SMC", "CD4 T cell", "NK Cell",
  "Dendritic cell", "Mixed myeloid cell", "Macrophage",
  "ProMacs", "B cell", "Eosinophil", "Mast cell"
)

obj <- readRDS(configured_path(
  PATHS, "DATA_PROCESSED_DIR", "8data_cca_3000.rds"
))
obj$Celltype <- factor(obj$Celltype, levels = final_celltype_levels)
DefaultAssay(obj) <- "RNA"

obj <- RunDEtest(
  srt = obj,
  group_by = "Celltype",
  fc.threshold = 1,
  only.pos = FALSE
)
saveRDS(
  obj,
  configured_path(
    PATHS, "DATA_PROCESSED_DIR", "8data_cca_3000_RunDEtest.rds"
  )
)

DEGs <- obj@tools$DEtest_Celltype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
write.csv(
  DEGs,
  configured_path(PATHS, "ATLAS_RESULTS_DIR", "8data_DEGs.csv"),
  row.names = FALSE
)

genes <- c(
  "VWF", "S100A8", "DCN", "MYH11", "CD8A",
  "POSTN", "CD4", "NCR1", "GPR183", "GCH1",
  "CCL3", "LYVE1", "CD79A", "SIGLEC8", "TPSAB1"
)
for (gene in genes) {
  p_gene <- FeatureDimPlot(obj, features = gene, theme_use = "theme_blank")
  ggsave(
    filename = file.path(marker_output_dir, paste0(gene, ".png")),
    plot = p_gene,
    width = 6,
    height = 6
  )
  ggsave(
    filename = file.path(marker_output_dir, paste0(gene, ".pdf")),
    plot = p_gene,
    width = 6,
    height = 6
  )
}
