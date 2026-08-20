# =============================================================================
# Title: Build the seven-sample downstream dataset
# Purpose: Construct the seven-sample calcified-artery object used by all
# downstream computational modules.
# Manuscript Figure/Table: Upstream input for Figures 2-4 and Figure 6
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `8data_cca_3000.rds`; seven-sample metadata
# Generated outputs: `7data_preannotation.rds`; seven-sample cluster markers
# and diagnostic plots
# Upstream dependencies: `01_umap.R`; `8data_cca_3000.rds`
# Downstream consumers: `05_annotate_7sample_dataset.R`
# Configuration keys: PROJECT_ROOT, DATA_RAW_DIR, DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, ATLAS_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/00_scRNA_atlas/04_build_7sample_dataset.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from the original workflow; only path/configuration wiring, output naming and
# comments were clarified.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(cowplot)
#library(SCP)
dir.create(PATHS$ATLAS_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
obj = readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "8data_cca_3000.rds"))
Idents(obj) <- 'sample_name'
sub <- subset(obj, idents = c('Mild_1','Mild_2','Mild_3','Mild_4','Severe_1','Severe_2','Severe_3'))
rm(obj)
obj <- CreateSeuratObject(counts = sub@assays$RNA@counts,meta.data= sub@meta.data)

obj.list <- SplitObject(object = obj,split.by = 'sample_name')
obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
features <- SelectIntegrationFeatures(object.list = obj.list)
immune.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
obj <- immune.combined
rm(immune.combined)

DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj,verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:27,n.neighbors =50,min.dist=0.4)
obj <- FindNeighbors(obj, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.5)
DefaultAssay(obj) <- "RNA"
DimPlot(obj, label = T)

# `Celltype` inherited from the eight-sample object is provisional here.
# The seven-sample annotation is rebuilt in 05_annotate_7sample_dataset.R.
saveRDS(
  obj,
  file = configured_path(
    PATHS,
    "DATA_PROCESSED_DIR",
    "7data_preannotation.rds"
  )
)
##################
dir.create(path = configured_path(PATHS, "ATLAS_RESULTS_DIR", "00_marker", trailing = TRUE), 
           recursive = TRUE, 
           showWarnings = TRUE)
ma <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue1", "pheatmap", "order.csv"),header = T)
ge <- ma$X

for (gene in ge) {
  try({
    if (gene %in% rownames(obj)) {
      p <- FeaturePlot(obj, features = gene)
      ggsave(
        filename = paste0(
          configured_path(PATHS, "ATLAS_RESULTS_DIR", "00_marker", trailing = TRUE),
          gene, '.png'
        ),
        plot = p
      )
    } else {
      message("跳过不存在基因:", gene)
    }
  }, silent = TRUE)
}
#################
marker = read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "atlas", "marker.csv"))
df = data.frame(rownames(obj@assays$RNA@counts))
colnames(df) = 'Marker'
marker = merge(marker,df,by = 'Marker',sort =F)

dir.create(
  configured_path(PATHS, "ATLAS_RESULTS_DIR", "Marker"),
  recursive = TRUE,
  showWarnings = FALSE
)
for (i in 1:nrow(marker)){
p = FeaturePlot(obj, features = marker[i,1])
ggsave(
  filename = paste0(
    configured_path(PATHS, "ATLAS_RESULTS_DIR", "Marker", trailing = TRUE),
    marker[i,2], '_', marker[i,1], '.png'
  ),
  plot = p,
  width = 6,
  height = 5
)
}

m <- FindAllMarkers(obj,only.pos = TRUE,  logfc.threshold = 1)
write.csv(
  m,
  configured_path(PATHS, "ATLAS_RESULTS_DIR", "7data_cluster_markers.csv"),
  quote = FALSE
)
###################
