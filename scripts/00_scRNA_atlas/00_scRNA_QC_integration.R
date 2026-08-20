# =============================================================================
# Title: 00_scRNA_QC_integration.R
# Purpose: Perform sample-level QC, doublet removal, CCA integration and initial
# clustering for the eight-sample Figure 1 atlas.
# Manuscript Figure/Table: Figure 1; Figure S1A-C
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: Eight de-identified single-sample `*_Seurat.rds` objects;
# command-line QC/integration arguments
# Generated outputs: Integrated Seurat object; `Annotation_Integrated.Find.rds`;
# QC, sample and cluster UMAP plots
# Upstream dependencies: Single-sample Seurat objects
# Downstream consumers: `01_umap.R`
# Configuration keys: DATA_PROCESSED_DIR, ATLAS_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/00_scRNA_atlas/00_scRNA_QC_integration.R [arguments]
# Reproducibility notes: QC thresholds, DoubletFinder settings, CCA integration
# dimensions, PCA/UMAP dimensions and clustering resolutions are unchanged.
# Input filenames must already use the eight public sample labels.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

parser <- argparse::ArgumentParser(description = "Script to QC and cluster scRNA data")
parser$add_argument("-I", "--input", help = "input directory")
parser$add_argument("-D", "--id", help = "sample ID pattern")
parser$add_argument("-G", "--filterGene", help = "gene number threshold")
parser$add_argument("-C", "--CutoffGene", help = "upper gene number threshold")
parser$add_argument("-M", "--filterMito", help = "mitochondrial percentage threshold")
parser$add_argument("-O", "--out", help = "output directory")
parser$add_argument("-F", "--double", help = "expected doublet ratio")
parser$add_argument("-P", "--pc", help = "number of principal components")
parser$add_argument("-V", "--vg", help = "number of variable genes")
parser$add_argument("-H", "--hb", help = "hemoglobin percentage threshold")
args <- parser$parse_args()

library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(cowplot)

public_sample_labels <- c(
  "Control",
  "Mild_1", "Mild_2", "Mild_3", "Mild_4",
  "Severe_1", "Severe_2", "Severe_3"
)

calcification_group_from_sample <- function(sample_name) {
  group <- rep(NA_character_, length(sample_name))
  group[sample_name == "Control"] <- "Non-calcified"
  group[grepl("^Mild_", sample_name)] <- "Mild"
  group[grepl("^Severe_", sample_name)] <- "Severe"
  factor(group, levels = c("Non-calcified", "Mild", "Severe"))
}

dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$ATLAS_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(args$input)

file_list <- list.files(pattern = args$id)
file_list <- gsub("_Seurat.rds", "", file_list)
unexpected_samples <- setdiff(file_list, public_sample_labels)
missing_samples <- setdiff(public_sample_labels, file_list)
if (length(unexpected_samples) > 0 || length(missing_samples) > 0) {
  stop(
    "The atlas input must contain exactly the eight de-identified public ",
    "sample labels."
  )
}
file_list <- public_sample_labels

scRNAlist <- list()
for (i in file_list) {
  obj <- readRDS(paste0(i, "_Seurat.rds"))
  obj$sample_name <- i
  obj$calcification_group <- calcification_group_from_sample(i)
  scRNAlist[[i]] <- obj
}

obj <- merge(
  scRNAlist[[1]],
  y = scRNAlist[2:length(file_list)],
  add.cell.ids = names(scRNAlist)
)
obj$orig.ident <- as.character(obj$sample_name)
ambiguous_metadata <- intersect(c("sample", "batch", "orig"), colnames(obj@meta.data))
if (length(ambiguous_metadata) > 0) {
  obj@meta.data <- obj@meta.data[, setdiff(
    colnames(obj@meta.data), ambiguous_metadata
  ), drop = FALSE]
}
HB.genes <- c(
  "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2",
  "HBM", "HBQ1", "HBZ"
)
HB_m <- match(HB.genes, rownames(obj@assays$RNA))
HB.genes <- rownames(obj@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj[["HB.genes"]] <- PercentageFeatureSet(obj, features = HB.genes)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(
  obj,
  subset =
    nFeature_RNA > as.numeric(args$filterGene) &
    nFeature_RNA < as.numeric(args$CutoffGene) &
    percent.mt < as.numeric(args$filterMito) &
    HB.genes < as.numeric(args$hb)
)
obj.list <- SplitObject(object = obj, split.by = "sample_name")

Find_doublet <- function(data) {
  sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(as.numeric(args$double) * ncol(data))
  p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
  data <- doubletFinder_v3(
    data,
    PCs = 1:10,
    pN = 0.25,
    pK = p,
    nExp = nExp_poi,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
  data
}

doublet <- list()
setwd(args$out)
for (i in 1:length(x = obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(
    object = obj.list[[i]],
    selection.method = "vst",
    nfeatures = as.numeric(args$vg),
    verbose = FALSE
  )
  obj.list[[i]] <- ScaleData(obj.list[[i]])
  obj.list[[i]] <- RunPCA(obj.list[[i]])
  obj.list[[i]] <- RunUMAP(obj.list[[i]], dims = 1:10)
  obj.list[[i]] <- Find_doublet(obj.list[[i]])
  doublet[[i]] <- obj.list[[i]]@meta.data[, -6]
  obj.list[[i]] <- subset(obj.list[[i]], subset = doublet_info == "Singlet")
}

reference.list <- obj.list
obj.anchors <- FindIntegrationAnchors(
  object.list = reference.list,
  dims = 1:30,
  reduction = "cca"
)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

DefaultAssay(object = obj.integrated) <- "integrated"
obj.integrated <- FindVariableFeatures(obj.integrated)
obj.integrated <- ScaleData(object = obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(object = obj.integrated, verbose = FALSE)
obj.integrated <- RunUMAP(
  object = obj.integrated,
  reduction = "pca",
  dims = 1:as.numeric(args$pc)
)
obj.integrated <- FindNeighbors(
  object = obj.integrated,
  dims = 1:as.numeric(args$pc)
)
obj.integrated.Find <- FindClusters(
  object = obj.integrated,
  resolution = 0.5
)
obj.integrated.Find <- RunTSNE(
  object = obj.integrated.Find,
  dims = 1:as.numeric(args$pc)
)
DefaultAssay(obj.integrated.Find) <- "RNA"

integration_path <- file.path(
  args$out,
  paste0("Seurat_", args$id, ".integrated.Find.rds")
)
saveRDS(obj.integrated.Find, file = integration_path)

colorlist <- c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
  "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
  "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
  "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
  "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
  "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
  "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
  "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
  "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
  "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
  "#191970", "#000080", "#6495ED", "#1E90FF", "#00BFFF",
  "#00FFFF", "#FF1493", "#FF00FF", "#A020F0", "#63B8FF",
  "#008B8B", "#54FF9F", "#00FF00", "#76EE00", "#FFF68F"
)

p_sample <- DimPlot(
  obj.integrated.Find,
  cols = colorlist,
  group.by = "sample_name"
)
ggsave(
  filename = file.path(args$out, paste0(args$id, ".sample_name.png")),
  plot = p_sample,
  width = 8,
  height = 8
)

p_sample_cluster <- DimPlot(
  obj.integrated.Find,
  cols = colorlist,
  group.by = "seurat_clusters",
  split.by = "sample_name"
)
ggsave(
  filename = file.path(args$out, paste0(args$id, ".sample_clusters.png")),
  plot = p_sample_cluster,
  width = 20,
  height = 8
)

p_cluster <- DimPlot(
  obj.integrated.Find,
  cols = colorlist,
  group.by = "seurat_clusters"
)
ggsave(
  filename = file.path(args$out, paste0(args$id, ".UMAP.png")),
  plot = p_cluster,
  width = 8,
  height = 8
)

obj <- readRDS(integration_path)
setwd(PATHS$ATLAS_RESULTS_DIR)
obj$orig.ident <- as.character(obj$sample_name)
ambiguous_metadata <- intersect(c("sample", "batch", "orig"), colnames(obj@meta.data))
if (length(ambiguous_metadata) > 0) {
  obj@meta.data <- obj@meta.data[, setdiff(
    colnames(obj@meta.data), ambiguous_metadata
  ), drop = FALSE]
}
DefaultAssay(obj) <- "integrated"
obj <- FindClusters(object = obj, resolution = 0.8)

p_resolution <- DimPlot(obj, cols = colorlist, group.by = "seurat_clusters")
ggsave(
  filename = configured_path(
    PATHS,
    "ATLAS_RESULTS_DIR",
    "resolution_0.8_UMAP.png"
  ),
  plot = p_resolution,
  width = 8,
  height = 8
)

cluster_dir <- configured_path(
  PATHS,
  "ATLAS_RESULTS_DIR",
  "resolution_0.8_clusters"
)
dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
meta <- obj@meta.data
for (i in unique(meta$seurat_clusters)) {
  cell <- rownames(meta[meta$seurat_clusters == i, ])
  p <- DimPlot(
    obj,
    group.by = "seurat_clusters",
    cells.highlight = cell,
    pt.size = 0.00001
  ) +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(
    filename = file.path(cluster_dir, paste0(i, ".png")),
    plot = p,
    width = 5,
    height = 4
  )
}

p_sample_cluster <- DimPlot(
  obj,
  cols = colorlist,
  group.by = "seurat_clusters",
  split.by = "sample_name"
)
ggsave(
  filename = configured_path(
    PATHS,
    "ATLAS_RESULTS_DIR",
    "resolution_0.8_clusters_by_sample.png"
  ),
  plot = p_sample_cluster,
  width = 50,
  height = 8,
  limitsize = FALSE
)

saveRDS(
  obj,
  configured_path(
    PATHS,
    "DATA_PROCESSED_DIR",
    "Annotation_Integrated.Find.rds"
  )
)
