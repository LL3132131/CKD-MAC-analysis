# =============================================================================
# Title: Eight-sample atlas annotation and UMAP
# Purpose: Reintegrate the final eight samples, preserve the manuscript cluster
# assignments, apply the final 15 major-cell labels and generate atlas panels.
# Manuscript Figure/Table: Figure 1B, Figure 1D, Figure S1B-C,
# Supplementary Figure 2 (bottom)
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: DATA_PROCESSED_DIR/Annotation_Integrated.Find.rds; atlas
# marker and ordering tables under DATA_EXTERNAL_DIR
# Generated outputs: DATA_PROCESSED_DIR/8data_cca_3000.rds; annotation markers;
# atlas UMAP/QC/composition source panels
# Upstream dependencies: 00_scRNA_QC_integration.R
# Downstream consumers: 02_major_cell_heatmap.R, 03_atlas_markers_DEG.R,
# 04_build_7sample_dataset.R
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, ATLAS_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/00_scRNA_atlas/01_umap.R
# Reproducibility notes: Integration, PCA/UMAP, clustering and reclustering
# parameters and cluster membership are unchanged. Only public metadata and final
# manuscript display labels are standardized.
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
library(colorspace)

dir.create(PATHS$ATLAS_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
marker_dir <- configured_path(PATHS, "ATLAS_RESULTS_DIR", "Marker")
order_marker_dir <- configured_path(PATHS, "ATLAS_RESULTS_DIR", "00_marker")
dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(order_marker_dir, recursive = TRUE, showWarnings = FALSE)

public_sample_labels <- c(
  "Control",
  "Mild_1", "Mild_2", "Mild_3", "Mild_4",
  "Severe_1", "Severe_2", "Severe_3"
)
final_celltype_levels <- c(
  "Neutrophil", "Fibroblast", "non-OFB-SMC", "Endothelial cell",
  "CD8 T cell", "OFB-SMC", "CD4 T cell", "NK Cell",
  "Dendritic cell", "Mixed myeloid cell", "Macrophage",
  "ProMacs", "B cell", "Eosinophil", "Mast cell"
)
celltype_colors <- c(
  "Neutrophil" = "#e366b6",
  "Fibroblast" = "#cd177b",
  "non-OFB-SMC" = "#3C234A",
  "Endothelial cell" = "#004d43",
  "CD8 T cell" = "#486eca",
  "OFB-SMC" = "#8b0254",
  "CD4 T cell" = "#5fb5f3",
  "NK Cell" = "#F7B46D",
  "Dendritic cell" = "#1f71b5",
  "Mixed myeloid cell" = "#0b0baa",
  "Macrophage" = "#435B95",
  "ProMacs" = "#8fb0ff",
  "B cell" = "#ff3a0e",
  "Eosinophil" = "#620d14",
  "Mast cell" = "#de545a"
)
sample_colors <- qualitative_hcl(n = length(public_sample_labels), palette = "Set3")
group_colors <- qualitative_hcl(n = 3, palette = "Set3")

calcification_group_from_sample <- function(sample_name) {
  group <- rep(NA_character_, length(sample_name))
  group[sample_name == "Control"] <- "Non-calcified"
  group[grepl("^Mild_", sample_name)] <- "Mild"
  group[grepl("^Severe_", sample_name)] <- "Severe"
  factor(group, levels = c("Non-calcified", "Mild", "Severe"))
}

obj <- readRDS(configured_path(
  PATHS, "DATA_PROCESSED_DIR", "Annotation_Integrated.Find.rds"
))
if (!"sample_name" %in% colnames(obj@meta.data)) {
  stop("Required metadata field is missing: sample_name")
}
unexpected_samples <- setdiff(unique(as.character(obj$sample_name)), public_sample_labels)
missing_samples <- setdiff(public_sample_labels, unique(as.character(obj$sample_name)))
if (length(unexpected_samples) > 0 || length(missing_samples) > 0) {
  stop("The atlas object must contain exactly the eight public sample labels.")
}
obj$sample_name <- factor(as.character(obj$sample_name), levels = public_sample_labels)
obj$calcification_group <- calcification_group_from_sample(as.character(obj$sample_name))

Idents(obj) <- "sample_name"
ob <- subset(obj, idents = public_sample_labels)
rm(obj)
obj <- CreateSeuratObject(counts = ob@assays$RNA@counts, meta.data = ob@meta.data)
obj$orig.ident <- as.character(obj$sample_name)
ambiguous_metadata <- intersect(c("sample", "batch", "orig"), colnames(obj@meta.data))
if (length(ambiguous_metadata) > 0) {
  obj@meta.data <- obj@meta.data[, setdiff(
    colnames(obj@meta.data), ambiguous_metadata
  ), drop = FALSE]
}

obj.list <- SplitObject(object = obj, split.by = "sample_name")
obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = obj.list)
immune.anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features
)
immune.combined <- IntegrateData(anchorset = immune.anchors)
obj <- immune.combined
rm(immune.combined)

set.seed(42)
DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = 1:34,
  n.neighbors = 50,
  min.dist = 0.5
)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:49)
obj <- FindClusters(obj, resolution = 0.8)

order_table <- read.csv(configured_path(
  PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue1",
  "pheatmap", "order.csv"
), header = TRUE)
for (gene in order_table$X) {
  if (gene %in% rownames(obj)) {
    p_gene <- FeaturePlot(obj, features = gene)
    ggsave(
      filename = file.path(order_marker_dir, paste0(gene, ".png")),
      plot = p_gene
    )
  }
}

marker_table <- read.csv(configured_path(
  PATHS, "DATA_EXTERNAL_DIR", "atlas", "marker.csv"
))
available_genes <- data.frame(Marker = rownames(obj@assays$RNA@counts))
marker_table <- merge(marker_table, available_genes, by = "Marker", sort = FALSE)
for (i in seq_len(nrow(marker_table))) {
  p_gene <- FeaturePlot(obj, features = marker_table[i, 1])
  ggsave(
    filename = file.path(
      marker_dir,
      paste0(marker_table[i, 2], "_", marker_table[i, 1], ".png")
    ),
    plot = p_gene,
    width = 6,
    height = 5
  )
}

cluster_markers <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 1)
write.csv(
  cluster_markers,
  configured_path(
    PATHS, "ATLAS_RESULTS_DIR", "8sample_cluster_markers_for_annotation.csv"
  ),
  quote = FALSE
)

# Final manuscript labels. Cluster membership is preserved exactly from the
# original analysis; only display names are standardized.
obj$Celltype <- ""
obj$Celltype[obj$seurat_clusters %in% c("1", "24")] <- "Fibroblast"
obj$Celltype[obj$seurat_clusters %in% "9"] <- "OFB-SMC"
obj$Celltype[obj$seurat_clusters %in% "20"] <- "B cell"
obj$Celltype[obj$seurat_clusters %in% "7"] <- "CD4 T cell"
obj$Celltype[obj$seurat_clusters %in% "12"] <- "CD8 T cell"
obj$Celltype[obj$seurat_clusters %in% "11"] <- "NK Cell"
obj$Celltype[obj$seurat_clusters %in% c("4", "13", "17", "18", "21")] <-
  "Endothelial cell"
obj$Celltype[obj$seurat_clusters %in% c("0", "3", "28", "16")] <-
  "non-OFB-SMC"
obj$Celltype[obj$seurat_clusters %in% "19"] <- "Mixed myeloid cell"
obj$Celltype[obj$seurat_clusters %in% c("15", "6")] <- "Macrophage"
obj$Celltype[obj$seurat_clusters %in% c("8", "2", "26", "5", "10", "14", "27")] <-
  "Neutrophil"
obj$Celltype[obj$seurat_clusters %in% "23"] <- "Eosinophil"
obj$Celltype[obj$seurat_clusters %in% "25"] <- "Dendritic cell"
obj$Celltype[obj$seurat_clusters %in% "22"] <- "Mast cell"
if (any(obj$Celltype == "")) {
  stop("At least one cluster has no final major-cell annotation.")
}

# Preserve the original cluster-15 reclustering parameters and assignments.
DefaultAssay(obj) <- "integrated"
Idents(obj) <- "seurat_clusters"
sub <- subset(obj, idents = "15")
sub <- FindNeighbors(sub, dims = 1:30)
sub <- FindClusters(sub, resolution = 0.1)
obj$recluster15 <- "missing"
sub_index <- match(colnames(sub), colnames(obj))
if (anyNA(sub_index)) {
  stop("Cluster-15 cell barcodes do not match the parent atlas exactly.")
}
obj$recluster15[sub_index] <- as.character(sub$seurat_clusters)
obj$Celltype[obj$recluster15 %in% c("0", "1")] <- "Macrophage"
obj$Celltype[obj$recluster15 == "2"] <- "ProMacs"
obj$Celltype <- factor(obj$Celltype, levels = final_celltype_levels)

DefaultAssay(obj) <- "RNA"
p_group <- CellDimPlot(
  srt = obj,
  reduction = "UMAP",
  palcolor = group_colors,
  group.by = "calcification_group",
  theme_args = list(base_size = 16)
)
p_sample <- CellDimPlot(
  srt = obj,
  reduction = "UMAP",
  cols = sample_colors,
  group.by = "sample_name",
  theme_args = list(base_size = 16)
)
p_cluster <- CellDimPlot(
  srt = obj,
  reduction = "UMAP",
  group.by = "seurat_clusters",
  label = TRUE,
  theme_args = list(base_size = 16)
)
p_celltype <- CellDimPlot(
  srt = obj,
  reduction = "UMAP",
  group.by = "Celltype",
  cols = unname(celltype_colors),
  theme_args = list(base_size = 16)
)
figure1_plot_obj <- obj
figure1_plot_obj$figure1_condition <- factor(
  ifelse(
    figure1_plot_obj$calcification_group == "Non-calcified",
    "Non-calcified",
    "Calcified"
  ),
  levels = c("Non-calcified", "Calcified")
)
p_split <- DimPlot(figure1_plot_obj, split.by = "figure1_condition")
rm(figure1_plot_obj)

ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "calcification_group.png"),
  plot = p_group, width = 10, height = 6
)
ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "sample_name.png"),
  plot = p_sample, width = 10, height = 6
)
ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "seurat_clusters.png"),
  plot = p_cluster, width = 10, height = 6
)
ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "Figure1B_celltypes.pdf"),
  plot = p_celltype, width = 5, height = 5, dpi = 900
)
ggsave(
  filename = configured_path(
    PATHS, "ATLAS_RESULTS_DIR", "Figure1D_calcification_groups.png"
  ),
  plot = p_split, width = 15, height = 6
)

p_nfeature <- FeatureStatPlot(
  obj, stat.by = "nFeature_RNA", group.by = "sample_name"
)
p_ncount <- FeatureStatPlot(
  obj, stat.by = "nCount_RNA", group.by = "sample_name"
)
p_mt <- FeatureStatPlot(
  obj, stat.by = "percent.mt", group.by = "sample_name"
)
ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "nFeature_RNA.png"),
  plot = p_nfeature, width = 10, height = 6
)
ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "nCount_RNA.png"),
  plot = p_ncount, width = 10, height = 6
)
ggsave(
  filename = configured_path(PATHS, "ATLAS_RESULTS_DIR", "percent.mt.png"),
  plot = p_mt, width = 10, height = 6
)

# Supplementary Figure 2 (bottom): final major-cell proportions per sample.
# The legend specifies one bar per individual sample. No group-level mean,
# pooled group proportion or error statistic is calculated.
observed_composition <- obj@meta.data %>%
  count(sample_name, Celltype, name = "cell_count")
composition <- tidyr::expand_grid(
  sample_name = public_sample_labels,
  Celltype = final_celltype_levels
) %>%
  left_join(
    observed_composition,
    by = c("sample_name", "Celltype")
  ) %>%
  mutate(
    cell_count = tidyr::replace_na(cell_count, 0L),
    sample_name = factor(sample_name, levels = public_sample_labels),
    calcification_group = calcification_group_from_sample(
      as.character(sample_name)
    ),
    Celltype = factor(Celltype, levels = final_celltype_levels)
  ) %>%
  group_by(sample_name) %>%
  mutate(
    sample_total_cells = sum(cell_count),
    proportion = cell_count / sample_total_cells
  ) %>%
  ungroup()
write.csv(
  composition,
  configured_path(
    PATHS, "ATLAS_RESULTS_DIR", "celltype_proportion_by_sample.csv"
  ),
  row.names = FALSE
)
p_composition <- ggplot(
  composition,
  aes(x = sample_name, y = proportion, fill = Celltype)
) +
  geom_col() +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = configured_path(
    PATHS,
    "ATLAS_RESULTS_DIR",
    "Supplementary_Figure_2_celltype_proportion_by_sample.pdf"
  ),
  plot = p_composition,
  width = 10,
  height = 6
)
ggsave(
  filename = configured_path(
    PATHS,
    "ATLAS_RESULTS_DIR",
    "Supplementary_Figure_2_celltype_proportion_by_sample.png"
  ),
  plot = p_composition,
  width = 10,
  height = 6,
  dpi = 300
)

saveRDS(
  obj,
  configured_path(PATHS, "DATA_PROCESSED_DIR", "8data_cca_3000.rds")
)
