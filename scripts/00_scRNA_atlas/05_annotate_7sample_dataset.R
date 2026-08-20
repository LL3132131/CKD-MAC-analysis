# =============================================================================
# Title: Rebuild the final seven-sample major-cell annotation
# Purpose: Apply the original seven-sample whole-atlas annotation and branch
# reclustering logic to the pre-annotation object.
# Manuscript Figure/Table: Figure 2A
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data_preannotation.rds`
# Generated outputs: `7data_umap.rds`; `SCP_umap.png`
# Upstream dependencies: `04_build_7sample_dataset.R`
# Downstream consumers: SMC/OFB-SMC, CellChat, scDRS, all-cell SCENIC and
# validation modules
# Configuration keys: PROJECT_ROOT, DATA_PROCESSED_DIR, ATLAS_RESULTS_DIR
# Expected environment: R with Seurat, SCP and ggplot2
# Example run command:
# Rscript scripts/00_scRNA_atlas/05_annotate_7sample_dataset.R
# Reproducibility notes: Cluster assignments, branch resolutions and plotting
# dimensions are restored from the original 00_zhushi.R workflow. Branch labels
# are written back to the full object by exact cell barcode. No marker-based
# reannotation is introduced.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(SCP)
library(ggplot2)

input_rds <- configured_path(
  PATHS,
  "DATA_PROCESSED_DIR",
  "7data_preannotation.rds"
)
output_rds <- configured_path(
  PATHS,
  "DATA_PROCESSED_DIR",
  "7data_umap.rds"
)
figure_2a_file <- configured_path(
  PATHS,
  "ATLAS_RESULTS_DIR",
  "SCP_umap.png"
)

dir.create(
  PATHS$ATLAS_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(PATHS$ATLAS_RESULTS_DIR)) {
  stop("Failed to create output directory: ", PATHS$ATLAS_RESULTS_DIR)
}

obj <- readRDS(input_rds)

required_metadata <- c(
  "sample_name",
  "calcification_group",
  "seurat_clusters"
)
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
if (length(missing_metadata) > 0) {
  stop(
    "7data_preannotation.rds is missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

expected_samples <- c(
  "Mild_1", "Mild_2", "Mild_3", "Mild_4",
  "Severe_1", "Severe_2", "Severe_3"
)
observed_samples <- unique(as.character(obj$sample_name))
if (!setequal(observed_samples, expected_samples)) {
  stop(
    "7data_preannotation.rds must contain exactly the seven documented ",
    "samples. Observed: ",
    paste(sort(observed_samples), collapse = ", ")
  )
}

cell_barcodes <- colnames(obj)
if (anyDuplicated(cell_barcodes)) {
  stop("7data_preannotation.rds contains duplicated cell barcodes.")
}

expected_initial_clusters <- as.character(0:23)
observed_initial_clusters <- unique(as.character(obj$seurat_clusters))
if (!setequal(observed_initial_clusters, expected_initial_clusters)) {
  stop(
    "The initial seven-sample clusters do not match the original 0-23 ",
    "annotation interface. Observed: ",
    paste(sort(observed_initial_clusters), collapse = ", ")
  )
}

write_branch_celltypes <- function(
  full_object,
  branch_object,
  cluster_to_celltype,
  branch_name
) {
  branch_cells <- colnames(branch_object)
  if (anyDuplicated(branch_cells)) {
    stop(branch_name, " contains duplicated cell barcodes.")
  }
  missing_cells <- setdiff(branch_cells, colnames(full_object))
  if (length(missing_cells) > 0) {
    stop(
      branch_name,
      " contains cell barcodes absent from the full object: ",
      paste(missing_cells, collapse = ", ")
    )
  }

  branch_clusters <- as.character(branch_object$seurat_clusters)
  unmapped_clusters <- setdiff(
    unique(branch_clusters),
    names(cluster_to_celltype)
  )
  if (length(unmapped_clusters) > 0) {
    stop(
      branch_name,
      " contains unmapped recluster identities: ",
      paste(sort(unmapped_clusters), collapse = ", ")
    )
  }

  full_object@meta.data[branch_cells, "Celltype"] <-
    unname(cluster_to_celltype[branch_clusters])
  full_object
}

# Original whole-object cluster annotation. Clusters 13, 15, 16, 18 and 20 are
# intentionally left blank until their original branch reclustering steps.
obj$Celltype <- ""
initial_cluster_to_celltype <- c(
  "0" = "NE",
  "1" = "SMC",
  "2" = "SMC",
  "3" = "NE",
  "4" = "NE",
  "5" = "CD4_help_T_Cell",
  "6" = "Macrophage",
  "7" = "SMC",
  "8" = "CD8_T_Cell",
  "9" = "EC",
  "10" = "FB",
  "11" = "FB_SMC",
  "12" = "NK",
  "14" = "FB",
  "17" = "EC",
  "19" = "EC",
  "21" = "Macrophage_NE",
  "22" = "B_Cell",
  "23" = "EC"
)
initial_clusters <- as.character(obj$seurat_clusters)
initial_labels <- unname(initial_cluster_to_celltype[initial_clusters])
known_initial_cells <- !is.na(initial_labels)
obj$Celltype[known_initial_cells] <- initial_labels[known_initial_cells]

DefaultAssay(obj) <- "integrated"

# Original clusters 16/20: resolution 0.2.
Idents(obj) <- "seurat_clusters"
branch_16_20 <- subset(obj, idents = c("16", "20"))
branch_16_20 <- FindNeighbors(branch_16_20, dims = 1:30)
branch_16_20 <- FindClusters(branch_16_20, resolution = 0.2)
obj <- write_branch_celltypes(
  obj,
  branch_16_20,
  c(
    "0" = "SMC",
    "1" = "unKnown",
    "2" = "unKnown",
    "3" = "FB"
  ),
  "Original clusters 16/20"
)

# Original cluster 13: resolution 0.2.
Idents(obj) <- "seurat_clusters"
branch_13 <- subset(obj, idents = "13")
branch_13 <- FindNeighbors(branch_13, dims = 1:30)
branch_13 <- FindClusters(branch_13, resolution = 0.2)
obj <- write_branch_celltypes(
  obj,
  branch_13,
  c(
    "0" = "SMC",
    "1" = "FB_SMC"
  ),
  "Original cluster 13"
)

# Original cluster 15: resolution 0.5.
Idents(obj) <- "seurat_clusters"
branch_15 <- subset(obj, idents = "15")
branch_15 <- FindNeighbors(branch_15, dims = 1:30)
branch_15 <- FindClusters(branch_15, resolution = 0.5)
obj <- write_branch_celltypes(
  obj,
  branch_15,
  c(
    "0" = "Macrophage",
    "1" = "Macrophage",
    "2" = "Macrophage",
    "3" = "Macrophage",
    "4" = "Macrophage",
    "5" = "Proliferation"
  ),
  "Original cluster 15"
)

# Original cluster 18: resolution 1.
Idents(obj) <- "seurat_clusters"
branch_18 <- subset(obj, idents = "18")
branch_18 <- FindNeighbors(branch_18, dims = 1:20)
branch_18 <- FindClusters(branch_18, resolution = 1)
obj <- write_branch_celltypes(
  obj,
  branch_18,
  c(
    "0" = "Eosinophils",
    "1" = "NE",
    "2" = "Mastocyte",
    "3" = "NE",
    "4" = "Macrophage"
  ),
  "Original cluster 18"
)

# Original NE branch. The first subset is explicitly taken from the Celltype
# identity. FindClusters then establishes the branch-local cluster identity,
# from which original branch cluster 5 is selected and reclustered.
Idents(obj) <- "Celltype"
ne_branch <- subset(obj, idents = "NE")
ne_branch <- FindNeighbors(ne_branch, dims = 1:50)
ne_branch <- FindClusters(ne_branch, resolution = 0.5)
Idents(ne_branch) <- "seurat_clusters"
ne_cluster_5 <- subset(ne_branch, idents = "5")
ne_cluster_5 <- FindNeighbors(ne_cluster_5, dims = 1:50)
ne_cluster_5 <- FindClusters(ne_cluster_5, resolution = 0.5)
obj <- write_branch_celltypes(
  obj,
  ne_cluster_5,
  c(
    "0" = "NE",
    "1" = "Macrophage",
    "2" = "Mastocyte"
  ),
  "NE branch cluster 5"
)

# Apply only the documented legacy-to-final display/interface mapping after all
# original cluster and subcluster assignments are complete.
legacy_to_final <- c(
  "NE" = "Neutrophil",
  "FB" = "Fibroblast",
  "SMC" = "non-OFB-SMC",
  "EC" = "Endothelial cell",
  "CD8_T_Cell" = "CD8 T cell",
  "FB_SMC" = "OFB-SMC",
  "CD4_help_T_Cell" = "CD4 T cell",
  "NK" = "NK Cell",
  "unKnown" = "Pericyte",
  "Macrophage_NE" = "Mixed myeloid cell",
  "Macrophage" = "Macrophage",
  "Proliferation" = "ProMacs",
  "TRMs" = "ProMacs",
  "B_Cell" = "B cell",
  "Eosinophils" = "Eosinophil",
  "Mastocyte" = "Mast cell"
)

legacy_celltypes <- as.character(obj$Celltype)
unmapped_labels <- setdiff(unique(legacy_celltypes), names(legacy_to_final))
if (length(unmapped_labels) > 0) {
  stop(
    "The restored seven-sample annotation contains unmapped legacy labels: ",
    paste(unmapped_labels, collapse = ", ")
  )
}
obj$Celltype <- unname(legacy_to_final[legacy_celltypes])

expected_7data_celltypes <- c(
  "Neutrophil",
  "Fibroblast",
  "non-OFB-SMC",
  "Endothelial cell",
  "CD8 T cell",
  "OFB-SMC",
  "CD4 T cell",
  "NK Cell",
  "Pericyte",
  "Mixed myeloid cell",
  "Macrophage",
  "ProMacs",
  "B cell",
  "Eosinophil",
  "Mast cell"
)
prohibited_7data_labels <- c(
  "Dendritic cell",
  "DC",
  "unKnown",
  "Natural killer cell",
  "Proliferating macrophage"
)

if (anyNA(obj$Celltype)) {
  stop("The final seven-sample Celltype annotation contains NA values.")
}
if (any(obj$Celltype == "")) {
  stop("The final seven-sample Celltype annotation contains empty labels.")
}
present_prohibited_labels <- intersect(
  unique(obj$Celltype),
  prohibited_7data_labels
)
if (length(present_prohibited_labels) > 0) {
  stop(
    "The final seven-sample annotation contains prohibited labels: ",
    paste(present_prohibited_labels, collapse = ", ")
  )
}
if (!setequal(unique(obj$Celltype), expected_7data_celltypes)) {
  stop(
    "The final seven-sample Celltype set does not match the documented ",
    "15-class interface. Observed: ",
    paste(sort(unique(obj$Celltype)), collapse = ", ")
  )
}
if (anyDuplicated(colnames(obj))) {
  stop("The final seven-sample object contains duplicated cell barcodes.")
}

obj$Celltype <- factor(
  obj$Celltype,
  levels = expected_7data_celltypes
)
Idents(obj) <- "Celltype"

saveRDS(obj, file = output_rds)

# Figure 2A producer restored from the original SCP plotting block.
DefaultAssay(obj) <- "RNA"
figure_2a <- CellDimPlot(
  srt = obj,
  group.by = "Celltype",
  reduction = "UMAP",
  theme_use = ggplot2::theme_classic,
  theme_args = list(base_size = 16)
)
ggsave(
  filename = figure_2a_file,
  plot = figure_2a,
  width = 10,
  height = 6
)
