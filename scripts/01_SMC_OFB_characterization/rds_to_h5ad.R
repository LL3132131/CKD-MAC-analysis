# =============================================================================
# Title: rds_to_h5ad.R
# Purpose: scDRS
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data_umap.rds`
# Generated outputs: `7data.h5ad`
# Upstream dependencies: `7data_umap.rds`
# Downstream consumers: scDRS
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/01_SMC_OFB_characterization/rds_to_h5ad.R
# Reproducibility notes: scDRS receives normalized RNA expression. The integrated
# assay is retained for atlas integration, dimensional reduction and clustering,
# but is not exported as the disease-scoring expression matrix.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(sceasy)
library(reticulate)
library(Seurat)

scdrs_dir <- configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs")
dir.create(scdrs_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(scdrs_dir)) {
  stop("Unable to create the scDRS output directory: ", scdrs_dir)
}
setwd(scdrs_dir)
python_bin <- Sys.which(PATHS$PYTHON_BIN)
if (!nzchar(python_bin)) {
  stop("Configured PYTHON_BIN is not available on PATH: ", PATHS$PYTHON_BIN)
}
use_python(python_bin)

data7 <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))

if (!"RNA" %in% names(data7@assays)) {
  stop("7data_umap.rds is missing the required RNA assay.")
}

required_metadata <- c(
  "Celltype",
  "sample_name",
  "calcification_group"
)
missing_metadata <- setdiff(required_metadata, colnames(data7@meta.data))
if (length(missing_metadata) > 0L) {
  stop(
    "7data_umap.rds is missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

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
prohibited_7data_labels <- c("Dendritic cell", "DC", "unKnown")
celltype_values <- as.character(data7$Celltype)
if (anyNA(celltype_values) || any(!nzchar(celltype_values))) {
  stop("7data_umap.rds contains missing or empty Celltype values.")
}
present_prohibited_labels <- intersect(
  unique(celltype_values),
  prohibited_7data_labels
)
if (length(present_prohibited_labels) > 0L) {
  stop(
    "7data_umap.rds contains prohibited seven-sample Celltype labels: ",
    paste(present_prohibited_labels, collapse = ", ")
  )
}
missing_celltypes <- setdiff(expected_7data_celltypes, unique(celltype_values))
unexpected_celltypes <- setdiff(unique(celltype_values), expected_7data_celltypes)
if (length(missing_celltypes) > 0L || length(unexpected_celltypes) > 0L) {
  details <- c(
    if (length(missing_celltypes) > 0L) {
      paste0("missing: ", paste(missing_celltypes, collapse = ", "))
    },
    if (length(unexpected_celltypes) > 0L) {
      paste0("unexpected: ", paste(unexpected_celltypes, collapse = ", "))
    }
  )
  stop(
    "7data_umap.rds Celltype values do not match the final seven-sample ",
    "interface (", paste(details, collapse = "; "), ")."
  )
}
if (!"Pericyte" %in% celltype_values) {
  stop("7data_umap.rds is missing the required Pericyte Celltype.")
}

# scDRS uses normalized RNA expression. The integrated assay is used only for
# atlas integration, dimensional reduction and clustering.
DefaultAssay(data7) <- "RNA"
rna_data <- GetAssayData(
  data7,
  assay = "RNA",
  slot = "data"
)
if (nrow(rna_data) == 0L || ncol(rna_data) == 0L) {
  stop("The normalized RNA data matrix is empty.")
}
if (ncol(rna_data) != ncol(data7)) {
  stop(
    "The normalized RNA data matrix cell count does not match the Seurat object."
  )
}
if (
  is.null(rownames(rna_data)) ||
  length(rownames(rna_data)) == 0L ||
  any(!nzchar(rownames(rna_data)))
) {
  stop("The normalized RNA data matrix is missing gene names.")
}
if (
  is.null(colnames(rna_data)) ||
  length(colnames(rna_data)) == 0L ||
  any(!nzchar(colnames(rna_data)))
) {
  stop("The normalized RNA data matrix is missing cell barcodes.")
}
if (!identical(colnames(rna_data), colnames(data7))) {
  stop(
    "The normalized RNA data matrix cell order does not match the Seurat object."
  )
}

sceasy::convertFormat(
  data7,
  from = "seurat",
  to = "anndata",
  assay = "RNA",
  main_layer = "data",
  outFile = configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs", "7data.h5ad")
)
