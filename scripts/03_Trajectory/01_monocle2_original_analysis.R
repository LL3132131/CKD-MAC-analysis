# =============================================================================
# Title: 01_monocle2_original_analysis.R
# Purpose: Provenance/optional only; not used by the default public plotting chain
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: PROVENANCE / OPTIONAL
# Required inputs: EXTERNAL_INPUT original processed `DATA_PROCESSED_DIR/sub3_m3.rds`
# Generated outputs: In-memory original Monocle2 CDS construction
# Upstream dependencies: EXTERNAL_INPUT provenance object; not used by default plotting
# Downstream consumers: Provenance/optional only; not used by the default public plotting chain
# Configuration keys: DATA_PROCESSED_DIR, TRAJECTORY_RESULTS_DIR, MONOCLE2_PRIVATE_RDS, MONOCLE2_PUBLIC_RDS
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/03_Trajectory/01_monocle2_original_analysis.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

# Provenance only: original Monocle2 CDS construction used for the manuscript.
# This script is not the default public plotting entry point.

library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle)

input_seurat <- configured_path(PATHS, "DATA_PROCESSED_DIR", "sub3_m3.rds")
ordering_gene_file <- configured_path(
  PATHS,
  "TRAJECTORY_RESULTS_DIR",
  "Promoter_celltype-fina-Fibro-like_SMC_top_gene_EX_1-way.csv"
)

obj <- readRDS(input_seurat)
SMCgene <- FindAllMarkers(
  obj,
  log2FC.threshold = 0.2,
  return.thresh = 0.05,
  test.use = "wilcox"
)

marker <- list()
id <- "EX"
for (i in id) {
  data <- SMCgene
  for (j in unique(data$cluster)) {
    data1 <- data[data$cluster == j, ]
    data1 <- data1[order(-data1$avg_log2FC, data1$p_val), ]
    if (length(data1[, 1]) > 100) {
      data2 <- data1[1:100, ]
      marker[[paste0(i, "_", j)]] <- data2
    } else {
      marker[[paste0(i, "_", j)]] <- data1
    }
  }
}

m_list <- do.call(rbind, marker)
m_list1 <- as.data.frame(unique(m_list$gene))
colnames(m_list1) <- "top_gene"
write.csv(m_list1, ordering_gene_file, quote = FALSE)
m_list1 <- read.csv(ordering_gene_file)

root <- "Contractile_SMC"

obj <- CreateSeuratObject(
  counts = obj@assays$RNA@counts,
  meta.data = obj@meta.data
)
Mono_matrix <- GetAssayData(obj, slot = "counts")
feature_ann <- data.frame(
  gene_id = rownames(Mono_matrix),
  gene_short_name = rownames(Mono_matrix)
)
rownames(feature_ann) <- rownames(Mono_matrix)
Mono_fd <- new("AnnotatedDataFrame", data = feature_ann)
Mono_pd <- new("AnnotatedDataFrame", data = obj@meta.data)

Mono.cds <- newCellDataSet(
  Mono_matrix,
  phenoData = Mono_pd,
  featureData = Mono_fd,
  expressionFamily = negbinomial.size()
)
gc()
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
disp_table <- dispersionTable(Mono.cds)
Mono.cds <- setOrderingFilter(Mono.cds, m_list1$top_gene)
Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 3,
  method = "DDRTree"
)

# The original script declared root = "Contractile_SMC" but called
# orderCells without a root_state argument. This exact behavior is retained.
Mono.cds <- orderCells(Mono.cds)

colnames(pData(Mono.cds))[
  colnames(pData(Mono.cds)) == "sample_name"
] <- "new_sample_name"

# The manuscript CDS was subsequently stored as sub3_m2.rds.
# This provenance script intentionally does not overwrite that object.
