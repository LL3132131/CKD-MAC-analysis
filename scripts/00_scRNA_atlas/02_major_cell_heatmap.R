# =============================================================================
# Title: Major-cell marker heatmap
# Purpose: Generate the Figure 1C heatmap and its source matrix using the final
# 15 manuscript cell types.
# Manuscript Figure/Table: Figure 1C
# Release role: FIGURE SOURCE AND PANEL GENERATION
# Required inputs: DATA_PROCESSED_DIR/8data_cca_3000.rds; manuscript marker and
# gene-order tables under DATA_EXTERNAL_DIR
# Generated outputs: all-cell marker table; Figure 1C PDF/PNG; source matrix;
# supplementary marker heatmap
# Upstream dependencies: 01_umap.R
# Downstream consumers: Figure 1C assembly
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, ATLAS_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/00_scRNA_atlas/02_major_cell_heatmap.R
# Reproducibility notes: FindAllMarkers thresholds, marker lists and heatmap
# scaling are unchanged. Column selection now uses explicit final cell-type names.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)
library(scop)
library(pheatmap)
library(grid)

output_dir <- configured_path(PATHS, "ATLAS_RESULTS_DIR", "03_pheatmap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
Idents(obj) <- "Celltype"
difgene <- FindAllMarkers(
  obj,
  logfc.threshold = 1.5,
  test.use = "wilcox",
  only.pos = TRUE
)
write.csv(
  difgene,
  file.path(output_dir, "all-gene-2.csv"),
  quote = FALSE
)

co <- c(
  "#3361A5", "#3164AB", "#3067B1", "#2F6AB7", "#2E6EBE", "#2C71C4",
  "#2B74CA", "#2A78D1", "#297BD7", "#287EDD", "#2682E4", "#2585EA",
  "#2488F0", "#238CF3", "#218FF4", "#2092F5", "#1F96F6", "#1E99F7",
  "#1C9CF8", "#1B9FF9", "#1AA3FA", "#18A6FB", "#17A9FC", "#16ADFD",
  "#14B0FE", "#16B3FE", "#1FB5FD", "#29B7FC", "#32BAFA", "#3BBCF9",
  "#45BEF8", "#4EC0F6", "#57C2F5", "#61C5F4", "#6AC7F3", "#74C9F1",
  "#7DCBF0", "#86CDEF", "#8CCEED", "#90CFEC", "#95CFEA", "#99D0E9",
  "#9ED0E7", "#A3D1E5", "#A7D1E4", "#ACD2E2", "#B0D3E1", "#B5D3DF",
  "#BAD4DE", "#BED4DC", "#C2D4D9", "#C5D4D3", "#C9D4CE", "#CCD4C8",
  "#CFD4C2", "#D3D4BD", "#D6D3B7", "#D9D3B2", "#DDD3AC", "#E0D3A7",
  "#E3D3A1", "#E7D39B", "#EAD295", "#EBD08B", "#EDCD81", "#EECA77",
  "#F0C86D", "#F1C563", "#F3C359", "#F4C04F", "#F6BD44", "#F8BB3A",
  "#F9B830", "#FBB626", "#FCB31C", "#FBAA1A", "#F99F1C", "#F7941D",
  "#F5891E", "#F37E20", "#F17321", "#EF6822", "#ED5D24", "#EB5225",
  "#E94726", "#E73B27", "#E53029", "#E22929", "#DC2828", "#D72727",
  "#D22626", "#CD2525", "#C72424", "#C22323", "#BD2222", "#B82121",
  "#B22020", "#AD1F1F", "#A81E1E", "#A31D1D", "#800026", "#800026",
  "#800026", "#800026", "#800026", "#800026", "#800026"
)

cluster.averages <- AverageExpression(obj)
average_data <- as.data.frame(cluster.averages$RNA, check.names = FALSE)
missing_columns <- setdiff(final_celltype_levels, colnames(average_data))
if (length(missing_columns) > 0) {
  stop(
    "AverageExpression is missing final cell-type columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
average_data <- average_data[, final_celltype_levels, drop = FALSE]

marker <- read.csv(configured_path(
  PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue1",
  "pheatmap", "pheatmap-marker.csv"
), header = FALSE)
data <- data.frame(V1 = rownames(average_data), average_data, check.names = FALSE)
data <- merge(marker, data, by = "V1", sort = FALSE)
rownames(data) <- data$V1
data <- data[, final_celltype_levels, drop = FALSE]

or2 <- read.csv(configured_path(
  PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue1",
  "pheatmap", "order2.csv"
))
da <- data[or2$X, final_celltype_levels, drop = FALSE]
write.csv(
  da,
  file.path(output_dir, "Figure1C_source_matrix.csv"),
  quote = FALSE
)
p <- pheatmap(
  da,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  color = colorRampPalette(c("MediumTurquoise", "white", "OrangeRed"))(100),
  angle_col = 90,
  display_numbers = FALSE,
  fontsize_col = 10,
  gaps_col = NULL,
  fontfamily_col = "Arial",
  fontface_col = "bold",
  fontfamily_row = "Arial",
  fontface_row = "bold"
)
ggsave(
  filename = file.path(output_dir, "Celltype_pheatmap.pdf"),
  plot = p,
  width = 5.8,
  height = 20.8,
  dpi = 900
)
ggsave(
  filename = file.path(output_dir, "Celltype_pheatmap.png"),
  plot = p,
  width = 5.8,
  height = 20.8,
  dpi = 300
)

ht <- GroupHeatmap(
  srt = obj,
  features = c(
    "FLT1", "VWF", "CSF3",
    "CCL3", "CTSB", "C1QC",
    "S100A8", "S100A9", "CXCR2",
    "CD8A", "TRGC2", "CD3E",
    "FBLN1", "APOD", "DCN",
    "OGN", "FN1", "POSTN",
    "MS4A2", "ADCYAP1", "CTSG",
    "DPYD", "LRMDA", "GCH1",
    "CD1C", "GPR183", "HLA-DPA1",
    "MYH11", "TAGLN", "ACTA2",
    "TCF7", "LEF1",
    "GNLY", "IL2RB", "GZMA",
    "MS4A1", "CD79A", "IGHM",
    "CENPF", "MKI67", "TOP2A",
    "SMPD3", "RNASE2", "ADGRE1"
  ),
  group.by = "Celltype",
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  limits = c(-3, 3),
  ht_params = list(
    heatmap_legend_param = list(
      labels_gp = gpar(fontface = "bold", fontsize = 20),
      title_gp = gpar(fontface = "bold", fontsize = 22)
    )
  )
)
ggsave(
  filename = file.path(output_dir, "GroupHeatmap.pdf"),
  plot = ht$plot,
  width = 5.8,
  height = 10.8,
  dpi = 300
)
