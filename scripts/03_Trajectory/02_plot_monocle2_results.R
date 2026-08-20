# =============================================================================
# Title: 02_plot_monocle2_results.R
# Purpose: Default Figure 4 and supplementary plotting entry
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `MONOCLE2_PUBLIC_RDS`; `CREB3L2_ACTIVITY_CSV` when the CDS
# does not already contain `pData$CREB3L2_activity`
# Generated outputs: Trajectory by SMC subtype/Pseudotime/State; CREB3L2, MYO1D and PLXDC2 pseudotime plots; CREB3L2 regulon-activity plot
# Upstream dependencies: `sub3_m2_public.rds` and the SMC-enriched subset
# SCENIC activity extraction chain
# Downstream consumers: Default Figure 4 and supplementary plotting entry
# Configuration keys: TRAJECTORY_RESULTS_DIR, MONOCLE2_PUBLIC_RDS,
# CREB3L2_ACTIVITY_CSV
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/03_Trajectory/02_plot_monocle2_results.R
# Reproducibility notes: This script does not calculate AUCell or alter cell
# barcodes. Activity is matched by exact CellID and added in memory only.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

# Default public Monocle2 plotting entry point.
# Reads the released CDS and never recalculates trajectory, State or Pseudotime.

library(monocle)
library(Biobase)
library(ggplot2)
library(viridis)

input_path <- PATHS$MONOCLE2_PUBLIC_RDS
activity_path <- PATHS$CREB3L2_ACTIVITY_CSV
output_dir <- PATHS$TRAJECTORY_RESULTS_DIR

cds <- readRDS(input_path)
if (!inherits(cds, "CellDataSet")) {
  stop("Input is not a Monocle2 CellDataSet.")
}
if (!all(c("Celltype", "Pseudotime", "State") %in% colnames(pData(cds)))) {
  stop("Celltype, Pseudotime and/or State are missing from the public CDS.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

final_celltype_levels <- c(
  "Contractile SMC",
  "Transitional SMC",
  "Fibromyocyte-like SMC",
  "Macrophage-like SMC-1",
  "Macrophage-like SMC-2",
  "OFB-SMC"
)
pData(cds)$Celltype <- factor(
  as.character(pData(cds)$Celltype),
  levels = final_celltype_levels
)

smc_colors <- c(
  "#B2DF8A",
  "#1F78B4",
  "#FDBF6F",
  "#FF7F00",
  "#A6CEE3",
  "#33A02C"
)

trajectory_subtype <- plot_cell_trajectory(
  cds,
  color_by = "Celltype",
  cell_size = 1
) +
  scale_colour_manual(values = smc_colors, drop = FALSE)
trajectory_pseudotime <- plot_cell_trajectory(
  cds,
  color_by = "Pseudotime",
  cell_size = 1
) +
  viridis::scale_color_viridis(option = "viridis")
trajectory_state <- plot_cell_trajectory(
  cds,
  color_by = "State",
  cell_size = 1
)

ggsave(
  file.path(output_dir, "trajectory_by_SMC_subtype.pdf"),
  trajectory_subtype,
  width = 6,
  height = 6,
  dpi = 900
)
ggsave(
  file.path(output_dir, "trajectory_by_Pseudotime.pdf"),
  trajectory_pseudotime,
  width = 6,
  height = 6,
  dpi = 900
)
ggsave(
  file.path(output_dir, "trajectory_by_State.pdf"),
  trajectory_state,
  width = 6,
  height = 6,
  dpi = 900
)

plot_gene_in_pseudotime <- function(cds, gene, output_file, cell_size = NULL) {
  if (!gene %in% rownames(cds)) {
    stop("Gene is absent from the CDS: ", gene)
  }
  gene_cds <- cds[gene, ]
  if (is.null(cell_size)) {
    gene_plot <- plot_genes_in_pseudotime(
      gene_cds,
      color_by = "Celltype"
    )
  } else {
    gene_plot <- plot_genes_in_pseudotime(
      gene_cds,
      color_by = "Celltype",
      cell_size = cell_size
    ) +
      labs(title = gene) +
      theme(
        strip.text = element_blank(),
        plot.title = element_text(
          face = "bold",
          size = 14,
          hjust = 0.5
        ),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
  ggsave(
    file.path(output_dir, output_file),
    gene_plot,
    width = 7,
    height = 5,
    dpi = 900
  )
}

plot_gene_in_pseudotime(
  cds,
  "CREB3L2",
  "CREB3L2_expression_along_pseudotime.pdf"
)
plot_gene_in_pseudotime(
  cds,
  "MYO1D",
  "MYO1D_expression_along_pseudotime.pdf",
  cell_size = 0.75
)
plot_gene_in_pseudotime(
  cds,
  "PLXDC2",
  "PLXDC2_expression_along_pseudotime.pdf",
  cell_size = 0.75
)

if ("CREB3L2_activity" %in% colnames(pData(cds))) {
  if (anyNA(pData(cds)$CREB3L2_activity)) {
    stop("Existing pData$CREB3L2_activity contains missing values.")
  }
} else {
  if (!file.exists(activity_path)) {
    stop("CREB3L2 activity file not found: ", activity_path)
  }
  activity_table <- read.csv(
    activity_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_activity_columns <- c("CellID", "CREB3L2_activity")
  missing_activity_columns <- setdiff(
    required_activity_columns,
    colnames(activity_table)
  )
  if (length(missing_activity_columns) > 0) {
    stop(
      "CREB3L2 activity file is missing columns: ",
      paste(missing_activity_columns, collapse = ", ")
    )
  }
  activity_cell_ids <- as.character(activity_table$CellID)
  cds_cell_ids <- rownames(pData(cds))
  if (anyNA(activity_cell_ids) || any(!nzchar(activity_cell_ids))) {
    stop("CREB3L2 activity file contains missing or empty CellID values.")
  }
  if (anyDuplicated(activity_cell_ids)) {
    stop("CREB3L2 activity file contains duplicated CellID values.")
  }
  if (anyDuplicated(cds_cell_ids)) {
    stop("Monocle2 pData contains duplicated cell row names.")
  }
  missing_from_activity <- setdiff(cds_cell_ids, activity_cell_ids)
  unmatched_activity <- setdiff(activity_cell_ids, cds_cell_ids)
  if (length(missing_from_activity) > 0 || length(unmatched_activity) > 0) {
    stop(
      "Exact CellID matching failed: ",
      length(missing_from_activity),
      " Monocle2 cells lack activity values; ",
      length(unmatched_activity),
      " activity rows lack Monocle2 cells."
    )
  }
  activity_match <- match(cds_cell_ids, activity_cell_ids)
  if (anyNA(activity_match)) {
    stop("Exact CellID matching produced missing indices.")
  }
  matched_activity <- activity_table$CREB3L2_activity[activity_match]
  if (anyNA(matched_activity)) {
    stop("Matched CREB3L2 activity values contain missing values.")
  }
  pData(cds)$CREB3L2_activity <- matched_activity
}
creb3l2_activity_matrix <- matrix(
  pData(cds)$CREB3L2_activity,
  nrow = 1,
  dimnames = list("CREB3L2", rownames(pData(cds)))
)
activity_feature_data <- new(
  "AnnotatedDataFrame",
  data = data.frame(
    gene_short_name = "CREB3L2",
    row.names = "CREB3L2"
  )
)
activity_cds <- newCellDataSet(
  creb3l2_activity_matrix,
  phenoData = phenoData(cds),
  featureData = activity_feature_data,
  expressionFamily = uninormal()
)
creb3l2_activity_plot <- plot_genes_in_pseudotime(
  activity_cds,
  color_by = "Celltype",
  cell_size = 1.5,
  trend_formula = "~sm.ns(Pseudotime, df = 3)"
) +
  labs(y = "CREB3L2 Activity (AUC)")
ggsave(
  file.path(
    output_dir,
    "CREB3L2_regulon_activity_along_pseudotime.pdf"
  ),
  creb3l2_activity_plot,
  width = 7,
  height = 5,
  dpi = 900
)
