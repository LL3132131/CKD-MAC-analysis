# =============================================================================
# Title: 00_prepare_public_monocle2_object.R
# Purpose: Maintainer-only processed-object preparation; does not recalculate trajectory
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAINTAINER-ONLY
# Required inputs: Private original `data/private/sub3_m2.rds`
# Generated outputs: De-identified `data/processed/sub3_m2_public.rds`
# Upstream dependencies: Private original `data/private/sub3_m2.rds`
# Downstream consumers: Maintainer-only processed-object preparation; does not recalculate trajectory
# Configuration keys: DATA_PROCESSED_DIR, TRAJECTORY_RESULTS_DIR, MONOCLE2_PRIVATE_RDS, MONOCLE2_PUBLIC_RDS
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/03_Trajectory/00_prepare_public_monocle2_object.R
# Reproducibility notes: Existing `CREB3L2_activity` is preserved when present.
# If absent, the plotting script joins the generated activity CSV in memory by
# exact CellID; this preparation script does not generate or remap activity.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

# Prepare a de-identified public copy of the final manuscript Monocle2 object.
# This script changes phenotype metadata only. It never recalculates trajectory.

library(monocle)
library(Biobase)

input_path <- PATHS$MONOCLE2_PRIVATE_RDS
output_path <- PATHS$MONOCLE2_PUBLIC_RDS

if (!file.exists(input_path)) {
  stop("Private input not found: ", input_path)
}
if (!dir.exists(dirname(output_path))) {
  stop("Create the output directory before running: ", dirname(output_path))
}
if (file.exists(output_path)) {
  stop("Refusing to overwrite an existing public object: ", output_path)
}

cds <- readRDS(input_path)
if (!inherits(cds, "CellDataSet")) {
  stop("Input is not a Monocle2 CellDataSet.")
}

required_trajectory_slots <- c(
  "reducedDimS",
  "reducedDimK",
  "minSpanningTree"
)
missing_slots <- setdiff(required_trajectory_slots, slotNames(cds))
if (length(missing_slots) > 0) {
  stop(
    "Required trajectory slots are missing: ",
    paste(missing_slots, collapse = ", ")
  )
}

metadata <- Biobase::pData(cds)
if (!all(c("Pseudotime", "State") %in% colnames(metadata))) {
  stop("Pseudotime and/or State are missing from phenotype metadata.")
}

preserved <- list(
  expression_values = Biobase::exprs(cds),
  Pseudotime = metadata$Pseudotime,
  State = metadata$State,
  reducedDimS = slot(cds, "reducedDimS"),
  reducedDimK = slot(cds, "reducedDimK"),
  trajectory_graph = slot(cds, "minSpanningTree")
)

patient_metadata_fields <- c(
  "patient",
  "patient_id",
  "patient_name",
  "Patient",
  "Patient_ID",
  "PatientName",
  "donor",
  "donor_id",
  "Donor",
  "Donor_ID",
  "subject",
  "subject_id",
  "Subject",
  "Subject_ID",
  "individual",
  "individual_id",
  "Individual",
  "Individual_ID",
  "name",
  "Name"
)
removed_metadata_fields <- intersect(
  patient_metadata_fields,
  colnames(metadata)
)
metadata <- metadata[
  ,
  setdiff(colnames(metadata), removed_metadata_fields),
  drop = FALSE
]

deidentify_sample_id <- function(x) {
  original_factor <- is.factor(x)
  x <- as.character(x)
  manuscript_pattern <- "^(Mild|Severe)_[0-9]+(_|$)"
  recognized <- !is.na(x) & grepl(manuscript_pattern, x)

  result <- x
  result[recognized] <- sub(
    "^((Mild|Severe)_[0-9]+)(_.*)?$",
    "\\1",
    x[recognized]
  )

  unrecognized_values <- unique(x[!recognized & !is.na(x)])
  if (length(unrecognized_values) > 0) {
    replacement <- setNames(
      sprintf("Sample_%02d", seq_along(unrecognized_values)),
      unrecognized_values
    )
    result[!recognized & !is.na(x)] <- replacement[
      x[!recognized & !is.na(x)]
    ]
  }

  if (original_factor) {
    result <- factor(result)
  }
  result
}

sample_id_fields <- intersect(
  c(
    "sample_name",
    "new_sample_name",
    "sample_id",
    "Sample_ID",
    "orig.ident",
    "Batch",
    "batch"
  ),
  colnames(metadata)
)
for (field in sample_id_fields) {
  metadata[[field]] <- deidentify_sample_id(metadata[[field]])
}

celltype_candidates <- c(
  "Celltype",
  "celltype",
  "CellType",
  "cell_type"
)
celltype_field <- intersect(celltype_candidates, colnames(metadata))
if (length(celltype_field) == 0) {
  stop("No cell-type metadata field was found.")
}
celltype_field <- celltype_field[[1]]

legacy_ofb_aliases <- c(
  "Fibro-like_SMC",
  paste0("F", "B-SMC"),
  paste0("F", "B_SMC"),
  "OFB_SMC",
  "OFB-SMC"
)
celltype_map <- c(
  "Contractile_SMC" = "Contractile SMC",
  "Contractile SMC" = "Contractile SMC",
  "Transitional_SMC" = "Transitional SMC",
  "Transitional SMC" = "Transitional SMC",
  "Fibromyocyte" = "Fibromyocyte-like SMC",
  "Fibromyocyte-like_SMC" = "Fibromyocyte-like SMC",
  "Fibromyocyte-like SMC" = "Fibromyocyte-like SMC",
  "Macrophage-like_SMC-1" = "Macrophage-like SMC-1",
  "Macrophage-like SMC-1" = "Macrophage-like SMC-1",
  "Macrophage-like_SMC-2" = "Macrophage-like SMC-2",
  "Macrophage-like SMC-2" = "Macrophage-like SMC-2",
  stats::setNames(rep("OFB-SMC", length(legacy_ofb_aliases)), legacy_ofb_aliases)
)
original_celltypes <- as.character(metadata[[celltype_field]])
unknown_celltypes <- setdiff(
  unique(original_celltypes[!is.na(original_celltypes)]),
  names(celltype_map)
)
if (length(unknown_celltypes) > 0) {
  stop(
    "Unmapped cell-type labels: ",
    paste(unknown_celltypes, collapse = ", ")
  )
}

final_celltype_levels <- c(
  "Contractile SMC",
  "Transitional SMC",
  "Fibromyocyte-like SMC",
  "Macrophage-like SMC-1",
  "Macrophage-like SMC-2",
  "OFB-SMC"
)
metadata[[celltype_field]] <- factor(
  unname(celltype_map[original_celltypes]),
  levels = final_celltype_levels
)
if (celltype_field != "Celltype") {
  metadata$Celltype <- metadata[[celltype_field]]
}

Biobase::pData(cds) <- metadata

stopifnot(
  identical(Biobase::exprs(cds), preserved$expression_values),
  identical(Biobase::pData(cds)$Pseudotime, preserved$Pseudotime),
  identical(Biobase::pData(cds)$State, preserved$State),
  identical(slot(cds, "reducedDimS"), preserved$reducedDimS),
  identical(slot(cds, "reducedDimK"), preserved$reducedDimK),
  identical(
    slot(cds, "minSpanningTree"),
    preserved$trajectory_graph
  )
)

saveRDS(cds, output_path)
message(
  "Removed patient metadata fields: ",
  if (length(removed_metadata_fields) > 0) {
    paste(removed_metadata_fields, collapse = ", ")
  } else {
    "none present"
  }
)
message(
  "De-identified sample metadata fields: ",
  if (length(sample_id_fields) > 0) {
    paste(sample_id_fields, collapse = ", ")
  } else {
    "none present"
  }
)
message("Saved public Monocle2 object: ", output_path)
