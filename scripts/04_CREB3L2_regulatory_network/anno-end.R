# =============================================================================
# Title: anno-end.R
# Purpose: Figure 6 / target integration
# Manuscript Figure/Table: Figure 6 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: CREB3L2 broadPeak; mm10 annotation packages
# Generated outputs: peak-annotation table and annotation pie plot
# Upstream dependencies: CREB3L2 broadPeak; mm10 annotation packages
# Downstream consumers: Figure 6 / target integration
# Configuration keys: CUTTAG_FASTQ_DIR, CUTTAG_RESULTS_DIR, REFERENCE_DIR, DATA_PROCESSED_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/04_CREB3L2_regulatory_network/anno-end.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPpeakAnno)
library(ChIPseeker)
library(clusterProfiler)
library(viridis)
library(ggplot2)

bedPeaksFile <- configured_path(PATHS, "CUTTAG_RESULTS_DIR", "CREB3L2-P", "CREB3L2-P_peaks.broadPeak")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- readPeakFile(bedPeaksFile, header = FALSE)
keepChr <- !grepl("_", seqlevels(peak))
seqlevels(peak, pruning.mode = "coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(
  peak,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Mm.eg.db"
)
peakAnno_df <- as.data.frame(peakAnno)

cdark_viridis_11 <- viridis(11, option = "inferno", begin = 0.1, end = 0.8)
warm_viridis_11 <- viridis(11, option = "magma", begin = 0.2, end = 0.9)

figure6_output_dir <- configured_path(
  PATHS,
  "RESULTS_DIR",
  "figures",
  "Figure6"
)
dir.create(
  figure6_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(figure6_output_dir)) {
  stop("Failed to create Figure 6 output directory: ", figure6_output_dir)
}

pdf(
  file = file.path(figure6_output_dir, "CREB3L2-AnnoPie.pdf"),
  width = 8,
  height = 5
)
plotAnnoPie(peakAnno, col = warm_viridis_11)
dev.off()

write.table(
  peakAnno_df,
  file = configured_path(PATHS, "CUTTAG_RESULTS_DIR", "CREB3L2-P", "CREB3L2-P_peaks.broadPeak_annotation.xls"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
