# =============================================================================
# Title: target_cellrepot_way.R
# Purpose: Figure 6 / supplementary
# Manuscript Figure/Table: Figure 6 and supplementary; Figure 5 is non-code wet-lab validation
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: Final SMC object and utility functions
# Generated outputs: Manuscript target-gene correlation PDFs
# Upstream dependencies: Final SMC object and utility functions
# Downstream consumers: Figure 6 / supplementary
# Configuration keys: DATA_PROCESSED_DIR, GSE234356_DIR, VALIDATION_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/05_Validation/target_cellrepot_way.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(tidyverse)
library(data.table)
library(broom)
library(ggrepel)
# Set seed for reproducibility
set.seed(1)
output_dir <- PATHS$VALIDATION_RESULTS_DIR
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(output_dir)) {
  stop(
    "Failed to create output directory: ",
    output_dir
  )
}
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
  stop(
    "Failed to create output directory: ",
    figure6_output_dir
  )
}
setwd(output_dir)
# Source our own scRNA analysis utils functions 
source("scRNA_processing_utils.R")
##########################################################################################
# This script will calculate correlations for a gene of interest with all of other genes #
# for cell types of interest.                                                            #
##########################################################################################
########
#2020 Yu G. (2020). Gene Ontology Semantic Similarity Analysis Using GOSemSim.
#2010 Yu G. et al. (2010). GOSemSim: an R package for measuring semantic similarity among GO terms and gene products.
# Load seurat object with pericytes, SMCs and Fibroblasts
rpca_smc_fibro_subset_v3 = read_rds(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
if (!"Celltype" %in% colnames(rpca_smc_fibro_subset_v3@meta.data)) {
  stop("smc_subclusters_final_seurat.rds is missing required metadata: Celltype")
}
DefaultAssay(rpca_smc_fibro_subset_v3) = "RNA"
Idents(rpca_smc_fibro_subset_v3) <- "Celltype"
m <- FindMarkers(rpca_smc_fibro_subset_v3, ident.1 = "OFB-SMC", ident.2 = "Contractile SMC",logfc.threshold = 0.5)
# Get sct counts from matrix
rpca_smc_fibrosct_counts = rpca_smc_fibro_subset_v3@assays$RNA@data
class(rpca_smc_fibrosct_counts)
rpca_smc_fibro_metadata = rpca_smc_fibro_subset_v3@meta.data
cell_types = c("Contractile SMC", "OFB-SMC")
cells_subset <- rownames(rpca_smc_fibro_metadata)[
  rpca_smc_fibro_metadata$Celltype %in% cell_types
]
min_cells <- 0.1 * length(cells_subset)
rpca_smc_fibrosct_counts <- rpca_smc_fibrosct_counts[
  rowSums(rpca_smc_fibrosct_counts[, cells_subset] > 0) >= min_cells,
  cells_subset
]

# Get metadata
rpca_smc_fibro_metadata = rpca_smc_fibro_subset_v3@meta.data

# Define cell types and target genes for correlations
cell_types = c("Contractile SMC", "OFB-SMC")
genes_of_interest = c("MYH11", "CNN1", "TAGLN", "LMOD1", "TPM2", "ACTA2",
                      "FN1","BGN","LUM","LTBP2","MGP","LTBP2", 
                      "DCN", "VCAM1","OMD","IBSP","SOX9","ALPL",
                      "RUNX2","BMP2","SPP1","GATA6","NR4A3","SIRT6")


# Calculate correlations for CRTAC1
# For this, we will source our custom functions (calc_gene_cors and plot_cors) 
# from the Utils script
plot_cors2 = function(cors_df, target_gene, target_cells, genes_to_label) {
  ggplot(cors_df, aes(x=gene_index, y=cor_estimate, 
                      color=cor_estimate, label=ifelse(gene %in% genes_to_label, gene, ""))) + 
    geom_point(size=0.9) +
    geom_text_repel(min.segment.length = 0, 
                    max.overlaps = Inf, color="black", 
                    force=1, fontface="italic", size=4, 
                    segment.color = "grey27") +
    ggtitle(paste0(target_gene, " Pearson correlations in ", target_cells)) +
    xlab("Genes") + 
    ylab("Cor coeff") + 
    labs(color="r") +
    theme_bw() +
    scale_colour_viridis_c(option="magma") +
    scale_x_continuous(limits=c(0,10000)) +
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(face="italic", size=20))
  #legend.position = "bottom")
}
rpca_smc_fibro_metadata$prelim_annotations <- rpca_smc_fibro_metadata$Celltype
colnames(rpca_smc_fibro_metadata)  # "prelim_annotations"
genes <- c("MYO1D", "PKD1", "SMAD7", "PLXDC2", "ANXA4", "AIG1", "PGRMC2", "SLC17A5", "TMEM50B")
for (gene in genes){
  smc_fibrochondro_crtac1 = calc_gene_cors(cell_types = cell_types,
                                         exp_matrix = rpca_smc_fibrosct_counts,
                                         metadata_df = rpca_smc_fibro_metadata, 
                                         target_gene = gene, 
                                         cor_method = "pearson")


# Plot gene correlations for CRTAC1
  crtac1_cors_plot = plot_cors2(smc_fibrochondro_crtac1, gene, "Contractile SMC and OFB-SMC", 
                             genes_of_interest) + 
    custom_theme +
    ggtitle(gene) + 
    theme(aspect.ratio = 1.5)
  ggsave(file=paste0(configured_path(PATHS, "RESULTS_DIR", "figures", "Figure6", trailing = TRUE),gene,'_correlations.pdf'),
         plot = crtac1_cors_plot, width = 5, height = 7)
}
