# =============================================================================
# Title: 01-2-2_m3-DEG.R
# Purpose: Trajectory, SMC SCENIC and validation
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `smc_subclusters_initial_seurat.rds`
# Generated outputs: `smc_subclusters_final_seurat.rds`; DEG/GSEA/enrichment plots
# Upstream dependencies: `smc_subclusters_initial_seurat.rds`
# Downstream consumers: Trajectory, SMC SCENIC and validation
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/01_SMC_OFB_characterization/01-2-2_m3-DEG.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(SCP)
library(org.Hs.eg.db)
library(DelayedMatrixStats)
library(scop)

output_dir <- configured_path(PATHS, "SMC_RESULTS_DIR", "03_DEG")
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
setwd(output_dir)
obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_initial_seurat.rds"))
if (!"Celltype" %in% colnames(obj@meta.data)) {
  stop("smc_subclusters_initial_seurat.rds is missing required metadata: Celltype")
}
smc_subtype_levels <- c(
  "Contractile SMC",
  "Transitional SMC",
  "Fibromyocyte-like SMC",
  "Macrophage-like SMC-1",
  "Macrophage-like SMC-2",
  "OFB-SMC"
)
Idents(obj) <- "Celltype"
#######################################################
sub <- subset(obj, idents = c('Macrophage-like SMC-2'))
sub <- FindNeighbors(sub, dims = 1:30)
sub <- FindClusters(sub, resolution = 0.2)
DimPlot(sub,label = T)

file1 = obj@meta.data
file2 = sub@meta.data
file2$cell_id = rownames(file2)
file1$cell_id = rownames(file1)

file1 = file1[,c('cell_id','Celltype')]
file2 = file2[,c('cell_id','seurat_clusters')]
file = left_join(file1,file2,by = "cell_id")

l = as.character(unique(file$seurat_clusters))
l = as.character(na.omit(l))
file$seurat_clusters = as.character(file$seurat_clusters)
rownames(file) = file$cell_id
file$seurat_clusters[is.na(file$seurat_clusters)] <- 'missing'
obj$recluster = file$seurat_clusters
meta = obj@meta.data

obj$Celltype = as.character(obj$Celltype)
obj$Celltype[which(obj$recluster == '0')] = 'Macrophage-like SMC-2'
obj$Celltype[which(obj$recluster == '1')] = 'Macrophage-like SMC-2'
obj$Celltype[which(obj$recluster == '2')] = 'Contractile SMC'

obj$Celltype <- factor(obj$Celltype, levels = smc_subtype_levels)
cols = c('#B2DF8A','#1F78B4','#FDBF6F','#FF7F00','#A6CEE3','#33A02C')
Idents(obj) <- "Celltype"
saveRDS(obj, file = configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
if (!identical(levels(obj$Celltype), smc_subtype_levels)) {
  stop("Final SMC metadata do not match the documented six-subtype order.")
}
if (!identical(as.character(Idents(obj)), as.character(obj$Celltype))) {
  stop("Final SMC identities do not match the documented Celltype metadata.")
}
#######################################################
DimPlot(obj, label = T)
DefaultAssay(obj) <- "RNA"

p <- CellDimPlot(srt = obj, group.by = "Celltype",reduction = "UMAP",theme_use = ggplot2::theme_classic, theme_args = list(base_size = 16),palcolor =cols)
ggsave(p,file = configured_path(PATHS, "SMC_RESULTS_DIR", "03_DEG", "sub_SCP_umap.png"),width = 6,height = 6,dpi = 900)
#ggsave(p,file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "sub_SCP_umap.pdf"),width = 6,height = 6,dpi = 900)

obj <- RunDEtest(srt = obj, group_by = "Celltype", fc.threshold = 1, only.pos = FALSE)
p <- VolcanoPlot(srt = obj, group_by = "Celltype", nlabel = 20)
#ggsave(p, file = 'sub-RunDEtest.png',width = 20, height = 15)
DEGs <- obj@tools$DEtest_Celltype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins差异表达分析
obj <- AnnotateFeatures(obj, species = "Homo_sapiens", db = c("Chromosome", "GeneType", "Enzyme", "VerSeDa"))
ht <- FeatureHeatmap(
  srt = obj, group.by = "Celltype", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Homo_sapiens", db = c("GO_BP", "KEGG"), anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4,
  border = FALSE, terms_fontsize = 10, GO_simplify = TRUE, topWord = 10, group_palcolor = cols, feature_split_palcolor = cols, cell_annotation_palcolor = cols
)
print(ht$plot)
p <- ht$plot

ggsave(p, file = 'sub_FeatureHeatmap.png',width = 25 , height = 20,limitsize = FALSE,dpi = 900)
###############富集分析
obj <- RunEnrichment(
  srt = obj, group_by = "Celltype", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)
p <- EnrichmentPlot(
  srt = obj, group_by = "Celltype", group_use = c("Contractile SMC","Transitional SMC","Macrophage-like SMC-1","Macrophage-like SMC-2","Fibromyocyte-like SMC","OFB-SMC"),
  plot_type = "bar"
)
ggsave(p, file = 'sub_EnrichmentPlot.png',width = 28 , height = 18)

EnrichmentPlot(
  srt = obj, group_by = "Celltype", group_use = c("OFB-SMC"),
  plot_type = "wordcloud",word_type = "feature"
)
p <- EnrichmentPlot(
  srt = obj, group_by = "Celltype", group_use = "OFB-SMC",
  plot_type = "network"
)
ggsave(p, file = 'sub_network.png',width = 28 , height = 18)

p <- EnrichmentPlot(srt = obj, group_by = "Celltype", plot_type = "comparison")
ggsave(p, file = 'sub_Dotplot.png',width = 20 , height = 10)
#############################GSEA
obj <- RunGSEA(
  srt = obj, group_by = "Celltype", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)
p1 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0001503")
p2 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0048741")
p3 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0006954")
p4 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0006979")
p5 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0034976")
p6 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0002062")
p7 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0032330")
p8 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0001501")
p9 <- GSEAPlot(srt = obj, group_by = "Celltype", group_use = "OFB-SMC", id_use = "GO:0048706")

#GO:0048741（血管平滑肌细胞分化）,GO:0001503（骨化), GO:0006954（炎症反应）, GO:0006979（氧化应激反应）, GO:0034976(内质网应激),GO:0002062：软骨细胞分化,GO:0032330：软骨细胞分化的正向调控,GO:0001501：骨骼系统发育,GO:0048706：胚胎骨骼系统发育
ggsave(p1, file = 'skeletal muscle cell differentiation.png',width = 6,height = 8)
ggsave(p2, file = 'ossification.png',width = 6,height = 8)
ggsave(p3, file = 'inflammatory response.png',width = 6,height = 8)
ggsave(p4, file = 'response to oxidative stress.png',width = 6,height = 8)
ggsave(p5, file = 'endoplasmic reticulum stress response.png',width = 6,height = 8)
#########################
p <- GSEAPlot(
  srt = obj, group_by = "Celltype", group_use = "OFB-SMC", plot_type = "bar",
  direction = "both", topTerm = 40
)
