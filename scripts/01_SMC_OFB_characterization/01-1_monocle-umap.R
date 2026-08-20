# =============================================================================
# Title: 01-1_monocle-umap.R
# Purpose: SMC refinement; Figure 2
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data_umap.rds`; SMC marker files
# Generated outputs: `smc_subclusters_initial_seurat.rds`; SMC plots/markers
# Upstream dependencies: `7data_umap.rds`; SMC marker files
# Downstream consumers: SMC refinement; Figure 2
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/01_SMC_OFB_characterization/01-1_monocle-umap.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(viridis)
library(Seurat)
library(dplyr)
library(ggplot2)

output_dir <- configured_path(PATHS, "SMC_RESULTS_DIR", "01_monocle")
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
##################Seurat_V4
ob <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
required_metadata <- c("sample_name", "calcification_group", "Celltype")
missing_metadata <- setdiff(required_metadata, colnames(ob@meta.data))
if (length(missing_metadata) > 0) {
  stop(
    "7data_umap.rds is missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}
expected_samples <- c(
  "Mild_1", "Mild_2", "Mild_3", "Mild_4",
  "Severe_1", "Severe_2", "Severe_3"
)
observed_samples <- unique(as.character(ob$sample_name))
if (!setequal(observed_samples, expected_samples)) {
  stop(
    "7data_umap.rds must contain exactly the seven documented samples. ",
    "Observed: ", paste(sort(observed_samples), collapse = ", ")
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
observed_7data_celltypes <- unique(as.character(ob$Celltype))
if (!setequal(observed_7data_celltypes, expected_7data_celltypes)) {
  stop(
    "7data_umap.rds must use the documented seven-sample major-cell ",
    "interface, including Pericyte and excluding Dendritic cell/unKnown. ",
    "Observed: ",
    paste(sort(observed_7data_celltypes), collapse = ", ")
  )
}
prohibited_7data_labels <- c("Dendritic cell", "DC", "unKnown")
present_prohibited_labels <- intersect(
  observed_7data_celltypes,
  prohibited_7data_labels
)
if (length(present_prohibited_labels) > 0) {
  stop(
    "7data_umap.rds contains prohibited seven-sample labels: ",
    paste(present_prohibited_labels, collapse = ", ")
  )
}
Idents(ob) <- "Celltype"
DefaultAssay(ob) <- 'RNA'
m <- FindAllMarkers(ob, logfc.threshold = 1.5, test.use = "wilcox",only.pos=T)
write.csv(m,file = configured_path(PATHS, "SMC_RESULTS_DIR", "marker.csv"))
sub <- subset(ob, idents = c("non-OFB-SMC", "OFB-SMC"))
smc_samples <- unique(as.character(sub$sample_name))
if (!setequal(smc_samples, expected_samples)) {
  stop(
    "The SMC/OFB-SMC subset must retain all seven samples independently. ",
    "Observed: ", paste(sort(smc_samples), collapse = ", ")
  )
}

#####检查亚群内部是否存在显著批次效应
#DimPlot(sub,group.by = 'sample_name',label = TRUE)
#DimPlot(sub,split.by = 'sample_name',label = TRUE)
##############
sub.list <- SplitObject(object = sub, split.by = "sample_name")
for (i in 1:length(x = sub.list)) {
    sub.list[[i]] <- NormalizeData(object = sub.list[[i]], verbose = FALSE)
    sub.list[[i]] <- FindVariableFeatures(object = sub.list[[i]],
        selection.method = "vst", nfeatures = as.numeric(2000), verbose = FALSE)
    sub.list[[i]] <- ScaleData(sub.list[[i]])
    sub.list[[i]] <- RunPCA(sub.list[[i]])
    sub.list[[i]] <- RunUMAP(sub.list[[i]], dims = 1:30)
}

reference.list <- sub.list
sub.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30,reduction='cca')
sub.integrated <- IntegrateData(anchorset = sub.anchors, dims = 1:30)
######################
DefaultAssay(object = sub.integrated) <- "integrated"
sub.integrated <- FindVariableFeatures(sub.integrated)
sub.integrated <- ScaleData(object = sub.integrated, verbose = FALSE)
sub.integrated <- RunPCA(object = sub.integrated, verbose = FALSE)
sub.integrated <- RunUMAP(object = sub.integrated, reduction = "pca",dims = 1:30,min.dist = 0.2)

sub.integrated = FindNeighbors(object = sub.integrated,dims = 1:30)
sub.integrated.Find <- FindClusters(object = sub.integrated,resolution = 0.5)
sub.integrated.Find <- RunUMAP(object = sub.integrated.Find,dims = 1:as.numeric(36), min.dist = 0.5)
DefaultAssay(sub.integrated.Find)="RNA"
DefaultAssay(sub.integrated.Find)="integrated"
DimPlot(sub.integrated.Find,label = TRUE)
DefaultAssay(object = sub.integrated.Find) <- "integrated"
Idents(sub.integrated.Find) <- 'Celltype'
#Idents(sub.integrated.Find) <- 'seurat_clusters'
sub <- subset(sub.integrated.Find, idents = "non-OFB-SMC")
sub <- FindNeighbors(object = sub,dims = 1:50)
sub <- FindClusters(object = sub,resolution = 0.5)
sub <- RunUMAP(object = sub, reduction = "pca", dims = 1:50, min.dist = 0.5)

DimPlot(sub,label = TRUE)

file1 = sub.integrated.Find@meta.data
file2 = sub@meta.data
file2$cell_id = rownames(file2)
file1$cell_id = rownames(file1)

file1 <- file1[, "cell_id", drop = FALSE]
file2 = file2[,c('cell_id','seurat_clusters')]
file = left_join(file1,file2,by = "cell_id")

l = as.character(unique(file$seurat_clusters))
l = as.character(na.omit(l))
file$seurat_clusters = as.character(file$seurat_clusters)
rownames(file) = file$cell_id
file$seurat_clusters[is.na(file$seurat_clusters)] <- 'missing'
sub.integrated.Find$recluster15 = file$seurat_clusters
meta = sub.integrated.Find@meta.data

sub.integrated.Find$Celltype = as.character(sub.integrated.Find$Celltype)
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '0')] = '0'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '1')] = '1'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '2')] = '2'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '3')] = '3'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '4')] = '4'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '5')] = '5'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '6')] = '6'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '7')] = '7'
sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '8')] = '8'
#sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '9')] = '9'
#sub.integrated.Find$Celltype[which(sub.integrated.Find$recluster15 == '10')] = '10'

Idents(sub.integrated.Find) <- 'Celltype'
DimPlot(sub.integrated.Find, label = T)
m <- FindAllMarkers(sub.integrated.Find, only.pos = TRUE)
write.csv(m, file = configured_path(PATHS, "SMC_RESULTS_DIR", "01_monocle", "sub_zhushi_markers.csv"))
###################
cellreport <- c("ID4", "FABP4", "NET1", "RERGL", "AHNAK", "NEAT1", "MYH10", "CARMN", "CCL19", "IGFBP3", "AGT", "TIMP1", "HLA-B", "HLA-A", "APOE", "APOC1", "CFD", "APOD", "SERPINF1", "C7", "FBLN1", "COMP", "CRTAC1", "CLU", "COL3A1", "COL1A2", "COL1A1", "COL6A2", "COL4A2", "LTBP1", "VCAN", "FN1", "BGN", "IGFBP2", "VCAM1", "ACTA2", "MYL9", "TPM2", "CNN1", "MYH11")
cellreport_marker_dir <- file.path(output_dir, "cellreport_Marker")
dir.create(
  cellreport_marker_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(cellreport_marker_dir)) {
  stop(
    "Failed to create output directory: ",
    cellreport_marker_dir
  )
}
for (i in 1:length(cellreport)){
p = FeaturePlot(sub.integrated.Find, features = cellreport[i])
ggsave(plot = p,file = paste0('./cellreport_Marker/',cellreport[i],'.png'),width = 6,height = 5)
}
p <- DotPlot(sub.integrated.Find, group.by = 'Celltype',features = cellreport,cols = c('blue','red'))  + RotatedAxis()
#ggsave(plot = p,file = paste0('./cellreport_Marker/','cellreport_Dotplot','.png'),width = 30,height = 5)

AP_1 <- c("MYH11", "ACTA2", "MYL9", "TAGLN", "CCL19", "CD9", "APOE", "APOC1", "CD36", "NOTCH3", "RGSS", "MCAM", "STEAP4", "COL1A1", "COL1A2", "COL3A1", "DCN", "FBLN1", "FBLN2", "HSPH1", "HSPD1", "HSP90AA1", "DNAJA1", "PECAM1", "VWF", "EGFL7", "ECSCR", "CXCR4", "CCL4", "C1QA", "C1QB", "IL1B", "MSLN", "MT1M", "MT2A", "MT1X", "MKI67", "CDC20", "TOP2A", "CCNB1", "CCND1")
ap1_output_dir <- file.path(output_dir, "AP_1")
dir.create(
  ap1_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(ap1_output_dir)) {
  stop(
    "Failed to create output directory: ",
    ap1_output_dir
  )
}
for (i in 1:length(AP_1)){

p = FeaturePlot(sub.integrated.Find, features = cellreport[i])

ggsave(plot = p,file = paste0('./AP_1/',AP_1[i],'.png'),width = 6,height = 5)
}
p <- DotPlot(sub.integrated.Find, group.by = 'Celltype',features = AP_1,cols = c('blue','red'))  + RotatedAxis()
ggsave(plot = p,file = paste0('./AP_1/','AP_1_Dotplot','.png'),width = 30,height = 5)


#########
ma <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "before_pheatmap.csv"))
DefaultAssay(object = sub.integrated.Find) <- "RNA"
marker_output_dir <- file.path(output_dir, "Marker")
dir.create(
  marker_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(marker_output_dir)) {
  stop(
    "Failed to create output directory: ",
    marker_output_dir
  )
}
for (i in 1:nrow(ma)){

p = FeaturePlot(sub.integrated.Find, features = ma[i,1])
ggsave(plot = p,file = paste0('./Marker/',ma[i,2],'_',ma[i,1],'.png'),width = 6,height = 5)
}
p <- DotPlot(sub.integrated.Find, group.by = 'Celltype',
        features = unique(ma[,1])) + RotatedAxis()
ggsave(p, file = './Dotplot.png',width = 24,height = 8)

Trans <- c('LGALS3','SIRT6','TUBA1B','KRT8','KRT18','SPARC','NOTCH3','BGN','ARID5B')
DotPlot(sub.integrated.Find, group.by = 'Celltype',features = Trans)  + RotatedAxis()
###################
sub.integrated.Find$Celltype[which(sub.integrated.Find$Celltype %in% c('0','6'))] = 'Contractile SMC'
sub.integrated.Find$Celltype[which(sub.integrated.Find$Celltype %in% c('3'))] = 'Fibromyocyte-like SMC'
sub.integrated.Find$Celltype[which(sub.integrated.Find$Celltype %in% c('1','2','4','7'))] = 'Transitional SMC'
sub.integrated.Find$Celltype[which(sub.integrated.Find$Celltype %in% c('OFB-SMC'))] = 'OFB-SMC'
sub.integrated.Find$Celltype[which(sub.integrated.Find$Celltype %in% c('8'))] = 'Macrophage-like SMC-1'
sub.integrated.Find$Celltype[which(sub.integrated.Find$Celltype %in% c('5'))] = 'Macrophage-like SMC-2'
smc_subtype_levels <- c(
  "Contractile SMC",
  "Transitional SMC",
  "Fibromyocyte-like SMC",
  "Macrophage-like SMC-1",
  "Macrophage-like SMC-2",
  "OFB-SMC"
)
sub.integrated.Find$Celltype <- factor(
  sub.integrated.Find$Celltype,
  levels = smc_subtype_levels
)
Idents(sub.integrated.Find) <- "Celltype"
saveRDS(sub.integrated.Find, file = configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_initial_seurat.rds"))

DimPlot(sub.integrated.Find, label = T)
