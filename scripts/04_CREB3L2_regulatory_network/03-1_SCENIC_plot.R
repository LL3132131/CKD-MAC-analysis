# =============================================================================
# Title: 03-1_SCENIC_plot.R
# Purpose: Figure 4
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: all-cell AUCell loom; `7data_umap.rds`
# Generated outputs: RSS, regulon-rank, heatmap and CREB3L2 activity plots
# Upstream dependencies: the single manuscript all-cell SCENIC run; `7data_umap.rds`
# Downstream consumers: Figure 4
# Configuration keys: DATA_PROCESSED_DIR, SCENIC_DB_DIR, SCENIC_TF_LIST, SCENIC_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/04_CREB3L2_regulatory_network/03-1_SCENIC_plot.R
# Reproducibility notes: This script visualizes the all-cell SCENIC result for
# Figure 4D-E. The independent SMC-enriched run supports Figure 4G.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

# Historical environment note: this script used the project R environment.
setwd(configured_path(PATHS, "SCENIC_RESULTS_DIR", "00.allCell", "7scenic"))
#加载分析包
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(Seurat)

#将loom文件读入R，提取数据
#sce_SCENIC <- open_loom("00.allCell.SCE-auc.loom")
sce_SCENIC <- open_loom(paste0(PATHS$SCENIC_ALL_CELL_PREFIX, ".SCE-auc.loom"))
# exprMat <- get_dgem(sce_SCENIC)#从sce_SCENIC文件提取表达矩阵
# exprMat_log <- log2(exprMat+1) # log处理
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
#RSS分析，查看细胞类型特异性转录因子，需要先加载seurat对象，提取metadata信息，并进行分析！默认是点图！
human_data <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
required_metadata <- c("sample_name", "calcification_group", "Celltype")
missing_metadata <- setdiff(required_metadata, colnames(human_data@meta.data))
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
observed_samples <- unique(as.character(human_data$sample_name))
if (!setequal(observed_samples, expected_samples)) {
  stop(
    "7data_umap.rds must contain exactly the seven documented MAC samples. ",
    "Observed: ", paste(sort(observed_samples), collapse = ", ")
  )
}

expected_groups <- c("Mild", "Severe")
observed_groups <- unique(as.character(human_data$calcification_group))
if (!setequal(observed_groups, expected_groups)) {
  stop(
    "calcification_group must contain exactly Mild and Severe. ",
    "Observed: ", paste(sort(observed_groups), collapse = ", ")
  )
}

atlas_celltype_levels <- c(
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
observed_celltypes <- unique(as.character(human_data$Celltype))
prohibited_7data_celltypes <- intersect(
  observed_celltypes,
  c("Dendritic cell", "DC", "unKnown")
)
if (length(prohibited_7data_celltypes) > 0) {
  stop(
    "7data_umap.rds contains prohibited seven-sample Celltype labels: ",
    paste(sort(prohibited_7data_celltypes), collapse = ", ")
  )
}
if (!"Pericyte" %in% observed_celltypes) {
  stop("7data_umap.rds is missing the required Pericyte Celltype.")
}
if (!setequal(observed_celltypes, atlas_celltype_levels)) {
  stop(
    "Celltype must contain exactly the documented 15 seven-sample cell types. ",
    "Missing: ",
    paste(setdiff(atlas_celltype_levels, observed_celltypes), collapse = ", "),
    "; unexpected: ",
    paste(setdiff(observed_celltypes, atlas_celltype_levels), collapse = ", ")
  )
}

# Preserve the original all-cell SCENIC colors by binding each historical
# color value to its final manuscript cell-type label.
celltype_colors <- c(
  "Neutrophil" = "#ff3a0e",
  "Fibroblast" = "#5fb5f3",
  "non-OFB-SMC" = "#486eca",
  "Endothelial cell" = "#1f71b5",
  "CD8 T cell" = "#004d43",
  "OFB-SMC" = "#cd177b",
  "CD4 T cell" = "#8b0254",
  "NK Cell" = "#435B95",
  "Pericyte" = "#0b0baa",
  "Mixed myeloid cell" = "#de545a",
  "Macrophage" = "#e366b6",
  "ProMacs" = "#F7B46D",
  "B cell" = "#8fb0ff",
  "Eosinophil" = "#3C234A",
  "Mast cell" = "#620d14"
)
if (!identical(names(celltype_colors), atlas_celltype_levels)) {
  stop(
    "SCENIC color names do not match the documented seven-sample ",
    "Celltype order."
  )
}

cellinfo <- human_data@meta.data[, c(
  "Celltype",
  "sample_name",
  "calcification_group",
  "nFeature_RNA",
  "nCount_RNA"
), drop = FALSE]
cellinfo$celltype <- cellinfo$Celltype

######计算细胞特异性TF
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)
rssPlot <- 
  plotRSS(rss,cluster_columns = FALSE,order_rows = TRUE,thr=0.1,varName = "cellType",col.low = '#330066',col.mid = '#66CC66',col.high = '#FFCC33')
rssPlot
##热图
library(reshape2)
library(ggheatmap)

rss_data <- rssPlot$plot$data
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
missing_rss_celltypes <- setdiff(atlas_celltype_levels, colnames(rss_data))
unexpected_rss_celltypes <- setdiff(colnames(rss_data), atlas_celltype_levels)
prohibited_rss_celltypes <- intersect(
  colnames(rss_data),
  c("Dendritic cell", "DC", "unKnown")
)
if (length(prohibited_rss_celltypes) > 0) {
  stop(
    "All-cell RSS contains prohibited seven-sample columns: ",
    paste(sort(prohibited_rss_celltypes), collapse = ", ")
  )
}
if (length(missing_rss_celltypes) > 0 || length(unexpected_rss_celltypes) > 0) {
  stop(
    "All-cell RSS columns do not match the documented atlas cell types. ",
    "Missing: ", paste(missing_rss_celltypes, collapse = ", "),
    "; unexpected: ", paste(unexpected_rss_celltypes, collapse = ", ")
  )
}
rss_data <- rss_data[, atlas_celltype_levels, drop = FALSE]
missing_color_labels <- setdiff(colnames(rss_data), names(celltype_colors))
if (length(missing_color_labels) > 0) {
  stop(
    "Missing SCENIC colors for: ",
    paste(missing_color_labels, collapse = ", ")
  )
}
col_ann <- data.frame(
  group = colnames(rss_data),
  row.names = colnames(rss_data)
)
annotation_colors <- celltype_colors[colnames(rss_data)]
if (!identical(names(annotation_colors), colnames(rss_data))) {
  stop("SCENIC annotation colors are not aligned to RSS columns by name.")
}
col <- list(group = annotation_colors)

text_columns <- sample(colnames(rss_data),0)#不显示列名

p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p
############第二个可视化：将转录因子分析结果与seurat对象结合，可视化类似于seurat
seurat_cell_ids <- colnames(human_data)
aucell_cell_ids <- colnames(regulonAUC)

if (anyDuplicated(seurat_cell_ids)) {
  stop("7data_umap.rds contains duplicated Seurat CellID values.")
}
if (anyDuplicated(aucell_cell_ids)) {
  stop("The all-cell AUCell loom contains duplicated CellID values.")
}
if (!setequal(seurat_cell_ids, aucell_cell_ids)) {
  stop(
    "Seurat and all-cell AUCell CellID sets are not identical. ",
    "Seurat-only: ", length(setdiff(seurat_cell_ids, aucell_cell_ids)),
    "; AUCell-only: ", length(setdiff(aucell_cell_ids, seurat_cell_ids))
  )
}

aucell_cell_order <- match(seurat_cell_ids, aucell_cell_ids)
if (anyNA(aucell_cell_order)) {
  stop("Failed to match every Seurat CellID to the all-cell AUCell loom.")
}
next_regulonAUC <- regulonAUC[, aucell_cell_order]
if (!identical(colnames(next_regulonAUC), seurat_cell_ids)) {
  stop("The reordered AUCell matrix does not follow the Seurat cell order.")
}
dim(next_regulonAUC)

library(SingleCellExperiment)  # 适用于 SingleCellExperiment 对象
regulon_AUC <- regulonAUC@NAMES
human_data@meta.data = cbind(human_data@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))

TF_plot <- c("ATF3(+)","ATF3(-)","ATF4(+)","ATF6(+)","ATF6(-)","BACH1(+)","BACH1(-)","CREB3L2(+)","CREB3L2(-)","FOXC2(+)","GATA6(+)","GATA6(-)","MSX2(+)","RUNX2(-)","SOX9(-)","TEAD1(+)","TEAD1(-)")
OFB_SMC_top_tfs <- c("ETV4(+)","DLX3(+)","TWIST1(+)","GATA6(+)","TBX5(+)","FOXL1(+)","FOXC2(+)","SOX8(+)","TEAD1(+)","ERF(+)","HOXC6(+)","FOXA2(+)","HAND2(+)","THRA(+)","HOXA10(+)","ETV1(+)","DLX5(+)","NR2F2(+)","MXI1(+)","PGR(+)","THRB(+)","EMX2(+)","PBX1(+)","FOXC1(+)","RARG(+)","HMGA2(+)","ZNF704(+)","AR(+)","MAFK(+)","SOX4(+)","CREB3L2(+)","HOXC9(+)","IRF6(+)")

p1 <- DotPlot(human_data, features = TF_plot, group.by ='Celltype' )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
p2 <- DotPlot(human_data, features = OFB_SMC_top_tfs, group.by ='Celltype' )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
#################################
library(cowplot)
library(magrittr)

graph1 <- plotRSS_oneSet(rss, setName = "OFB-SMC", n = 10)
###########only positive
keep_rows <- grepl("\\(\\+\\)$", rownames(rss))  # 正则表达式匹配 "(+)" 结尾
rss_positive <- rss[keep_rows, ]                  # 子集化矩阵
head(rss_positive)  # 查看前几行
dim(rss_positive)   # 查看筛选后的行数-mask_dropouts \ #会变为默认参数

rss_positivePlot <- plotRSS(rss_positive,zThreshold = 1,cluster_columns = FALSE,order_rows = TRUE,thr=0.1,varName = "cellType",col.low = '#330066',col.mid = '#66CC66',col.high = '#FFCC33')
rss_positivePlot

regulon_AUC_positive <- regulon_AUC[grepl("\\(\\+\\)$", regulon_AUC)]
human_data_positive <- human_data
human_data_positive@meta.data = cbind(human_data@meta.data ,t(assay(next_regulonAUC[regulon_AUC_positive,])))

graph2 <- plotRSS_oneSet(rss_positive, setName = "OFB-SMC", n = 5)
ggsave(
  filename = configured_path(
    PATHS,
    "SCENIC_RESULTS_DIR",
    "00.allCell",
    "7scenic",
    "allCellrss_Plot.png"
  ),
  plot = rss_positivePlot$plot
)
#
#########美化plotRSS_oneSet
B_rss <- as.data.frame(rss_positive)
target_celltypes <- c(
  "OFB-SMC",
  "non-OFB-SMC",
  "Fibroblast"
)
missing_target_celltypes <- setdiff(target_celltypes, colnames(B_rss))
if (length(missing_target_celltypes) > 0) {
  stop(
    "Required RSS cell-type columns are missing: ",
    paste(missing_target_celltypes, collapse = ", ")
  )
}
celltype <- target_celltypes
rssRanklist <- list()

for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    ggrepel::geom_text_repel(data= data_rank_plot[1:6,],
                            aes(label=TF), color="black", size=3, fontface="italic", 
                            arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                            point.padding = 0.3, segment.color = 'black', 
                            segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[i]] <- p
}

ggsave(
  filename = "allCell-OFB-SMC-rank.png",
  plot = rssRanklist[[1]],
  width = 5,
  height = 6
)
##############################
##热图
library(reshape2)
library(ggheatmap)

rss_data <- rss_positivePlot$plot$data
rss_data<-dcast(rss_data,
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
missing_rss_celltypes <- setdiff(atlas_celltype_levels, colnames(rss_data))
unexpected_rss_celltypes <- setdiff(colnames(rss_data), atlas_celltype_levels)
prohibited_rss_celltypes <- intersect(
  colnames(rss_data),
  c("Dendritic cell", "DC", "unKnown")
)
if (length(prohibited_rss_celltypes) > 0) {
  stop(
    "Positive-regulon RSS contains prohibited seven-sample columns: ",
    paste(sort(prohibited_rss_celltypes), collapse = ", ")
  )
}
if (length(missing_rss_celltypes) > 0 || length(unexpected_rss_celltypes) > 0) {
  stop(
    "Positive-regulon RSS columns do not match the documented atlas cell types. ",
    "Missing: ", paste(missing_rss_celltypes, collapse = ", "),
    "; unexpected: ", paste(unexpected_rss_celltypes, collapse = ", ")
  )
}
rss_data <- rss_data[, atlas_celltype_levels, drop = FALSE]
missing_color_labels <- setdiff(colnames(rss_data), names(celltype_colors))
if (length(missing_color_labels) > 0) {
  stop(
    "Missing SCENIC colors for: ",
    paste(missing_color_labels, collapse = ", ")
  )
}
col_ann <- data.frame(
  group = colnames(rss_data),
  row.names = colnames(rss_data)
)
annotation_colors <- celltype_colors[colnames(rss_data)]
if (!identical(names(annotation_colors), colnames(rss_data))) {
  stop(
    "Positive-regulon SCENIC annotation colors are not aligned to RSS ",
    "columns by name."
  )
}
col <- list(group = annotation_colors)

text_columns <- sample(colnames(rss_data),0)#不显示列名

p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p
ggsave(
  filename = configured_path(
    PATHS,
    "SCENIC_RESULTS_DIR",
    "00.allCell",
    "7scenic",
    "allCell_pheatmap.png"
  ),
  plot = p
)
# ggsave(filename = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure4", "pheat_SCENIC.pdf"), plot = p, dpi = 900)
######
p2 = RidgePlot(human_data_positive, features = 'CREB3L2(+)') 
ggsave(
  filename = configured_path(
    PATHS,
    "SCENIC_RESULTS_DIR",
    "00.allCell",
    "7scenic",
    "allCell_RidgePlot.png"
  ),
  plot = p2
)
