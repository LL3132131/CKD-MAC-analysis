# =============================================================================
# Title: 02_AddModule.R
# Purpose: Figure 2 / supplementary
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data_umap.rds`; three manuscript gene lists
# Generated outputs: Module-score metadata and plots
# Upstream dependencies: `7data_umap.rds`; three manuscript gene lists
# Downstream consumers: Figure 2 / supplementary
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/01_SMC_OFB_characterization/02_AddModule.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(Seurat)
library(scCustomize) # 需要Seurat版本4.3.0
library(ggpubr)

obj = readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
required_metadata <- c("sample_name", "calcification_group", "Celltype")
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
if (length(missing_metadata) > 0) {
  stop(
    "7data_umap.rds is missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}
Idents(obj) <- "Celltype"
#obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
DefaultAssay(obj) <- 'RNA'
output_dir <- configured_path(PATHS, "SMC_RESULTS_DIR", "02_AddModule")
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
Contraction_score <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue2", "AddModuleScore", "Contraction_EHJ.csv"), header = T) #, col.names = "Gene")
calcification_score <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue2", "AddModuleScore", "calficited.csv"), header = T)
Synthesis_score <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue2", "AddModuleScore", "Synthesis.csv"), header = T)

col = inferno(10)##从 inferno 调色板中均匀提取 10 种颜色
co <- c("#4B0C6BFF","#ED6925FF")
sub = subset(obj,idents = c('non-OFB-SMC','OFB-SMC'))

obj <- AddModuleScore(
    object = obj,
    features = Contraction_score,
    ctrl = 100, #默认值是100
    name = 'Contraction_score',
    seed = 1
)
VlnPlot(obj, features = 'Contraction_score1')

p1 <- FeaturePlot(obj, features = 'Contraction_score1', cols = col,pt.size = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),         # 隐藏坐标轴刻度文本
    axis.ticks = element_blank(),        # 隐藏坐标轴刻度线
    panel.border = element_blank(),      # 移除绘图区域边框
    text = element_text(size = 14, face = "bold"),  # 全局字体设置
    legend.position = c(0.98, 0.02),      # 图例位置：右下角
    legend.justification = c(1, 0),       # 图例对齐方式
    # 调整坐标轴标题位置到左下角
  ) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Contraction_score")

meta = obj@meta.data
m = meta[which(meta$Celltype %in% c('non-OFB-SMC','OFB-SMC')),]
m$Celltype <- factor(m$Celltype, levels = c("non-OFB-SMC", "OFB-SMC"))
my_comparisons = list(c('non-OFB-SMC','OFB-SMC'))
p = ggviolin(m, x="Celltype", y="Contraction_score1", fill = "Celltype",
         palette = co,
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(method = "wilcox.test", label = "p.signif",hide.ns = F,comparisons =my_comparisons,test.args = list(alternative = "greater"))
ggsave(p1, file = 'Contraction_score1.png',width = 6,height = 6)
ggsave(p, file = 'Contraction_score_Vlnplot.png')
#ggsave(p1, file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "Contraction_score1.pdf"),width = 6,height = 6,dpi = 900)
#ggsave(p, file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "Contraction_score_Vlnplot.pdf"),dpi = 900)

############
obj <- AddModuleScore(
    object = obj,
    features = calcification_score,
    ctrl = 100, #默认值是100
    name = 'Calcification_score'
)
VlnPlot(obj, features = 'Calcification_score1')

p2 <- FeaturePlot(obj, features = 'Calcification_score1', cols = col,pt.size = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),         # 隐藏坐标轴刻度文本
    axis.ticks = element_blank(),        # 隐藏坐标轴刻度线
    panel.border = element_blank(),      # 移除绘图区域边框
    text = element_text(size = 14, face = "bold"),  # 全局字体设置
    legend.position = c(0.98, 0.02),      # 图例位置：右下角
    legend.justification = c(1, 0),       # 图例对齐方式
    # 调整坐标轴标题位置到左下角
    axis.title.x = element_text(hjust = 0.05, vjust = -1, margin = margin(t = -10)),  # X轴标题
    axis.title.y = element_text(hjust = 0.05, vjust = 2, margin = margin(r = -10))   # Y轴标题
  ) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Calcification_score")

meta = obj@meta.data
m = meta[which(meta$Celltype %in% c('non-OFB-SMC','OFB-SMC')),]
m$Celltype <- factor(m$Celltype, levels = c("non-OFB-SMC", "OFB-SMC"))
my_comparisons = list(c('non-OFB-SMC','OFB-SMC'))
p = ggviolin(m, x="Celltype", y="Calcification_score1", fill = "Celltype",
         palette = co,
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(method = "wilcox.test", hide.ns = F,label = "p.signif",comparisons =my_comparisons,test.args = list(alternative = "greater"))

ggsave(p2, file = 'Calcification_score1.png',width = 6,height = 6)
ggsave(p, file = 'Calcification_score_Vlnplot.png')
#ggsave(p2, file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "Calcification_score1.pdf"),width = 6,height = 6,dpi = 900)
#ggsave(p, file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "Calcification_score_Vlnplot.pdf"),width = 6,height = 6,dpi = 900)

######################################################Synthesis_score
obj <- AddModuleScore(
    object = obj,
    features = Synthesis_score,
    ctrl = 100, #默认值是100
    name = 'Synthesis_score'
)
VlnPlot(obj, features = 'Synthesis_score1')

p3 <- FeaturePlot(obj, features = 'Synthesis_score1', cols = col,pt.size = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),         # 隐藏坐标轴刻度文本
    axis.ticks = element_blank(),        # 隐藏坐标轴刻度线
    panel.border = element_blank(),      # 移除绘图区域边框
    text = element_text(size = 14, face = "bold"),  # 全局字体设置
    legend.position = c(0.98, 0.02),      # 图例位置：右下角
    legend.justification = c(1, 0),       # 图例对齐方式
    # 调整坐标轴标题位置到左下角
    axis.title.x = element_text(hjust = 0.05, vjust = -1, margin = margin(t = -10)),  # X轴标题
    axis.title.y = element_text(hjust = 0.05, vjust = 2, margin = margin(r = -10))   # Y轴标题
  ) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Synthesis_score")
ggsave(plot = p3,file = 'Synthesis_score1.png',width = 6,height = 6)

meta = obj@meta.data
m = meta[which(meta$Celltype %in% c('non-OFB-SMC','OFB-SMC')),]
m$Celltype <- factor(m$Celltype, levels = c("non-OFB-SMC", "OFB-SMC"))
my_comparisons = list(c('non-OFB-SMC','OFB-SMC'))
p = ggviolin(m, x="Celltype", y="Synthesis_score1", fill = "Celltype",
         palette = co,
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(method = "wilcox.test", hide.ns = F,label = "p.signif",comparisons =my_comparisons,test.args = list(alternative = "greater"))
#ggsave(plot = p3,file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "Synthesis_score1.pdf"),width = 6,height = 6,dpi = 900)
#ggsave(plot = p,file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "Synthesis_score1_Vlnplot.pdf"),width = 6,height = 6,dpi = 900)
###伪时间的AddModuleScore-sub
